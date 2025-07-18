#pragma once

#include <iostream>
#include <numeric>

#include "tree_summary.hpp"
#include "opvertex.hpp"
#include <boost/format.hpp>
#include <boost/program_options.hpp>

using namespace std;
using namespace boost;

namespace op {

class OP {
    public:
                            OP();
                            ~OP();

        void                clear();
        void                processCommandLineOptions(int argc, const char * argv[]);
        double              calcBHVDistance(unsigned ref_index, unsigned test_index) const;
        double              calcKFDistance(unsigned ref_index, unsigned test_index) const;
        void                run();

    private:

        static double opCalcTreeIDLength(
            const Split::treeid_t & splits);
        double opCalcLeafContribution(
            const Split::treeid_t & Alvs,
            const Split::treeid_t & Blvs) const;
        double opFindCommonEdges(
            const Split::treeid_t & A,
            const Split::treeid_t & B,
            vector<Split> & common_edges) const;
        void opSplitAtCommonEdges(
            const vector<Split> & common_edges,
            vector<pair<Split::treeid_t,Split::treeid_t> > & in_pairs) const;
#if defined(OP_SAVE_DOT_FILE)
        static void opSaveIncompatibilityGraph(
            vector<OPVertex> & avect,
            vector<OPVertex> & bvect);
#endif
        static void opEdmondsKarp(
            vector<OPVertex> & avect,
            vector<OPVertex> & bvect,
            Split::treeid_t & C1,
            Split::treeid_t & C2,
            Split::treeid_t & D1,
            Split::treeid_t & D2,
            bool quiet);
        bool opRefineSupport(
            const Split::treeid_pair_t & AB,
            Split::treeid_pair_t & AB1,
            Split::treeid_pair_t & AB2) const;
        double opCalcGeodesicDist(
            vector<Split::treeid_pair_t> & ABpairs) const;

        bool                    _quiet;
        string                  _tree_file_name;
        string                  _distance_measure;
        TreeSummary::SharedPtr  _tree_summary;

        static string           _program_name;
        static unsigned        _major_version;
        static unsigned        _minor_version;

    };

inline OP::OP() : _quiet(true) {
    //cout << "Constructing a SStrom" << endl;
    clear();
}

inline OP::~OP() = default;

inline void OP::clear() {
    _quiet = true;
    _tree_file_name = "";
    _distance_measure = "geodesic";
    _tree_summary   = nullptr;
}

inline void OP::processCommandLineOptions(int argc, const char * argv[]) {
    program_options::variables_map       vm;
    program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("treefile,t",  program_options::value(&_tree_file_name)->required(), "name of data file in NEXUS format (required, no default)")
        ("dist", program_options::value(&_distance_measure), "specify either kf or geodesic (default: geodesic)")
        ("quiet,q", program_options::value(&_quiet), "suppress all output except for errors (default: yes)")
        ;
    program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
    try {
        const program_options::parsed_options & parsed = program_options::parse_config_file< char >("op.conf", desc, false);
        program_options::store(parsed, vm);
    }
    catch(program_options::reading_file &) {
        cout << "Note: configuration file (kfdist.conf) not found" << endl;
    }
    program_options::notify(vm);

    // If the user specified --help on the command line, output usage summary and quit
    if (vm.count("help") > 0) {
        cout << desc << "\n";
        exit(1);
    }

    // If the user specified --version on the command line, output the version and quit
    if (vm.count("version") > 0) {
        cout << str(format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << endl;
        exit(1);
    }
}

inline double OP::opCalcTreeIDLength(const Split::treeid_t & splits) {
    double length = 0.0;
    for (auto & split : splits) {
        length += pow(split.getEdgeLen(),2);
    }
    return sqrt(length);
}

inline double OP::opCalcLeafContribution(const Split::treeid_t & Alvs, const Split::treeid_t & Blvs) const {
    if (!_quiet) {
        cout << "Leaves from starting tree:" << endl;
        for (auto & a : Alvs) {
            cout << "  " << a.createPatternRepresentation() << endl;
        }

        cout << "Leaves from ending tree:" << endl;
        for (auto & b : Blvs) {
            cout << "  " << b.createPatternRepresentation() << endl;
        }
    }

    // Compute leaf contribution
    double leaf_contribution_squared = 0.0;
    for (auto & b : Blvs) {
        auto it = find(Alvs.begin(), Alvs.end(), b);
        assert(it != Alvs.end());
        double leafa = it->getEdgeLen();
        double leafb = b.getEdgeLen();
        leaf_contribution_squared += pow(leafa-leafb, 2);
    }

    if (!_quiet)
        cout << str(format("\nLeaf contribution (squared) = %.9f") % leaf_contribution_squared) << endl;
    return leaf_contribution_squared;
}

inline double OP::opFindCommonEdges(const Split::treeid_t & A, const Split::treeid_t & B, vector<Split> & common_edges) const {
    set_intersection(
        A.begin(), A.end(),
        B.begin(), B.end(),
        back_inserter(common_edges)
    );
    double common_edge_contribution_squared = 0.0;
    for (auto & s : common_edges) {
        auto itA = find(A.begin(), A.end(), s);
        double edgeA = itA->getEdgeLen();
        auto itB = find(B.begin(), B.end(), s);
        double edgeB = itB->getEdgeLen();
        common_edge_contribution_squared += pow(edgeA-edgeB, 2);
    }

    if (!_quiet) {
        cout << "\nCommon edges:" << endl;
        for (auto & s : common_edges) {
            cout << "  " << s.createPatternRepresentation() << endl;
        }
        cout << str(format("Common edge contribution (squared) = %.9f") % common_edge_contribution_squared) << endl;
    }
    return common_edge_contribution_squared;
}

inline void OP::opSplitAtCommonEdges(const vector<Split> & common_edges, vector<pair<Split::treeid_t,Split::treeid_t> > & in_pairs) const {
    vector<pair<Split::treeid_t,Split::treeid_t> > out_pairs;
    for (auto & common : common_edges) {
        //cout << "\ncommon: " << common.createPatternRepresentation() << endl;

        // Create a mask that can be used to zero out all bits in common except the first
        Split mask = common;
        mask.invertBits();
        unsigned first_common_bit = common.findFirstSetBit();
        mask.setBitAt(first_common_bit);
        //cout << "  mask: " << mask.createPatternRepresentation() << endl;

        for (auto & inpair : in_pairs) {
            //cout << "\n***** new tree pair *****" << endl;
            // Separate out splits in starting (a_splits) vs. ending (b_splits) (sub)trees
            Split::treeid_t & a_splits = inpair.first;
            Split::treeid_t & b_splits = inpair.second;

            // Create split sets to hold splits subsumed in s vs. other splits
            Split::treeid_t a_common_splits, b_common_splits;
            Split::treeid_t a_other_splits, b_other_splits;

            // Divvy up a_splits to a_common_splits and a_other_splits
            //cout << "  Divvying up a_splits:" << endl;
            for (auto & asplit : a_splits) {
                //cout << "    asplit: " << asplit.createPatternRepresentation();
                bool is_common = (asplit == common);
                if (is_common) {
                    //cout << " (common)" << endl;
                }
                else {
                    if (asplit.subsumedIn(common)) {
                        //cout << " (subsumed in common)" << endl;
                        a_common_splits.insert(asplit);
                    }
                    else {
                        //cout << " (other)" << endl;
                        Split masked = asplit;
                        masked.bitwiseAnd(mask);
                        //cout << "    masked: " << masked.createPatternRepresentation() << endl;
                        a_other_splits.insert(masked);
                    }
                }
            }

            // Divvy up b_splits to b_common_splits and b_other_splits
            //cout << "  Divvying up b_splits:" << endl;
            for (auto & bsplit : b_splits) {
                //cout << "  bsplit: " << bsplit.createPatternRepresentation();
                bool is_common = (bsplit == common);
                if (is_common) {
                    //cout << " (common)" << endl;
                }
                else {
                    if (bsplit.subsumedIn(common)) {
                        //cout << " (subsumed in common)" << endl;
                        b_common_splits.insert(bsplit);
                    }
                    else {
                        //cout << " (other)" << endl;
                        Split masked = bsplit;
                        masked.bitwiseAnd(mask);
                        //cout << "  masked: " << masked.createPatternRepresentation() << endl;
                        b_other_splits.insert(masked);
                    }
                }
            }

            // Create two new tree pairs (a_common_splits, b_common_splits) and (a_other_splits, b_other_splits)
            out_pairs.emplace_back(a_common_splits, b_common_splits);
            out_pairs.emplace_back(a_other_splits, b_other_splits);

            if (!_quiet) {
                cout << "\nSplitting trees at common split: " << common.createPatternRepresentation() << endl;
                cout << "  Left subtree above common split:" << endl;
                for (auto & asplit : a_common_splits) {
                    cout << "    " << asplit.createPatternRepresentation() << endl;
                }
                cout << "  Right subtree above common split:" << endl;
                for (auto & bsplit : b_common_splits) {
                    cout << "    " << bsplit.createPatternRepresentation() << endl;
                }
                cout << "  Left subtree below common split:" << endl;
                for (auto & asplit : a_other_splits) {
                    cout << "    " << asplit.createPatternRepresentation() << endl;
                }
                cout << "  Right subtree below common split:" << endl;
                for (auto & bsplit : b_other_splits) {
                    cout << "    " << bsplit.createPatternRepresentation() << endl;
                }
            }
        }

        // Swap in_pairs and out_pairs
        in_pairs.swap(out_pairs);
        out_pairs.clear();
    }
}

#if defined(OP_SAVE_DOT_FILE)
inline void OP::opSaveIncompatibilityGraph(vector<OPVertex> & avect, vector<OPVertex> & bvect) {
    unsigned minsize = min(avect.size(), bvect.size());
    unsigned maxsize = max(avect.size(), bvect.size());

    // Save all the entities needed
    // tuple key: <0> name, <1>capacity, <2>edgelen, <3>split, <4>shape, <5>color
    typedef vector<tuple<string, string, string, string, string, string> > gnode_t;
    gnode_t anodes, bnodes;
    for (unsigned i = 0; i < avect.size(); i++) {
        string name     = str(format("a%d") % (i+1));
        string capacity = (avect[i]._capacity == 0.0 ? "\"0\"" : str(format("\"%.3f\"") % avect[i]._capacity));
        string edgelen = str(format("\"%.3f\"") % avect[i]._split->getEdgeLen());
        string split   = str(format("\"%.3f\"") % avect[i]._split->createPatternRepresentation(true));
        string shape   = (avect[i]._capacity == 0.0 ? "circle" : "box");
        string color   = (avect[i]._capacity == 0.0 ? "red" : "black");
        anodes.emplace_back(name, capacity, edgelen, split, shape, color);
    }

    for (unsigned i = 0; i < bvect.size(); i++) {
        string name     = str(format("b%d") % (i+1));
        string capacity = (bvect[i]._capacity == 0.0 ? "\"0\"" : str(format("\"%.3f\"") % bvect[i]._capacity));
        string edgelen = str(format("\"%.3f\"") % bvect[i]._split->getEdgeLen());
        string split   = str(format("\"%.3f\"") % bvect[i]._split->createPatternRepresentation(true));
        string shape   = (bvect[i]._capacity == 0.0 ? "circle" : "box");
        string color   = (bvect[i]._capacity == 0.0 ? "red" : "black");
        bnodes.emplace_back(name, capacity, edgelen, split, shape, color);
    }

    ofstream dotf("graph.dot");

    dotf << "digraph G {\n";
    dotf << "\trankdir=LR;\n\n";

    for (unsigned i = 0; i < minsize; i++) {
        dotf << str(format("\t{ rank = %d; %s [shape = plain]; %s [label = %s, shape = %s, color = %s]; %s [label = %s, shape = %s, color = %s]; %s [shape = plain]}\n")
            % (i+1)
            % get<3>(anodes[i])
            % get<0>(anodes[i])
            % get<1>(anodes[i])
            % get<4>(anodes[i])
            % get<5>(anodes[i])
            % get<0>(bnodes[i])
            % get<1>(bnodes[i])
            % get<4>(bnodes[i])
            % get<5>(bnodes[i])
            % get<3>(bnodes[i])
        );
    }

    for (unsigned i = minsize; i < maxsize; i++) {
        if (i < avect.size()) {
            dotf << str(format("\t{ rank = %d; %s [label = %s, shape = %s, color = %s]; }\n")
                % (i+1)
                % get<0>(anodes[i])
                % get<1>(anodes[i])
                % get<4>(anodes[i])
                % get<5>(anodes[i])
            );
        }
        else {
            dotf << str(format("\t{ rank = %d; %s [label = %s, shape = %s, color = %s]; }\n")
                % (i+1)
                % get<0>(bnodes[i])
                % get<1>(bnodes[i])
                % get<4>(bnodes[i])
                % get<5>(bnodes[i])
            );
        }
    }

    // dotf << "\t{\n";

    for (auto & anode : anodes) {
        dotf << str(format("\t\t%s -> %s [style = invis];\n") % get<3>(anode) % get<0>(anode));
    }

    for (unsigned i = 0; i < avect.size(); ++i) {
        auto a = avect[i]._split;
        for (unsigned j = 0; j < bvect.size(); ++j) {
            auto b = bvect[j]._split;
            if (!a->compatibleWith(*b)) {
                dotf << str(format("\t\t%s -> %s;\n") % get<0>(anodes[i]) % get<0>(bnodes[j]));
            }
        }
    }

    for (auto & bnode : bnodes) {
        dotf << str(format("\t\t%s -> %s [style = invis, dir = back];\n") % get<0>(bnode) % get<3>(bnode));
    }

    // dotf << "\t}\n";
    dotf << "}\n";
    dotf.close();
}
#endif

inline void OP::opEdmondsKarp(
        vector<OPVertex> & avect,
        vector<OPVertex> & bvect,
        Split::treeid_t & C1,
        Split::treeid_t & C2,
        Split::treeid_t & D1,
        Split::treeid_t & D2,
        bool quiet) {
    // Assumes avect and bvect form an incompatibility graph that is solvable
    // (i.e., some vertices in avect are compatible with some vertices in bvect)
    if (!quiet) {
        cout << "\nEdmonds-Karp" << endl;

#if defined(OP_SAVE_DOT_FILE)
        // Uncomment the line below to save graph.dot (visualization of incompatibility graph in dot language)
        opSaveIncompatibilityGraph(avect, bvect);
#endif

    }

    double max_flow = 0.0;
    bool done = false;
    while (!done) {
        vector<OPVertex *> route;

        // Make sure none of the "b" vertices are marked as visited
        for (auto & b : bvect) {
            b._visited = false;
        }

        // Add all avect vertices to the route if they have residual capacity > 0
        for (auto & a : avect) {
            if (a._capacity > 0.0) {
                route.push_back(&a);
            }
        }

        // Add all children of the "a" vertices already in the route if they have nonzero capacity
        // and if they haven't already been added
        auto route_size = static_cast<unsigned>(route.size());
        for (unsigned aindex = 0; aindex < route_size; aindex++) {
            OPVertex * a = route[aindex];
            for (auto & b : a->_children) {
                if (b->_capacity > 0.0 && !b->_visited) {
                    b->_visited = true;
                    b->_parent_index = static_cast<int>(aindex);
                    route.push_back(b);
                }
            }
        }

        // Find "b" vertex that first reaches the sink
        OPVertex * last = nullptr;
        for (auto & r : route) {
            if (r->_parent_index > -1) {
                last = r;
                break;
            }
        }

        if (!last)
            // If we did not reach the sink, we're done
            done = true;
        else {
            // Find minimum capacity along the route
            double min_capacity = last->_capacity;
            if (route[last->_parent_index]->_capacity < min_capacity) {
                min_capacity = route[last->_parent_index]->_capacity;
            }
            max_flow += min_capacity;

            if (!quiet) {
                cout << "\nRoute (asterisks show path obtained by following parents from sink to source):" << endl;
                for (unsigned i = 0; i < route.size(); i++) {
                    if (i == last->_parent_index || route[i] == last)
                        cout << str(format("    %s (capacity = %.3f) *") % route[i]->_split->createPatternRepresentation() % route[i]->_capacity) << endl;
                    else
                        cout << str(format("    %s (capacity = %.3f)") % route[i]->_split->createPatternRepresentation() % route[i]->_capacity) << endl;
                }
                cout << "  Min capacity along route: " << min_capacity << endl;
            }

            // Reduce capacities along the route by an amount min_capacity
            last->_capacity -= min_capacity;
            if (fabs(last->_capacity) < 1e-10) {
                last->_capacity = 0.0;
            }

            route[last->_parent_index]->_capacity -= min_capacity;
            if (fabs(route[last->_parent_index]->_capacity) < 1e-10) {
                route[last->_parent_index]->_capacity = 0.0;
            }

#if defined(OP_SAVE_DOT_FILE)
            if (!quiet) {
                // Uncomment the line below to save graph.dot (visualization of incompatibility graph in dot language)
                opSaveIncompatibilityGraph(avect, bvect);
            }
#endif
        }
    }

    // Identify C1, C2
    for (auto & a : avect) {
        if (a._capacity > 0.0) {
            C2.insert(*(a._split));
        }
        else {
            C1.insert(*(a._split));
        }
    }

    // Identify D1, D2
    for (auto & b : bvect) {
        if (b._capacity > 0.0) {
            D1.insert(*(b._split));
        }
        else {
            D2.insert(*(b._split));
        }
    }

    if (!quiet) {
        double C1len = opCalcTreeIDLength(C1);
        double C2len = opCalcTreeIDLength(C2);
        double D1len = opCalcTreeIDLength(D1);
        double D2len = opCalcTreeIDLength(D2);
        // double C1sq = pow(C1len, 2);
        // double D2sq = pow(D2len, 2);
        // double Asq = pow(C1len + C2len, 2);
        // double Bsq = pow(D1len + D2len, 2);
        // double C1_div_D2 = C1sq/Asq + D2sq/Bsq;
        cout << "\nResults:" << endl;
        cout << str(format("  ||C1|| = %.9f") % C1len) << endl;
        cout << str(format("  ||C2|| = %.9f") % C2len) << endl;
        cout << str(format("  ||D1|| = %.9f") % D1len) << endl;
        cout << str(format("  ||D2|| = %.9f") % D2len) << endl;
        cout << "  Check whether ||C1||/||D1|| < ||C2||/||D2||" << endl;
        cout << str(format("    %.9f < %.9f") % (C1len/D1len) % (C2len/D2len));
        cout << endl;
    }
}

inline bool OP::opRefineSupport(const Split::treeid_pair_t & AB, Split::treeid_pair_t & AB1, Split::treeid_pair_t & AB2) const {
    // Create a vector of incompatibility graph vertices
    vector<OPVertex> avect(AB.first.size());
    vector<OPVertex> bvect(AB.second.size());

    unsigned aindex = 0;
    double asum = 0.0;
    for (auto & a : AB.first) {
        asum += pow(a.getEdgeLen(),2);
    }
    for (auto & a : AB.first) {
        avect[aindex]._split = &a;
        avect[aindex]._capacity = pow(a.getEdgeLen(),2)/asum;
        aindex++;
    }

    unsigned bindex = 0;
    double bsum = 0.0;
    for (auto & b : AB.second) {
        bsum += pow(b.getEdgeLen(),2);
    }
    for (auto & b : AB.second) {
        bvect[bindex]._split = &b;
        bvect[bindex]._capacity = pow(b.getEdgeLen(),2)/bsum;
        bindex++;
    }

    // Create the incompatibility graph
    unsigned nincompatible = 0;
    auto asize = static_cast<unsigned>(avect.size());
    auto bsize = static_cast<unsigned>(bvect.size());
    for (unsigned i = 0; i < asize; i++) {
        for (unsigned j = 0; j < bsize; j++) {
            const Split * a = avect[i]._split;
            assert(a);
            const Split * b = bvect[j]._split;
            assert(b);
            if (!a->compatibleWith(*b)) {
                // Create an edge in the incompatibility graph
                avect[i]._children.push_back(&bvect[j]);
                nincompatible++;
            }
        }
    }

    bool success = false;
    if (nincompatible < asize*bsize) {
        // At least one independent pair of vertices exists
        // Carry out Edmonds-Karp algorithm to find min-weight vertex cover (identifies max weight independent set)
        // In Owens-Provan terminology,
        //   C1 = A's contribution to vertex cover     (equals AB1.first)
        //   C2 = A's contribution to independent set  (equals AB2.first)
        //   D1 = B's contribution to independent set  (equals AB1.second)
        //   D2 = B's contribution to vertex cover     (equals AB2.second)
        // C2 and D1 (AB2.first and AB1.second) are compatible sets of splits
        // C1 and D2 (AB1.first and AB2.second) compose the minimum weight vertex cover
        // ||AB1.first||/||Ab1.second|| <= ||AB2.first||/||Ab2.second||
        // Length of this segment of the geodesic is
        //   L = sqrt{ (||AB1.first|| + ||AB1.second||)^2 +  (||AB2.first|| + ||AB2.second||)^2 }
        // Orthants crossed:
        //   start:        AB1.first,  AB1.second
        //   intermediate: AB1.second, AB2.first
        //   finish:       AB1.second, AB2.second
        opEdmondsKarp(avect, bvect, AB1.first, AB2.first, AB1.second, AB2.second, _quiet);
        success = true;

        if (!_quiet) {
            cout << "\nSuccessfully refined support:" << endl;
            // cout << "  Input A vertices:" << endl;
            // for (auto & a : AB.first) {
            //     cout << "    " << a.createPatternRepresentation() << endl;
            // }
            // cout << "  Input B vertices:" << endl;
            // for (auto & b : AB.second) {
            //     cout << "    " << b.createPatternRepresentation() << endl;
            // }
            cout << "  Output A1 vertices:" << endl;
            for (auto & a : AB1.first) {
                cout << "    " << a.createPatternRepresentation(true) << endl;
            }
            cout << "  Output B1 vertices:" << endl;
            for (auto & b : AB1.second) {
                cout << "    " << b.createPatternRepresentation(true) << endl;
            }
            cout << "  Output A2 vertices:" << endl;
            for (auto & a : AB2.first) {
                cout << "    " << a.createPatternRepresentation(true) << endl;
            }
            cout << "  Output B2 vertices:" << endl;
            for (auto & b : AB2.second) {
                cout << "    " << b.createPatternRepresentation(true) << endl;
            }
            cout << endl;
        }
    }
    return success;
}

inline double OP::opCalcGeodesicDist(vector<Split::treeid_pair_t> & ABpairs) const {
    // Assumes a_splits and b_splits have no common edges
    vector<Split::treeid_pair_t> support;
    bool done = false;
    while (!done) {
        unsigned nrefinements = 0;
        for (auto & ABpair : ABpairs) {
            Split::treeid_pair_t AB1;
            Split::treeid_pair_t AB2;
            bool success = opRefineSupport(ABpair, AB1, AB2);
            if (success) {
                // ABpair was successfully refined, so add AB1 and AB2 to support
                support.push_back(AB1);
                support.push_back(AB2);
                nrefinements++;
            }
            else {
                // ABpair was not successfully refined, so add ABpair to support
                support.push_back(ABpair);
            }
        }
        done = (nrefinements == 0);
        if (!done) {
            ABpairs = support;
            support.clear();
        }
    }

    // Calculate geodesic distance
    unsigned ratio_index = 1;
    double geodesic_distance = 0.0;
    for (auto & AB : support) {
        double dropped_length = opCalcTreeIDLength(AB.first);
        double added_length   = opCalcTreeIDLength(AB.second);
        double ratio = dropped_length/added_length;
        geodesic_distance += pow(dropped_length + added_length, 2);

        if (!_quiet) {
            cout << str(format("\nRatio %d: %.9f") % ratio_index % ratio) << endl;
            cout << "  Edges dropped:" << endl;
            for (auto & a : AB.first) {
                cout << "    " << a.createPatternRepresentation() << endl;
            }
            cout << "  Edges added:" << endl;
            for (auto & b : AB.second) {
                cout << "    " << b.createPatternRepresentation() << endl;
            }
        }
        ++ratio_index;
    }
    geodesic_distance = sqrt(geodesic_distance);
    return geodesic_distance;
}

inline double OP::calcBHVDistance(unsigned ref_index, unsigned test_index) const {
    // Get the reference tree
    string ref_newick = _tree_summary->getNewick(ref_index);
    bool ref_isrooted = _tree_summary->isRooted(ref_index);

    // Get the test tree
    string test_newick = _tree_summary->getNewick(test_index);
    bool test_isrooted = _tree_summary->isRooted(ref_index);

    // Ensure both trees are rooted
    if (!ref_isrooted || !test_isrooted) {
        throw Xop("Trees must be rooted in this version");
    }

    // Build the reference tree
    TreeManip starttm;
    starttm.buildFromNewick(ref_newick, /*rooted*/ref_isrooted, /*allow_polytomies*/false);
    //TODO: get rooted status from treeManip object

    // Store splits from the reference tree
    Split::treeid_t A;
    Split::treeid_t Alvs;
    starttm.storeSplits(A, Alvs);

    if (!_quiet) {
        cout << "Internal splits from starting tree:" << endl;
        for (const auto& a : A) {
            cout << "  " << a.createPatternRepresentation() << endl;
        }
    }

    // Build the test tree
    TreeManip endtm;
    endtm.buildFromNewick(test_newick, /*rooted*/test_isrooted, /*allow_polytomies*/false);

    // Store splits from the reference tree
    Split::treeid_t B;
    Split::treeid_t Blvs;
    endtm.storeSplits(B, Blvs);

    if (!_quiet) {
        cout << "Internal splits from ending tree:" << endl;
        for (const auto& b : B) {
            cout << "  " << b.createPatternRepresentation() << endl;
        }
    }

    // Calculate the contribution of leaf edges to the geodesic
    double leaf_contribution_squared = opCalcLeafContribution(Alvs, Blvs);

    // Find common edges and calculate the contribution of common edge lengths to the geodesic
    vector<Split> common_edges;
    double common_edge_contribution_squared = opFindCommonEdges(A, B, common_edges);

    // Create a vector of paired subtrees by splitting at common edges
    vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
    in_pairs.emplace_back(A,B);
    if (!common_edges.empty())
        opSplitAtCommonEdges(common_edges, in_pairs);

    unsigned pair_index = 1;
    vector<double> geodesic_distances;
    for (auto & inpair : in_pairs) {
        if (!_quiet)
            cout << str(format("\nTree pair %d (of %d)") % pair_index % in_pairs.size()) << endl;

        vector<Split::treeid_pair_t> ABpairs;
        ABpairs.push_back(inpair);

        if (!_quiet) {
            cout << "  A splits:" << endl;
            for (const auto& a : ABpairs[0].first) {
                cout << "    " << a.createPatternRepresentation() << endl;
            }

            cout << "  B splits:" << endl;
            for (const auto& b : ABpairs[0].second) {
                cout << "    " << b.createPatternRepresentation() << endl;
            }
        }

        double L = opCalcGeodesicDist(ABpairs);

        if (!_quiet)
            cout << str(format("  L for tree pair %d = %.9f") % pair_index % L) << endl;

        geodesic_distances.push_back(L);
        ++pair_index;
    }

    if (!_quiet)
        cout << endl;

    // Calculate total geodesic distance
    double total_geodesic_distance = 0.0;
    for (double x : geodesic_distances) {
        total_geodesic_distance += pow(x, 2);
    }
    total_geodesic_distance += leaf_contribution_squared;
    total_geodesic_distance += common_edge_contribution_squared;
    total_geodesic_distance = sqrt(total_geodesic_distance);

    if (!_quiet)
        cout << str(format("Total geodesic distance = %.9f") % total_geodesic_distance) << endl;

    return total_geodesic_distance;
}

inline double OP::calcKFDistance(unsigned ref_index, unsigned test_index) const {
    // Get the reference tree
    string ref_newick = _tree_summary->getNewick(ref_index);
    bool ref_isrooted = _tree_summary->isRooted(ref_index);

    // Get the test tree
    string test_newick = _tree_summary->getNewick(test_index);
    bool test_isrooted = _tree_summary->isRooted(ref_index);

    // Ensure both trees are rooted
    if (!ref_isrooted || !test_isrooted) {
        throw Xop(format("Trees must be rooted in this version of %s") % OP::_program_name);
    }

    // Build the reference tree
    TreeManip reftm;
    reftm.buildFromNewick(ref_newick, /*rooted*/ref_isrooted, /*allow_polytomies*/false);
    //TODO: get rooted status from treeManip object

    // Store splits from the reference tree
    Split::treeid_t refsplits;
    Split::treeid_t reflvs;
    reftm.storeSplits(refsplits, reflvs);

    // Build the test tree
    TreeManip testtm;
    testtm.buildFromNewick(test_newick, /*rooted*/test_isrooted, /*allow_polytomies*/false);

    // Store splits from the reference tree
    Split::treeid_t testsplits;
    Split::treeid_t testlvs;
    testtm.storeSplits(testsplits, testlvs);

    // Store union of refsplits and testsplits in allsplits
    Split::treeid_t allsplits;
    set_union(
        refsplits.begin(), refsplits.end(),
        testsplits.begin(), testsplits.end(),
        inserter(allsplits, allsplits.begin()));
    
    // Traverse allsplits, storing squared branch length differences in KLinternals
    vector<double> KLinternals(allsplits.size());
    unsigned i = 0;
    for (auto s : allsplits) {
        Node * nd0 = reftm.getNodeWithSplit(s);
        Node * nd  = testtm.getNodeWithSplit(s);
        assert(!(nd0 == nullptr && nd == nullptr));
        if (nd0 == nullptr) {
            double edge_length = nd->getEdgeLength();
            double square = pow(edge_length, 2.0);
            KLinternals[i++] = square;
        }
        else if (nd == nullptr) {
            double edge_length = nd0->getEdgeLength();
            double square = pow(edge_length, 2.0);
            KLinternals[i++] = square;
        }
        else {
            double edge_length0 = nd0->getEdgeLength();
            double edge_length  = nd->getEdgeLength();
            double square = pow(edge_length0 - edge_length, 2.0);
            KLinternals[i++] = square;
        }
    }
        
    // Create the map in which keys are taxon names and values are Node pointers
    // for the reference tree
    map<string, Node *> leafmap0;
    reftm.createLeafNodeMap(leafmap0);

    // Create a map in which keys are taxon names and values are Node pointers
    // for the test tree
    map<string, Node *> leafmap;
    testtm.createLeafNodeMap(leafmap);

    // The two trees should have the same number of leaves
    assert(leafmap0.size() == leafmap.size());

    // Get taxon names from the reference tree (assuming the taxon names
    // in the test tree are the same)
    vector<string> names;
    names.reserve(leafmap0.size());
for (const auto& p : leafmap0) {
        names.push_back(p.first);
    }
    sort(names.begin(), names.end());
    
    // Now calculate squares for leaf nodes, storing in KLleaves
    vector<double> KLleaves(names.size());
    i = 0;
    for (const auto& nm : names) {
        Node * nd0 = leafmap0[nm];
        Node * nd  = leafmap[nm];
        double edge_length0 = nd0->getEdgeLength();
        double edge_length  = nd->getEdgeLength();
        double square = pow(edge_length0 - edge_length, 2.0);
        KLleaves[i++] = square;
    }
    
    // Calculate KL distance
    double KLdist = 0.0;
    for (auto square : KLinternals) {
        KLdist += square;
    }
    for (auto square : KLleaves) {
        KLdist += square;
    }
    
    return KLdist;
}

inline void OP::run() {
    try {
        // Read in a tree
        _tree_summary = make_shared<TreeSummary>();
        _tree_summary->readTreefile(_tree_file_name, 0);
        unsigned ntrees = _tree_summary->getNumTrees();
        if (ntrees < 2) {
            throw Xop("Must input at least 2 trees to compute tree distances");
        }

        if (_distance_measure == "geodesic") {
            if (!_quiet)
                cout << "Writing geodesic distances to file \"bhvdists.txt\"" << endl;
            ofstream outf("bhvdists.txt");
            outf << "tree	distance to tree 1" << endl;
            for (unsigned i = 1; i < ntrees; i++) {
                double bhvdist = calcBHVDistance(0, i);
                outf << str(format("%d\t%.5f") % (i+1) % bhvdist) << endl;
            }
            outf.close();
        }
        else if (_distance_measure == "kf") {
            if (!_quiet)
                cout << "Writing KF distances to file \"kfdists.txt\"" << endl;
            ofstream outf("kfdists.txt");
            outf << "tree	distance to tree 1" << endl;
            for (unsigned i = 1; i < ntrees; i++) {
                double kfss = calcKFDistance(0, i);
                double kfdist = sqrt(kfss);
                outf << str(format("%d\t%.5f") % (i+1) % kfdist) << endl;
            }
            outf.close();
        }
    }
    catch (Xop & x) {
        cerr << "Strom encountered a problem:\n  " << x.what() << endl;
        }
}

} // namespace strom

#pragma once

using namespace std;
using namespace boost;

namespace op {

class OP {
public:
    OP();
    ~OP();

    void                clear();
    void                processCommandLineOptions(int argc, const char * argv[]);
    void                run();

private:

    void buildTree(unsigned tree_index, TreeManip & tm) const;
    double calcBHVDistance(TreeManip & starttm, TreeManip & endtm, vector<Split::treeid_pair_t> & ABpairs) const;
    double calcKFDistance(unsigned ref_index, unsigned test_index) const;
    void chooseRandomTree(TreeManip & tm, Lot & lot) const;
    void displaceTreeAlongGeodesic(TreeManip & starttree, TreeManip & endtree, double displacement) const;
    bool frechetCloseEnough(vector<TreeManip> & mu, unsigned lower, unsigned upper, double epsilon) const;
    void computeFrechetMean(TreeManip & meantree) const;
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
    static string opCreateVertexLabel(string name, string capacity, string edgelen, string bipartition);
    static string opCreateEdgeLabel(double capacity, double reverse_flow);
    static void opSaveIncompatibilityGraph(
        OPVertex & source,
        OPVertex & sink,
        vector<OPVertex> & avect,
        vector<OPVertex> & bvect);
#endif
#if 1
    static void opEdmondsKarp(
        OPVertex & source,
        OPVertex & sink,
        vector<OPVertex> & avect,
        vector<OPVertex> & bvect,
        Split::treeid_t & C1,
        Split::treeid_t & C2,
        Split::treeid_t & D1,
        Split::treeid_t & D2,
        bool quiet);
#else
    static void opEdmondsKarp(
        vector<OPVertex> & avect,
        vector<OPVertex> & bvect,
        edgemap_t & edgemap,
        Split::treeid_t & C1,
        Split::treeid_t & C2,
        Split::treeid_t & D1,
        Split::treeid_t & D2,
        bool quiet);
#endif
        bool opRefineSupport(
            const Split::treeid_pair_t & AB,
            Split::treeid_pair_t & AB1,
            Split::treeid_pair_t & AB2) const;
        double opCalcGeodesicDist(
            vector<Split::treeid_pair_t> & ABpairs) const;

        bool                    _quiet;
        bool                    _output_for_gtp;
        bool                    _frechet_mean;
        unsigned                _precision;
        unsigned                _rnseed;
        string                  _tree_file_name;
        string                  _distance_measure;
        TreeSummary::SharedPtr  _tree_summary;

        static string           _program_name;
        static unsigned         _major_version;
        static unsigned         _minor_version;

#if defined(OP_SAVE_DOT_FILE)
    static unsigned            _graph_number;
#endif

    };

inline OP::OP() : _quiet(true), _output_for_gtp(false), _precision(9) {
    //cout << "Constructing a SStrom" << endl;
    clear();
}

inline OP::~OP() = default;

inline void OP::clear() {
    _quiet = true;
    _output_for_gtp = false;
    _frechet_mean = false;
    _tree_file_name = "";
    _distance_measure = "geodesic";
    _tree_summary   = nullptr;
    _precision = 9;
    _rnseed = 1;
}

inline void OP::processCommandLineOptions(int argc, const char * argv[]) {
    program_options::variables_map       vm;
    program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("version,v", "show program version")
        ("treefile,t",  program_options::value(&_tree_file_name)->required(), "name of data file in NEXUS format (required, no default)")
        ("dist", program_options::value(&_distance_measure), "specify either kf or geodesic (default: geodesic)")
        ("precision", program_options::value(&_precision)->default_value(9), "number of digits precision to use in outputting distances (default: 9)")
        ("frechet", program_options::value(&_frechet_mean)->default_value(false), "compute Frechet mean tree using Sturm's algorithm and geodesic distances")
        ("gtptest", program_options::value(&_output_for_gtp)->default_value(false), "output treefile that can be read by Owens-Provan GTP program")
        ("quiet,q", program_options::value(&_quiet), "suppress all output except for errors (default: yes)")
        ("rnseed", program_options::value(&_rnseed), "pseudorandom number generator seed (used only when estimating mean tree)")
        ;
    program_options::store(program_options::parse_command_line(argc, argv, desc), vm);
    try {
        const program_options::parsed_options & parsed = program_options::parse_config_file< char >("op.conf", desc, false);
        program_options::store(parsed, vm);
    }
    catch(program_options::reading_file &) {
        cout << "Note: configuration file (op.conf) not found" << endl;
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
            cout << "  " << a.createPatternRepresentation(true) << endl;
        }

        cout << "Leaves from ending tree:" << endl;
        for (auto & b : Blvs) {
            cout << "  " << b.createPatternRepresentation(true) << endl;
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
    // Find splits in the intersection of A and B
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

    // Count the number of splits in B compatible with each split in A (and vice versa)
    map<const Split *, unsigned> acompatibilities;
    map<const Split *, unsigned> bcompatibilities;
    for (auto & a : A) {
        for (auto & b : B) {
            if (a.compatibleWith(b)) {
                acompatibilities[&a]++;
                bcompatibilities[&b]++;
            }
        }
    }

    // Add splits in A that are compatible with all splits in B to common_edges
    for (auto & apair : acompatibilities) {
        if (apair.second == B.size()) {
            // This split is compatible with every split in the other tree, so add it to the vector
            // of common edges if it is not already in the vector of common edges
            if (find(common_edges.begin(), common_edges.end(), *apair.first) == common_edges.end()) {
                common_edges.push_back(*apair.first);
                common_edge_contribution_squared += pow(apair.first->getEdgeLen(), 2);
            }
        }
    }

    // Add splits in B that are compatible with all splits in A to common_edges
    for (auto & bpair : bcompatibilities) {
        if (bpair.second == A.size()) {
            // This split is compatible with every split in the other tree, so add it to the vector
            // of common edges if it is not already in the vector of common edges
            if (find(common_edges.begin(), common_edges.end(), *bpair.first) == common_edges.end()) {
                common_edges.push_back(*bpair.first);
                common_edge_contribution_squared += pow(bpair.first->getEdgeLen(), 2);
            }
        }
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
        // //temporary!
        // cout << "\ncommon: " << common.createPatternRepresentation() << endl;

        // Create a mask that can be used to zero out all bits in common except the first
        Split mask = common;
        mask.invertBits();
        unsigned first_common_bit = common.findFirstSetBit();
        mask.setBitAt(first_common_bit);

        // //temporary!
        // cout << "  mask: " << mask.createPatternRepresentation() << endl;

        for (auto & inpair : in_pairs) {
            // //temporary!
            // cout << "\n***** new tree pair *****" << endl;
            // Separate out splits in starting (a_splits) vs. ending (b_splits) (sub)trees
            Split::treeid_t & a_splits = inpair.first;
            Split::treeid_t & b_splits = inpair.second;

            // Create split sets to hold splits subsumed in s vs. other splits
            Split::treeid_t a_common_splits, b_common_splits;
            Split::treeid_t a_other_splits, b_other_splits;

            // Divvy up a_splits to a_common_splits and a_other_splits
            // //temporary!
            // cout << "  Divvying up a_splits:" << endl;
            for (auto & asplit : a_splits) {
                // //temporary!
                // cout << "    asplit: " << asplit.createPatternRepresentation();
                bool is_common = (asplit == common);
                if (is_common) {
                    // //temporary!
                    // cout << " (common)" << endl;
                }
                else {
                    if (asplit.subsumedIn(common)) {
                        // //temporary!
                        // cout << " (subsumed in common)" << endl;
                        a_common_splits.insert(asplit);
                    }
                    else {
                        // //temporary!
                        // cout << " (other)" << endl;
                        Split masked = asplit;
                        masked.bitwiseAnd(mask);
                        // //temporary!
                        // cout << "    masked: " << masked.createPatternRepresentation() << endl;
                        a_other_splits.insert(masked);
                    }
                }
            }

            // Divvy up b_splits to b_common_splits and b_other_splits
            // cout << "  Divvying up b_splits:" << endl;
            for (auto & bsplit : b_splits) {
                // cout << "  bsplit: " << bsplit.createPatternRepresentation();
                bool is_common = (bsplit == common);
                if (is_common) {
                    // cout << " (common)" << endl;
                }
                else {
                    if (bsplit.subsumedIn(common)) {
                        // cout << " (subsumed in common)" << endl;
                        b_common_splits.insert(bsplit);
                    }
                    else {
                        // cout << " (other)" << endl;
                        Split masked = bsplit;
                        masked.bitwiseAnd(mask);
                        // cout << "  masked: " << masked.createPatternRepresentation() << endl;
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
inline string OP::opCreateVertexLabel(string name, string capacity, string edgelen, string bipartition) {
    string s = "";
    s += "<<table bgcolor=\"white\" border=\"0\">";
    s += "<tr><td><b><font color=\"blue\" face=\"Courier\" point-size=\"16\">%s</font></b></td></tr>";
    s += "<tr><td><b><font color=\"blue\" face=\"Courier\" point-size=\"16\">%s</font></b></td></tr>";
    s += "<tr><td><font color=\"black\" face=\"Courier\" point-size=\"10\">%s</font></td></tr>";
    s += "<tr><td><font color=\"black\" face=\"Courier\" point-size=\"10\">%s</font></td></tr>";
    s += "</table>>";
    return str(format(s) % name % capacity % bipartition % edgelen);
}

inline string OP::opCreateEdgeLabel(double capacity, double reverse_flow) {
    string s = "";
    s += "\tlabeldistance=8\n";
    s += "\tlabelangle=0\n";
#if 1
    s += "\theadlabel=<<font color=\"red\" face=\"Verdana\" point-size=\"12\">%.3f</font>>\n";
#else
    s += "\theadlabel=<\n";
    s += "\t\t<table bgcolor=\"white\" border=\"0\">\n";
    s += "\t\t\t<tr>\n";
    s += "\t\t\t\t<td>\n";
    s += "\t\t\t\t\t<font color=\"blue\" face=\"Courier\" point-size=\"10\">%.3f</font>\n";
    s += "\t\t\t\t</td>\n";
    s += "\t\t\t</tr>\n";
    s += "\t\t\t<tr>\n";
    s += "\t\t\t\t<td>\n";
    s += "\t\t\t\t\t<font color=\"red\" face=\"Courier\" point-size=\"10\">%.3f</font>\n";
    s += "\t\t\t\t</td>\n";
    s += "\t\t\t</tr>\n";
    s += "\t\t</table>\n";
    s += "\t>\n";
#endif
    return str(format(s) % reverse_flow);
}

inline void OP::opSaveIncompatibilityGraph(OPVertex & source, OPVertex & sink, vector<OPVertex> & avect, vector<OPVertex> & bvect) {
    // Example of the kind of dot file generated by this function:
    // digraph G {
    //     rankdir=LR;
    //     graph [ranksep=2];
    //
    //     subgraph Avertices {
    //         label="A";
    //             {
    //             rank=same;
    //             a1 [label = "0.349\na1\n----**-", shape = box, color = black];
    //             a2 [label = "0.100\na2\n--*-**-", shape = box, color = black];
    //             a3 [label = "0.240\na3\n-**-**-", shape = box, color = black];
    //             }
    //         a1 -> a2 [dir=none, style=invisible];
    //         a2 -> a3 [dir=none, style=invisible];
    //     }
    //
    //     subgraph Bvertices {
    //         label="B";
    //             {
    //             rank=same;
    //             b1 [label = "0.016\nb1\n--*-*--", shape = box, color = black];
    //             b2 [label = "0.524\nb2\n-**-*--", shape = box, color = black];
    //             b3 [label = "0.122\nb3\n-----**", shape = box, color = black];
    //             b4 [label = "0.339\nb4\n-**-***", shape = box, color = black];
    //             }
    //         b1 -> b2 [dir=none, style=invisible];
    //         b2 -> b3 [dir=none, style=invisible];
    //         b3 -> b4 [dir=none, style=invisible];
    //     }
    //
    //     a1 -> b1;
    //     a1 -> b2;
    //     a1 -> b3;
    //     a2 -> b2;
    //     a2 -> b3;
    //     a3 -> b3;
    //     a3 -> b4;
    // }
    //
    //

    auto asize = static_cast<unsigned>(avect.size());
    auto bsize = static_cast<unsigned>(bvect.size());
    unsigned minsize = min(asize, bsize);
    unsigned maxsize = max(asize, bsize);

    // Save all the entities needed
    // tuple key: <0> name, <1>capacity, <2>edgelen, <3>split, <4>shape, <5>color
    typedef std::tuple<string, string, string, string, string, string> gnode_element_t;
    typedef vector<gnode_element_t> gnode_t;

    // Save everything needed for the A nodes
    gnode_t anodes;
    for (unsigned i = 0; i < source._edges.size(); i++) {
        OPEdge * edge   = source._edges[i];
        string name     = edge->_to->_name;
        string capacity = (edge->_capacity == 0.0 ? "0" : str(format("%.3f") % edge->_capacity));
        string edgelen  = str(format("%.3f") % edge->_to->_split->getEdgeLen());
        string split    = str(format("%.3f") % edge->_to->_split->createPatternRepresentation(false));
        string shape    = (edge->_capacity == 0.0 ? "circle" : "box");
        string color    = (edge->_capacity == 0.0 ? "red" : "black");
        anodes.emplace_back(name, capacity, edgelen, split, shape, color);
    }
    if (asize < maxsize) {
        for (unsigned i = asize; i < maxsize; i++) {
            anodes.emplace_back(str(format("adummy%d") % (i+1)), "0", "0", "0", "box", "black");
        }
    }

    // Save everything needed for the B nodes
    gnode_t bnodes;
    for (unsigned i = 0; i < bsize; i++) {
        assert(bvect[i]._edges.size() == 1);
        OPEdge * edge   = bvect[i]._edges[0];
        string name     = edge->_from->_name;
        string capacity = (edge->_capacity == 0.0 ? "0" : str(format("%.3f") % edge->_capacity));
        string edgelen  = str(format("%.3f") % bvect[i]._split->getEdgeLen());
        string split    = str(format("%.3f") % bvect[i]._split->createPatternRepresentation(false));
        string shape    = (edge->_capacity == 0.0 ? "circle" : "box");
        string color    = (edge->_capacity == 0.0 ? "red" : "black");
        bnodes.emplace_back(name, capacity, edgelen, split, shape, color);
    }
    if (bsize < maxsize) {
        for (unsigned i = bsize; i < maxsize; i++) {
            bnodes.emplace_back(str(format("bdummy%d") % (i+1)), "0", "0", "0", "box", "black");
        }
    }

    // Create (or append to) a rundot.sh file containing commands to compile all the incompatibility graphs
    if (OP::_graph_number == 1) {
        ofstream shf("rundot.sh");
        shf << "#!/bin/bash\n\n";
        shf << "# This file requires previous installation of dot and is intended\n";
        shf << "# to be run on a Mac (because of the use of the open command)\n\n";
        shf << "dot -Tpng graph-1.dot > graph-1.png; open graph-1.png\n";
        shf.close();
    }
    else {
        ofstream shf("rundot.sh", ios::out | ios::app);
        shf << str(format("dot -Tpng graph-%d.dot > graph-%d.png; open graph-%d.png\n")
            % OP::_graph_number
            % OP::_graph_number
            % OP::_graph_number);
        shf.close();
    }

    // Open the dot file
    cout << "Saving incompatibility graph to file: " << str(format("graph-%d.dot") % OP::_graph_number) << endl;
    ofstream dotf(str(format("graph-%d.dot") % OP::_graph_number++));

    // Create the opening preamble
    dotf << "digraph G {\n";
    dotf << "\trankdir=LR;\n";
    dotf << "\tgraph [ranksep=2];\n\n";

    // Create a subgraph containing only the A vertices in a vertical rank
    dotf << "\tsubgraph Avertices {\n";
    dotf << "\t    label=\"A\";\n";
    dotf << "\t    {\n";
    dotf << "\t        rank=same;\n";
    for (unsigned i = 0; i < maxsize; i++) {
        if (i < asize) {
            // tuple key: <0> name, <1>capacity, <2>edgelen, <3>split, <4>shape, <5>color

            dotf << str(format("\t        %s [label = %s, shape = %s, color = %s];\n")
                % get<0>(anodes[i])
                % opCreateVertexLabel(get<0>(anodes[i]), get<1>(anodes[i]), get<2>(anodes[i]), get<3>(anodes[i]))
                % get<4>(anodes[i])
                % get<5>(anodes[i])
            );
        }
        else {
            dotf << str(format("\t        %s [style = invisible];\n")
                % get<0>(anodes[i])
            );
        }
    }
    dotf << "\t    }\n";
    for (unsigned i = 0; i < maxsize - 1; i++) {
        dotf << str(format("\t    %s -> %s [dir=none, style=invisible];\n")
                % get<0>(anodes[i])
                % get<0>(anodes[i+1])
        );
    }
    dotf << "\t}\n";

    // Create a subgraph containing only the B vertices in a vertical rank
    dotf << "\tsubgraph Bvertices {\n";
    dotf << "\t    label=\"B\";\n";
    dotf << "\t    {\n";
    dotf << "\t        rank=same;\n";
    for (unsigned i = 0; i < maxsize; i++) {
        if (i < bsize) {
            dotf << str(format("\t        %s [label = %s, shape = %s, color = %s];\n")
                % get<0>(bnodes[i])
                % opCreateVertexLabel(get<0>(bnodes[i]), get<1>(bnodes[i]), get<2>(bnodes[i]), get<3>(bnodes[i]))
                % get<4>(bnodes[i])
                % get<5>(bnodes[i])
            );
        }
        else {
            dotf << str(format("\t        %s [style = invisible];\n")
                % get<0>(bnodes[i])
            );
        }
    }
    dotf << "\t    }\n";
    for (unsigned i = 0; i < maxsize - 1; i++) {
        dotf << str(format("\t    %s -> %s [dir=none, style=invisible];\n")
                % get<0>(bnodes[i])
                % get<0>(bnodes[i+1])
        );
    }
    if (bsize < maxsize) {
        dotf << str(format("\t    %s -> bdummy%d [dir=none, style=invisible];\n")
            % get<0>(anodes[bsize-1])
            % (bsize+1));
        for (unsigned i = bsize; i < maxsize - 1; i++) {
            dotf << str(format("\t    bdummy%d -> bdummy%d [dir=none, style=invisible];\n") % (i+1) % (i+2));
        }
    }
    dotf << "\t}\n";

    // Create the edges connecting A with B vertices
    for (unsigned i = 0; i < asize; ++i) {
        for (unsigned j = 0; j < avect[i]._edges.size(); ++j) {
            OPEdge * edge = avect[i]._edges[j];
            dotf << str(format("\t\t%s -> %s [%s];\n") % edge->_from->_name % edge->_to->_name % opCreateEdgeLabel(edge->_capacity, edge->_reverse_flow));
        }
    }

    // dotf << "\t}\n";
    dotf << "}\n";
    dotf.close();
}
#endif

#if 1
    inline void OP::opEdmondsKarp(
            OPVertex & source,
            OPVertex & sink,
            vector<OPVertex> & avect,
            vector<OPVertex> & bvect,
            Split::treeid_t & C1,
            Split::treeid_t & C2,
            Split::treeid_t & D1,
            Split::treeid_t & D2,
            bool quiet) {
#if defined(OP_SAVE_DOT_FILE)
    opSaveIncompatibilityGraph(source, sink, avect, bvect);
#endif

        double cumulative_flow = 0.0;
        bool done_augmenting_path = false;
        while (!done_augmenting_path) {
            vector<OPVertex *> route;

            // All vertices begin as unvisited, have residual capacity 0, and all "A" to "B" edges begin as non-reversed
            source._parent_edge = nullptr;
            for (auto & a : avect) {
                a._parent_edge = nullptr;
                a._residual_capacity = 0.0;
                for (auto & a_edge : a._edges) {
                    a_edge->_edge_is_reversed = false;
                }
            }
            for (auto & b : bvect) {
                b._parent_edge = nullptr;
            }
            sink._parent_edge = nullptr;

            //


            // Add the source to the route
            route.push_back(&source);
            unsigned route_cursor = 0;

            // Find the next augmenting route
            bool sink_found = false;
            while (!sink_found) {
                OPVertex * current = route[route_cursor];
                for (auto & edge : current->_edges) {
                    // //temporary!
                    // cout << "checking " << edge->_from->_name << " -> " << edge->_to->_name << endl;
                    if (edge->_capacity > 0.0 && edge->_to->_parent_edge == nullptr) {
                        edge->_to->_parent_edge = edge;

                        // //temporary!
                        // cout << "  adding edge to route" << endl;

                        route.push_back(edge->_to);
                        if (edge->_to == &sink) {
                            sink_found = true;
                            break;
                        }
                    }
                    else {
                        bool already_visited = edge->_to->_parent_edge != nullptr;
                        bool capacity_zero = edge->_capacity == 0.0;
                        if (already_visited && capacity_zero) {
                            // //temporary!
                            // cout << "  rejected edge because capacity is zero and to-vertex already visited" << endl;
                        }
                        else if (already_visited) {
                            // //temporary!
                            // cout << "  rejected edge because to-vertex already visited" << endl;
                        }
                        else {
                            // //temporary!
                            // cout << "  rejected edge because capacity is zero" << endl;

                            // If it is a "B" vertex that has zero capacity, then the route cannot
                            // go from thix "B" vertex to the sink but it may be able to go back to
                            // an "A" vertex if that "A" vertex was not accessible from the source
                            // and there is residual flow on the edge
                            if (edge->_to == &sink) {
                                OPVertex * B_vertex = edge->_from;
                                for (unsigned j = 0; j < source._edges.size(); j++) {
                                    OPEdge * s_to_A_edge = source._edges[j];  // edge from source to an "A" vertex
                                    if (s_to_A_edge->_capacity == 0) {
                                        // See if any of the "A" vertex edges lead to the "B" vertex in question
                                        OPVertex * A_vertex = s_to_A_edge->_to;
                                        for (unsigned k = 0; k < A_vertex->_edges.size(); k++) {
                                            OPEdge * A_to_B_edge = A_vertex->_edges[k];  // edge from "A" to a "B" vertex
                                            if (A_vertex->_parent_edge == nullptr && A_to_B_edge->_to == B_vertex && A_to_B_edge->_reverse_flow > 0.0) {
                                                // If we are here, we know that:
                                                // 1. the "A" vertex has not been visited
                                                // 2. the "A" vertex is joined to the "B" vertex in question
                                                // 3. the "A" vertex is not accessible from the source, and
                                                // 4. there is residual flow on the A_to_B_edge

                                                // //temporary!
                                                // cout << "  adding REVERSE edge to route (" << B_vertex->_name << " -> " << A_vertex->_name << ")" << endl;

                                                // //temporary!
                                                // if (A_to_B_edge->_from->_name == "a2" && A_to_B_edge->_to->_name == "b4") {
                                                //     cerr << "*** debug breakpoint ***" << endl;
                                                // }

                                                A_vertex->_parent_edge = A_to_B_edge;
                                                A_to_B_edge->_edge_is_reversed = true;
                                                A_vertex->_residual_capacity = A_to_B_edge->_reverse_flow;
                                                route.push_back(A_vertex);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                route_cursor++;
                if (route_cursor >= route.size()) {
                    // Not able to reach sink, so all augmenting routes have been found
                    done_augmenting_path = true;
                    break;
                }
            }

            if (!done_augmenting_path) {
                // Follow _from pointers from sink back to source to identify
                // the bottleneck and determine flow through this route
                double min_capacity = 1.0;
                OPVertex * current = &sink;
                while (current != &source) {
                    // //temporary!
                    //cerr << "current = " << current->_name << endl;
                    OPEdge * edge = current->_parent_edge;
                    if (edge->_edge_is_reversed) {
                        if (edge->_reverse_flow < min_capacity) {
                            min_capacity = edge->_reverse_flow;
                        }
                        current = edge->_to;
                    }
                    else {
                        if (edge->_capacity < min_capacity) {
                            min_capacity = edge->_capacity;
                        }
                        current = edge->_from;
                    }
                }

                // Adjust capacities and flows along the route
                cumulative_flow += min_capacity;
                current = &sink;
                while (current != &source) {
                    OPEdge * edge = current->_parent_edge;
                    if (edge->_edge_is_reversed) {
                        edge->_reverse_flow -= min_capacity;
                        if (fabs(edge->_reverse_flow) < 1e-10)
                            edge->_reverse_flow = 0.0;
                        edge->_capacity += min_capacity;
                        current = edge->_to;
                    }
                    else {
                        edge->_capacity -= min_capacity;
                        if (fabs(edge->_capacity) < 1e-10)
                            edge->_capacity = 0.0;
                        edge->_reverse_flow += min_capacity;
                        current = edge->_from;
                    }
                }
            }
#if defined(OP_SAVE_DOT_FILE)
        opSaveIncompatibilityGraph(source, sink, avect, bvect);

        // //temporary!
        // cerr << "just saved dot graph" << endl;
#endif
        }   // while (!done_augmenting_path)

    // Identify C1, C2, D1, and D2
    // C1 and D2 compose the min weight vertex cover
    // C2 and D` compose the independent set
    for (auto & source_edge : source._edges) {
        OPVertex * avertex = source_edge->_to;
        if (source_edge->_capacity > 0.0 || avertex->_residual_capacity > 0.0) {
            // Because this source edge has remaining capacity, its distal vertex is part of the independent set
            C2.insert(*(avertex->_split));

            // This A vertex allows access to the B side, so any connected B vertices with
            // zero capacity are part of the vertex cover
            for (const auto central_edge : avertex->_edges) {
                OPVertex * bvertex = central_edge->_to;
                assert(bvertex->_edges.size() == 1);
                OPEdge * sink_edge = bvertex->_edges[0];
                if (sink_edge->_capacity == 0.0) {
                    D2.insert(*(bvertex->_split));
                }
            }
        }
        else {
            // Because this A-vertex has zero capacity, it is part of the vertex cover
            C1.insert(*(avertex->_split));

            // No need to consider connected B vertices because this A-vertex already
            // covers all connected edges
        }
    }

    // D1 includes every split not in D2
    D1.clear();
    for (auto & b : bvect) {
        if (D2.count(*(b._split)) == 0) {
            D1.insert(*(b._split));
        }
    }
}
#else
inline void OP::opEdmondsKarp(
        vector<OPVertex> & avect,
        vector<OPVertex> & bvect,
        edgemap_t & edgemap,
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

    bool done = false;
    while (!done) {
        vector<OPVertex *> route;

        // Make sure none of the "b" vertices are marked as visited
        for (auto & b : bvect) {
            b._parent_index = -1;
        }

        // Add all avect vertices to the route if they have capacity > 0
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
                if (b->_capacity > 0.0 && b->_parent_index == -1) {
                    b->_parent_index = static_cast<int>(aindex);
                    route.push_back(b);
                }
                else if (b->_parent_index == -1) {
                    // B is accessible but has outgoing capacity zero; see if any viable route exists to the left
                    // Check to see if this ever happens; if it does happen, more coding needs to happen
                    for (unsigned ai = 0; ai < route_size; ai++) {
                        if (route[ai]->_capacity == 0.0 && edgemap.at(make_pair(route[ai], b)) > 0.0) {
                            throw Xop("error: need to handle reverse flow in Edmonds-Karp algorithm");
                        }
                    }
                }
            }
        }

        // Find B vertex that first reaches the sink
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
            // Identify the route
            auto avertex = route[last->_parent_index];
            auto bvertex = last;
            auto abpair = make_pair(avertex, bvertex);

            // Find minimum capacity along the route
            double min_capacity = last->_capacity;
            if (route[last->_parent_index]->_capacity < min_capacity) {
                min_capacity = route[last->_parent_index]->_capacity;
            }

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

            // Reduce capacity of the leftmost edge on the route by an amount min_capacity
            avertex->_capacity -= min_capacity;
            if (fabs(avertex->_capacity) < 1e-10) {
                avertex->_capacity = 0.0;
            }

            // Increase the flow on the focal edge by min_capacity
            edgemap.at(abpair) += min_capacity;

            // Reduce capacity of the rightmost edge on the route by an amount min_capacity
            bvertex->_capacity -= min_capacity;
            if (fabs(bvertex->_capacity) < 1e-10) {
                bvertex->_capacity = 0.0;
            }

#if defined(OP_SAVE_DOT_FILE)
            if (!quiet) {
                // Uncomment the line below to save graph.dot (visualization of incompatibility graph in dot language)
                opSaveIncompatibilityGraph(avect, bvect);
            }
#endif
        }
    }

    // Identify C1, C2, D1, and D2
    // C1 and D2 compose the min weight vertex cover
    // C2 and D` compose the independent set
    for (auto & a : avect) {
        if (a._capacity > 0.0) {
            // Because this A-vertex has residual capacity, it is part of the independent set
            C2.insert(*(a._split));

            // This A vertex allows access to the B side, so any connected B vertices with
            // zero capacity are part of the vertex cover
            for (const auto b : a._children) {
                if (b->_capacity == 0.0) {
                    D2.insert(*(b->_split));
                }
            }
        }
        else {
            // Because this A-vertex has zero capacity, it is part of the vertex cover
            C1.insert(*(a._split));

            // No need to consider connected B vertices because this A-vertex already
            // covers all connected edges
        }
    }

    // D1 includes every split not in D2
    D1.clear();
    for (auto & b : bvect) {
        if (D2.count(*(b._split)) == 0) {
            D1.insert(*(b._split));
        }
    }
}
#endif

inline bool OP::opRefineSupport(const Split::treeid_pair_t & AB, Split::treeid_pair_t & AB1, Split::treeid_pair_t & AB2) const {
    // Create a vector of incompatibility graph vertices
    vector<OPVertex> avect(AB.first.size());
    vector<OPVertex> bvect(AB.second.size());

    //                                      23  123  456  123456
    //vector<unsigned> a_splits_in_order = { 6,   7,  56,     63};

    // Calculate weights for the "A" vertices.
    // For example, if A = {a1,a2,a3}, then the weight of a1 is
    // weight a1 = a1^2 / (a1^2 + a2^2 + a3^2)
    unsigned aindex = 0;
    double asum = 0.0;
    for (auto & a : AB.first) {
        asum += pow(a.getEdgeLen(),2);
    }
    for (auto & a : AB.first) {
        avect[aindex]._split = &a;
        avect[aindex]._weight = pow(a.getEdgeLen(),2)/asum;
        aindex++;
    }

    // Calculate weights for the "B" vertices.
    // For example, if B = {b1,b2,b3}, then the weight of b1 is
    // weight b1 = b1^2 / (b1^2 + b2^2 + b3^2)
    unsigned bindex = 0;
    double bsum = 0.0;
    for (auto & b : AB.second) {
        bsum += pow(b.getEdgeLen(),2);
    }
#if 0 //temporary!
    //                                    25  67  2567  34  234567
    vector<unsigned> b_splits_in_order = {18, 96,  114, 12,    126};
    for (auto & b : AB.second) {
        unsigned long bvalue = b.getBits()[0];
        bindex = 99;
        for (unsigned z = 0; z < b_splits_in_order.size(); z++) {
            if (bvalue == b_splits_in_order[z]) {
                bindex = z;
                break;
            }
        }
        assert(bindex < 99);
        bvect[bindex]._split = &b;
        bvect[bindex]._weight = pow(b.getEdgeLen(),2)/bsum;
    }
#else
    for (auto & b : AB.second) {
        bvect[bindex]._split = &b;
        bvect[bindex]._weight = pow(b.getEdgeLen(),2)/bsum;
        bindex++;
    }
#endif

    // Create the incompatibility graph
    // A vertices go on the left, B vertices go on the right, and edges connect an A vertex
    // to a B vertex only if the two vertices are incompatible.
    vector<OPEdge> all_edges;
    all_edges.reserve(avect.size() * bvect.size() + avect.size() + bvect.size());
    unsigned nincompatibilities = 0;
    auto asize = static_cast<unsigned>(avect.size());
    auto bsize = static_cast<unsigned>(bvect.size());

    // Create edges from source to the A-vertices
    OPVertex source;
    source._name = "source";
    for (unsigned i = 0; i < asize; i++) {
        // Assign a name to avect[i]
        avect[i]._name = str(format("a%d") % i);

        // //temporary!
        // cerr << "avect[" << i << "] split: \n";
        // for (auto v : avect[i]._split->getBits()) {
        //     cerr << " " << v << endl;
        // }

        // Create the forward edge
        all_edges.emplace_back();
        OPEdge & source_forward_edge = all_edges.back();
        source_forward_edge._from = &source;
        source_forward_edge._to = &avect[i];
        source_forward_edge._capacity = avect[i]._weight;
        source_forward_edge._flow = 0.0;
        source_forward_edge._reverse_flow = 0.0;
        source_forward_edge._open = true;
        source._edges.push_back(&source_forward_edge);
    }

    // Create edges from A-vertices to B-vertices
    for (unsigned i = 0; i < asize; i++) {
        // Get split associated with avect[i]
        const Split * a = avect[i]._split;
        assert(a);
        for (unsigned j = 0; j < bsize; j++) {
            // Get split associated with bvect[j]
            const Split * b = bvect[j]._split;
            assert(b);
            if (!a->compatibleWith(*b)) {
                nincompatibilities++;

                // Create a forward edge in the incompatibility graph
                all_edges.emplace_back();
                OPEdge & forward_edge = all_edges.back();
                forward_edge._from = &avect[i];
                forward_edge._to = &bvect[j];
                forward_edge._capacity = 1.0;
                forward_edge._flow = 0.0;
                forward_edge._reverse_flow = 0.0;
                forward_edge._open = true;
                avect[i]._edges.emplace_back(&forward_edge);
            }
        }
    }

    // Create edges from the B-vertices to the sink
    OPVertex sink;
    sink._name = "sink";
    for (unsigned j = 0; j < bsize; j++) {
        // Assign a name to bvect[j]
        bvect[j]._name = str(format("b%d") % j);

        // //temporary!
        // cerr << "bvect[" << j << "] split: \n";
        // for (auto v : bvect[j]._split->getBits()) {
        //     cerr << " " << v << endl;
        // }

        // Create the forward edge
        all_edges.emplace_back();
        OPEdge & sink_forward_edge = all_edges.back();
        sink_forward_edge._from = &bvect[j];
        sink_forward_edge._to = &sink;
        sink_forward_edge._capacity = bvect[j]._weight;
        sink_forward_edge._flow = 0.0;
        sink_forward_edge._reverse_flow = 0.0;
        sink_forward_edge._open = true;
        bvect[j]._edges.emplace_back(&sink_forward_edge);
    }

    bool success = false;
    if (nincompatibilities < asize*bsize) {
        // At least one independent pair of vertices exists
        // Carry out Edmonds-Karp algorithm to find min-weight vertex cover (identifies max weight independent set)
        // In Owens-Provan terminology,
        //   C1 = A's contribution to vertex cover     (equals AB1.first)
        //   C2 = A's contribution to independent set  (equals AB2.first)
        //   D1 = B's contribution to independent set  (equals AB1.second)
        //   D2 = B's contribution to vertex cover     (equals AB2.second)
        // A1 and A2 must represent a non-trivial partition of A
        // B1 and B2 must represent a non-trivial partition of B
        // C2 and D1 (AB2.first and AB1.second) are compatible sets of splits
        // C1 and D2 (AB1.first and AB2.second) compose the minimum weight vertex cover
        // ||C1||/||D1|| < ||C2||/||D2|| must be true
        // If all of the above conditions hold, then the support can be refined:
        // A = (A1,A2) and B = (B1,B2)
        // where A1 = C1, A2 = C2, B1 = D1, and B2 = D2
        // Length of this segment of the geodesic is
        //   L = sqrt{ (||A1|| + ||B1||)^2 +  (||A2|| + ||B2||)^2 }
        // Orthants crossed:
        //   start:        A1, A2
        //      edges in A1 decreasing, edges in B1 increasing
        //   intermediate: B1, A2
        //      edges in A2 decreasing, edges in B2 increasing
        //   finish:       B1, B2
        // The point at which A1 edges become 0 is
        //   ||A1||/(||A1|| + ||B1||)
        // The point at which A2 edges become 0 is
        //   ||A2||/(||A2|| + ||B2||)
        opEdmondsKarp(source, sink, avect, bvect, AB1.first, AB2.first, AB1.second, AB2.second, _quiet);

        // if (!quiet) {
        //     cout << "\nResults:" << endl;
        //     cout << str(format("  ||C1|| = %.9f") % C1len) << endl;
        //     cout << str(format("  ||C2|| = %.9f") % C2len) << endl;
        //     cout << str(format("  ||D1|| = %.9f") % D1len) << endl;
        //     cout << str(format("  ||D2|| = %.9f") % D2len) << endl;
        //     cout << "  Check whether ||C1||/||D1|| < ||C2||/||D2||" << endl;
        //     cout << str(format("    %.9f < %.9f") % (C1len/D1len) % (C2len/D2len));
        //     cout << endl;
        // }

        // Check conditions for successful refinement
        bool A_is_trivial = (AB1.first.empty()) || (AB2.first.empty());
        bool B_is_trivial = (AB1.second.empty())|| (AB2.second.empty());
        double C1len = opCalcTreeIDLength(AB1.first);
        double C2len = opCalcTreeIDLength(AB2.first);
        double D1len = opCalcTreeIDLength(AB1.second);
        double D2len = opCalcTreeIDLength(AB2.second);
        bool P2 = C1len/D1len < C2len/D2len;
        success = !A_is_trivial && !B_is_trivial && P2;

        if (!_quiet) {
            if (success) {
                cout << "\nSuccessfully refined support:" << endl;
            }
            else {
                cout << "\nSupport not refined because:" << endl;
                if (A_is_trivial) {
                    cout << "  A is a trivial partition:" << endl;
                }
                if (B_is_trivial) {
                    cout << "  B is a trivial partition:" << endl;
                }
                if (!P2) {
                    cout << "  ||C1||/||D1|| is not strictly less than ||C2||/||D2||" << endl;
                }
            }
            // cout << "  Input A vertices:" << endl;
            // for (auto & a : AB.first) {
            //     cout << "    " << a.createPatternRepresentation() << endl;
            // }
            // cout << "  Input B vertices:" << endl;
            // for (auto & b : AB.second) {
            //     cout << "    " << b.createPatternRepresentation() << endl;
            // }
            cout << "  Output A1 vertices:" << endl;
            if (AB1.first.empty()) {
                cout << "    empty set" << endl;
            }
            else {
                for (auto & a : AB1.first) {
                    cout << "    " << a.createPatternRepresentation(true) << endl;
                }
            }
            cout << "  Output B1 vertices:" << endl;
            if (AB1.second.empty()) {
                cout << "    empty set" << endl;
            }
            else {
                for (auto & b : AB1.second) {
                    cout << "    " << b.createPatternRepresentation(true) << endl;
                }
            }
            cout << "  Output A2 vertices:" << endl;
            if (AB2.first.empty()) {
                cout << "    empty set" << endl;
            }
            else {
                for (auto & a : AB2.first) {
                    cout << "    " << a.createPatternRepresentation(true) << endl;
                }
            }
            cout << "  Output B2 vertices:" << endl;
            if (AB2.second.empty()) {
                cout << "    empty set" << endl;
            }
            else {
                for (auto & b : AB2.second) {
                    cout << "    " << b.createPatternRepresentation(true) << endl;
                }
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

inline void OP::buildTree(unsigned tree_index, TreeManip & tm) const {
    string newick = _tree_summary->getNewick(tree_index);

    bool isrooted = _tree_summary->isRooted(tree_index);
    if (!isrooted) {
        throw Xop("Trees must be rooted in this version");
    }
    tm.buildFromNewick(newick, /*rooted*/isrooted, /*allow_polytomies*/true);
}

inline double OP::calcBHVDistance(TreeManip & starttm, TreeManip & endtm, vector<Split::treeid_pair_t> & ABpairs) const {
    ABpairs.clear();

    // Store splits from the reference tree
    Split::treeid_t A0;
    Split::treeid_t Alvs;
    starttm.storeSplits(A0, Alvs);

    // Only save splits with non-zero edge length
    Split::treeid_t A;
    for (auto & a : A0) {
        if (a.getEdgeLen() > 0.0) {
            A.insert(a);
        }
    }

    if (!_quiet) {
        cout << "Internal splits from starting tree:" << endl;
        for (const auto& a : A) {
            cout << "  " << a.createPatternRepresentation(true) << endl;
        }
    }

    // Store splits from the reference tree
    Split::treeid_t B0;
    Split::treeid_t Blvs;
    endtm.storeSplits(B0, Blvs);

    // Only save splits with non-zero edge length
    Split::treeid_t B;
    for (auto & b : B0) {
        if (b.getEdgeLen() > 0.0) {
            B.insert(b);
        }
    }

    if (!_quiet) {
        cout << "Internal splits from ending tree:" << endl;
        for (const auto& b : B) {
            cout << "  " << b.createPatternRepresentation(true) << endl;
        }
    }

    // Calculate the contribution of leaf edges to the geodesic
    double leaf_contribution_squared = opCalcLeafContribution(Alvs, Blvs);

    // Find common edges and calculate the contribution of common edge lengths to the geodesic
    vector<Split> common_edges;
    double common_edge_contribution_squared = opFindCommonEdges(A, B, common_edges);

#if 0
    // Create a vector of paired subtrees by splitting at common edges
    vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
    in_pairs.emplace_back(A,B);
    if (!common_edges.empty())
        opSplitAtCommonEdges(common_edges, in_pairs);
#else
    if (!common_edges.empty()) {
        // Remove splits representing common edge from both A and B
        for (auto & c : common_edges) {
            A.erase(c);
            B.erase(c);
        }
    }
    vector<pair<Split::treeid_t, Split::treeid_t> > in_pairs;
    in_pairs.emplace_back(A,B);
#endif

    unsigned pair_index = 1;
    vector<double> geodesic_distances;
    for (auto & inpair : in_pairs) {
        if (!_quiet)
            cout << str(format("\nTree pair %d (of %d)") % pair_index % in_pairs.size()) << endl;

        ABpairs.push_back(inpair);

        if (!_quiet) {
            cout << "  A splits:" << endl;
            for (const auto& a : ABpairs[0].first) {
                cout << "    " << a.createPatternRepresentation(true) << endl;
            }

            cout << "  B splits:" << endl;
            for (const auto& b : ABpairs[0].second) {
                cout << "    " << b.createPatternRepresentation(true) << endl;
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
    reftm.buildFromNewick(ref_newick, /*rooted*/ref_isrooted, /*allow_polytomies*/true);
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

inline void OP::chooseRandomTree(TreeManip & tm, Lot & lot) const {
    unsigned n = _tree_summary->getNumTrees();
    unsigned index = lot.randint(0, n-1);
    string newick = _tree_summary->getNewick(index);
    bool rooted = _tree_summary->isRooted(index);
    assert(rooted);
    tm.buildFromNewick(newick, rooted, /*allow_polytomies*/true);
}

inline void OP::displaceTreeAlongGeodesic(TreeManip & starttree, TreeManip & endtree, double displacement) const {
    // Move starttree a distance displacement along the geodesic from starttree to endtree
    // First, get geodesic
    vector<Split::treeid_pair_t> ABpairs;
    double geodesic = calcBHVDistance(starttree, endtree, ABpairs);

    // Support A = (A1, A2, ..., Ak) and B = (B1, B2, ..., Bk)
    // Path has i-1,2,...,k "legs"
    //    G0 if lambda/(1-lambda) <= length(A1)/length(B1)
    //    Gi if length(Ai)/length(Bi) = lambda/(1-lambda) <= length(Ai+1)/length(Bi+1)
    //    Gk if length(Ak)/length(Bk) <= lambda/(1-lambda)
    // and tree Ti along path has edge sets
    //    B1 U ... U Bi U Ai+1 U ... U Ak
    // and edge sets
    //    length edge e in Ti = [(1-lambda) length(Aj) - lambda length(Bj)]/length(Aj) if e in Aj
    //    length edge e in Ti = [lambda length(Bj) - (1-lambda) length(Aj)]/length(Bj) if e in Bj

    // Precalculate lengths of all segments
    auto support_size = (unsigned)ABpairs.size();
    vector<double> lenA(support_size);
    vector<double> lenB(support_size);
    for (unsigned i = 0; i < support_size; ++i) {
        lenA[i] = opCalcTreeIDLength(ABpairs[i].first);
        lenB[i] = opCalcTreeIDLength(ABpairs[i].second);
    }

    // First step is to determine the leg
    double lambda_ratio = displacement/(1.0 - displacement);
    unsigned leg = 0;
    double ratioi = lenA[0]/lenB[0];
    if (lambda_ratio > ratioi) {
        for (unsigned iplus1 = 1; iplus1 < ABpairs.size(); ++iplus1) {
            double ratioiplus1 = lenA[iplus1]/lenB[iplus1];
            if (lambda_ratio >= ratioi && lambda_ratio <= ratioiplus1) {
                leg = iplus1 - 1;
                break;
            }
        }
    }

    // Walk down ABpairs, dropping and adding splits as needed from starttree until we arrive at the destination leg
    double lambda = displacement;
    Split::treeid_t splitset;
    for (unsigned i = 0; i <= leg; ++i) {
        double edgelen_multiplicative_factor = (lambda*lenB[i] - (1.0 - lambda)*lenA[i])/lenB[i];
        for (auto & asplit : ABpairs[i].first) {
            // Drop asplit from starttree
            starttree.dropSplit(asplit);
        }
        for (auto & bsplit : ABpairs[i].second) {
            // Add bsplit to starttree
            double edgelen = bsplit.getEdgeLen();
            edgelen *= edgelen_multiplicative_factor;
            starttree.addSplit(bsplit, edgelen);
        }
    }
}

inline bool OP::frechetCloseEnough(vector<TreeManip> & mu, unsigned lower, unsigned upper, double epsilon) const {
    // Compute pairwise distances between trees in mu with index > lower and index <= upper and return
    // true iff all pairwise distances are less than epsilon
    bool is_close_enough = true;
    vector<Split::treeid_pair_t> ABpairs;
    assert(lower < upper);
    assert (upper < mu.size());
    for (unsigned i = lower; i < upper - 1; ++i) {
        for (unsigned j = i+1; j < upper; ++j) {
            double bhvdist = calcBHVDistance(mu[i], mu[j], ABpairs);
            if (bhvdist > epsilon) {
                is_close_enough = false;
                break;
            }
        }
    }
    return is_close_enough;
}

inline void OP::computeFrechetMean(TreeManip & meantree) const {
    double   epsilon = 0.001; // successive mean estimates must be at least this close to stop iterating
    unsigned       N = 5;     // number of previous mean estimates that must be as close as epsilon
    unsigned       K = 100;   // maximum number of iterations
    unsigned k = 1;           // keeps track of iterations
    vector<TreeManip> mu(k);  // the trail of estimated mean trees (always has length k)
    Lot lot;
    lot.setSeed(_rnseed);
    chooseRandomTree(mu[k], lot);
    bool close_enough = false;
    while (k < K || close_enough) {
        mu.emplace_back();
        chooseRandomTree(mu[k+1], lot);
        displaceTreeAlongGeodesic(mu[k+1], mu[k], 1.0/(k+1));
        ++k;
        if (k >= N)
            close_enough = frechetCloseEnough(mu, k-N, k, epsilon);
        }
    meantree.setTree(mu[k].getTree());
}

inline void OP::run() {
    try {
        // Read in trees
        _tree_summary = std::make_shared<TreeSummary>();
        _tree_summary->readTreefile(_tree_file_name, 0);
        unsigned ntrees = _tree_summary->getNumTrees();
        if (ntrees < 2) {
            throw Xop("Must input at least 2 trees to compute tree distances");
        }

        //temporary!
        TreeManip starttm;
        buildTree(0, starttm);
        cerr << "start: " << starttm.makeNewick(5) << endl;
        TreeManip endtm;
        buildTree(1, endtm);
        cerr << "end: " << starttm.makeNewick(5) << endl;
        displaceTreeAlongGeodesic(starttm, endtm, 1.0);
        cerr << starttm.makeNewick(5) << endl;
        cerr << "start moved to end: " << starttm.makeNewick(5) << endl;
        exit(0);

        if (_output_for_gtp) {
            if (!_quiet)
                cout << "Writing trees in newick format to file \"trees-for-gtp.txt\"" << endl;

            ofstream gtpf("trees-for-gtp.txt");
            for (unsigned i = 0; i < ntrees; ++i) {
                gtpf << _tree_summary->getNewick(i) << ";\n";
            }
            gtpf.close();

            // If output for gtp was requested, then that is all
            // we try to do on this run
            if (!_quiet)
                cout << "Done." << endl;
            return;
        }

        if (_frechet_mean) {
            TreeManip meantree;
            if (!_quiet)
                cout << "Computing Frechet mean tree..." << endl;
            computeFrechetMean(meantree);

            if (!_quiet) {
                cout << "Mean tree:" << endl;
                cout << meantree.makeNewick(9) << endl;
            }

            if (!_quiet)
                cout << "Done." << endl;
            return;
        }

        if (_distance_measure == "geodesic") {
            if (!_quiet)
                cout << "Writing geodesic distances to file \"bhvdists.txt\"" << endl;
            ofstream outf("bhvdists.txt");
            outf << "tree	distance to tree 1" << endl;
            TreeManip starttm;
            buildTree(0, starttm);
            for (unsigned i = 1; i < ntrees; i++) {
                TreeManip endtm;
                buildTree(i, endtm);
                vector<Split::treeid_pair_t> ABpairs;
                double bhvdist = calcBHVDistance(starttm, endtm, ABpairs);
                outf << (i+1) << '\t' << setprecision(static_cast<int>(_precision)) << bhvdist << '\n';
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
        cerr << "OP encountered a problem:\n  " << x.what() << endl;
        }
}

} // namespace op

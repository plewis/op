#pragma once

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "split.hpp"
#include "tree_manip.hpp"
#include "xop.hpp"

#include "ncl/nxsmultiformat.h"

using namespace std;
using namespace boost;

namespace op
    {

    class TreeSummary
        {
        public:
                                        TreeSummary();
                                        ~TreeSummary();

            void                        readTreefile(const string filename, unsigned skip);
            void                        showSummary() const;
            unsigned                    getNumTrees() const;
            typename Tree::SharedPtr    getTree(unsigned index);
            string                 getNewick(unsigned index);
            bool                        isRooted(unsigned index);
            void                        clear();

        private:

            Split::treemap_t            _treeIDs;
            vector<string>    _newicks;
            vector<bool>           _is_rooted;

        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
        };

inline TreeSummary::TreeSummary()
    {
    //cout << "Constructing a TreeSummary" << endl;
    }

inline TreeSummary::~TreeSummary()
    {
    //cout << "Destroying a TreeSummary" << endl;
    }

inline unsigned TreeSummary::getNumTrees() const {
    assert(_newicks.size() > 0);
    return (unsigned)_newicks.size();
}

inline Tree::SharedPtr TreeSummary::getTree(unsigned index)
    {
    if (index >= _newicks.size())
        throw Xop("getTree called with index >= number of stored trees");

    TreeManip tm;

    // build the tree
    tm.buildFromNewick(_newicks[index], false, true);

    return tm.getTree();
    }

inline bool TreeSummary::isRooted(unsigned index) {
    if (index >= _is_rooted.size())
        throw Xop("isRooted called with index >= number of stored trees");

        return _is_rooted[index];
}

inline string TreeSummary::getNewick(unsigned index) {
    if (index >= _newicks.size())
        throw Xop("getNewick called with index >= number of stored trees");

        return _newicks[index];
}

inline void TreeSummary::clear()
    {
    _treeIDs.clear();
    _newicks.clear();
    }

inline void TreeSummary::readTreefile(const string filename, unsigned skip)
    {
    TreeManip tm;
    Split::treeid_t splitset;
    Split::treeid_t leafsplitset;

    // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

    //MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
    MultiFormatReader nexusReader(-1, NxsReader::IGNORE_WARNINGS);

    // Both of these needed to suppress "storing read block" messages
    // see NxsReader::statusMessage in nxsreader.cpp
    nexusReader.SetAlwaysReportStatusMessages(false);
    nexusReader.SetWarningOutputLevel(NxsReader::SUPPRESS_WARNINGS_LEVEL);

    try {
        nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
    catch(...)
        {
        nexusReader.DeleteBlocksFromFactories();
        throw;
        }

    int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
    for (int i = 0; i < numTaxaBlocks; ++i)
        {
        clear();
        NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
        string taxaBlockTitle = taxaBlock->GetTitle();

        const unsigned nTreesBlocks = nexusReader.GetNumTreesBlocks(taxaBlock);
        for (unsigned j = 0; j < nTreesBlocks; ++j)
            {
            const NxsTreesBlock * treesBlock = nexusReader.GetTreesBlock(taxaBlock, j);
            unsigned ntrees = treesBlock->GetNumTrees();
            if (skip < ntrees)
                {
                //cout << "Trees block contains " << ntrees << " tree descriptions.\n";
                for (unsigned t = skip; t < ntrees; ++t)
                    {
                    // NxsFullTreeDescription::TreeDescFlags possibilities
                    // (stored in flags data member):
                    // NXS_IS_ROOTED_BIT					= 0x0001,
                    // NXS_HAS_SOME_EDGE_LENGTHS_BIT		= 0x0002,
                    // NXS_MISSING_SOME_EDGE_LENGTHS_BIT	= 0x0004,
                    // NXS_EDGE_LENGTH_UNION 				= 0x0006,
                    // NXS_INT_EDGE_LENGTHS_BIT 			= 0x0008,
                    // NXS_HAS_ALL_TAXA_BIT				= 0x0010,
                    // NXS_HAS_NHX_BIT 					= 0x0020,
                    // NXS_HAS_DEG_TWO_NODES_BIT			= 0x0040,
                    // NXS_HAS_POLYTOMY_BIT				= 0x0080,
                    // NXS_HAS_INTERNAL_NAMES_BIT			= 0x0100,
                    // NXS_HAS_NEW_INTERNAL_NAMES_BIT		= 0x0200,
                    // NXS_KNOWN_INTERNAL_NAMES_BIT		= 0x0400,
                    // NXS_SOME_ZERO_EDGE_LEN_BIT			= 0x0800,
                    // NXS_SOME_NEGATIVE_EDGE_LEN_BIT		= 0x1000,
                    // NXS_TREE_PROCESSED 					= 0x2000
                    
                    const NxsFullTreeDescription & d = treesBlock->GetFullTreeDescription(t);
                    
                    bool is_rooted = d.IsRooted();
                    _is_rooted.push_back(is_rooted);
                    
                    // store the newick tree description
                    string newick = d.GetNewick();;
                    _newicks.push_back(newick);
                    unsigned tree_index = (unsigned)_newicks.size() - 1;

                    // build the tree
                    tm.buildFromNewick(newick, /*rooted*/is_rooted, /*allow_polytomies*/true);

                    // store set of splits
                    splitset.clear();
                    leafsplitset.clear();
                    tm.storeSplits(splitset, leafsplitset);

                    // iterator iter will point to the value corresponding to key splitset
                    // or to end (if splitset is not already a key in the map)
                    Split::treemap_t::iterator iter = _treeIDs.lower_bound(splitset);

                    if (iter == _treeIDs.end() || iter->first != splitset)
                        {
                        // splitset key not found in map, need to create an entry
                        vector<unsigned> v(1, tree_index);  // vector of length 1 with only element set to tree_index
                        _treeIDs.insert(iter, Split::treemap_t::value_type(splitset, v));
                        }
                    else
                        {
                        // splitset key was found in map, need to add this tree's index to vector
                        iter->second.push_back(tree_index);
                        }
                    } // trees loop
                } // if skip < ntrees
            } // TREES block loop
        } // TAXA block loop

    // No longer any need to store raw data from nexus file
    nexusReader.DeleteBlocksFromFactories();
    }

inline void TreeSummary::showSummary() const
    {
    // Produce some output to show that it works
    cout << str(format("\nRead %d trees from file") % _newicks.size()) << endl;

    // Show all unique topologies with a list of the trees that have that topology
    // Also create a map that can be used to sort topologies by their sample frequency
    typedef pair<unsigned, unsigned> sorted_pair_t;
    vector< sorted_pair_t > sorted;
    int t = 0;
    for (auto & key_value_pair : _treeIDs)
        {
        unsigned topology = ++t;
        unsigned ntrees = (unsigned)key_value_pair.second.size();
        sorted.push_back(pair<unsigned, unsigned>(ntrees,topology));
        cout << "Topology " << topology << " seen in these " << ntrees << " trees:" << endl << "  ";
        copy(key_value_pair.second.begin(), key_value_pair.second.end(), ostream_iterator<unsigned>(cout, " "));
        cout << endl;
        }

    // Show sorted histogram data
    sort(sorted.begin(), sorted.end());
    //unsigned npairs = (unsigned)sorted.size();
    cout << "\nTopologies sorted by sample frequency:" << endl;
    cout << str(format("%20s %20s") % "topology" % "frequency") << endl;
    for (auto & ntrees_topol_pair : adaptors::reverse(sorted))
        {
        unsigned n = ntrees_topol_pair.first;
        unsigned t = ntrees_topol_pair.second;
        cout << str(format("%20d %20d") % t % n) << endl;
        }
    }

    }

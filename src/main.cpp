// Uncomment the line below to save incompatibility graphs as dot files for debugging purposes
//#define OP_SAVE_DOT_FILE

#include <string>
#include <vector>
#include <iostream>
#include <numeric>

using namespace std;

#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/exception/all.hpp> // Include for boost::exception
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "split.hpp"
#include "node.hpp"
#include "tree.hpp"
#include "tree_manip.hpp"
#include "tree_summary.hpp"
#include "opvertex.hpp"
#include "xop.hpp"
#include "lot.hpp"
#include "op.hpp"

using namespace op;

// static data member initializations
string  OP::_program_name        = "kfdist";
unsigned     OP::_major_version       = 1;
unsigned     OP::_minor_version       = 0;
const double Node::_smallest_edge_length = 1.0e-12;
vector<string> TreeManip::_taxon_names;
map<string, unsigned> TreeManip::_taxon_map;

#if defined(OP_SAVE_DOT_FILE)
unsigned OP::_graph_number          = 1;
#endif

int main(int argc, const char * argv[]) {
    try {
        OP strom;
        strom.processCommandLineOptions(argc, argv);
        strom.run();
    }
    catch(std::exception & x) {
        cerr << "Exception: " << x.what() << endl;
        cerr << "Aborted." << endl;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

    return 0;
}

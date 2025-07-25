#include <iostream>

// Uncomment the line below to save incompatibility graphs as dot files for debugging purposes
//#define OP_SAVE_DOT_FILE

#include "op.hpp"

using namespace op;

// static data member initializations
string  OP::_program_name        = "kfdist";
unsigned     OP::_major_version       = 1;
unsigned     OP::_minor_version       = 0;
const double Node::_smallest_edge_length = 1.0e-12;

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

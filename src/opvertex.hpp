#pragma once

#include <boost/format.hpp>

#include "split.hpp"

using namespace std;
using namespace boost;

namespace op {

    struct OPVertex {
        OPVertex() : _split(nullptr), _capacity(0), _parent_index(-1) {}
        const Split * _split;
        double _capacity;
        int _parent_index;
        vector<OPVertex *> _children;
    };

    typedef map<pair<OPVertex *, OPVertex *>, double> edgemap_t;
}
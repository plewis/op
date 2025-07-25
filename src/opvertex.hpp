#pragma once

#include <sstream>
#include <boost/format.hpp>

#include "split.hpp"

using namespace std;
using namespace boost;

namespace op {
    string memoryAddressAsString(const void * ptr) {
        if (ptr == nullptr) {
            return "nullptr";
        }
        ostringstream memory_address;
        memory_address << ptr;
        return memory_address.str();
    }

    struct OPEdge;

    struct OPVertex {
        OPVertex() : _split(nullptr), _weight(0), _residual_capacity(0.0), _parent_edge(nullptr) {
            //cerr << "in OPVertex default constructor" << endl;
            //cerr << "  _split:        " << memoryAddressAsString(_split) << endl;
            //cerr << "  _weight:       " << _weight << endl;
            //cerr << "  _parent_edge:  " << memoryAddressAsString(_parent_edge) << endl;
            //cerr << "  _edges.size(): " << _edges.size() << endl;
        }
        OPVertex(const OPVertex & other) : _split(nullptr), _weight(0), _residual_capacity(0.0), _parent_edge(nullptr) {
            _name = other._name;
            _split = other._split;
            _weight = other._weight;
            _residual_capacity = other._residual_capacity;
            _parent_edge = other._parent_edge;
            _edges = other._edges;
            //cerr << "in OPVertex copy constructor" << endl;
            //cerr << "  other._split:        " << memoryAddressAsString(other._split) << endl;
            //cerr << "  other._weight:       " << other._weight << endl;
            //cerr << "  other._parent_edge:  " << memoryAddressAsString(other._parent_edge) << endl;
            //cerr << "  other._edges.size(): " << other._edges.size() << endl;
        }
        string _name;
        const Split * _split;
        double _weight;
        double _residual_capacity;
        OPEdge * _parent_edge;
        vector<OPEdge *> _edges;
    };

    struct OPEdge {
        OPEdge() : _from(nullptr), _to(nullptr), _open(false), _capacity(0), _flow(0), _reverse_flow(0.0), _edge_is_reversed(false) {
            //cerr << "in OPEdge default constructor:" << endl;
            //cerr << "  _from:         " << memoryAddressAsString(_from) << endl;
            //cerr << "  _to:           " << memoryAddressAsString(_to) << endl;
            //cerr << "  _capacity:     " << _capacity << endl;
            //cerr << "  _flow:         " << _flow << endl;
            //cerr << "  _reverse_flow: " << _reverse_flow << endl;
        }
        OPEdge(const OPEdge & other) : _from(nullptr), _to(nullptr), _capacity(0), _reverse_flow(0.0), _edge_is_reversed(false) {
            _from = other._from;
            _to = other._to;
            _capacity = other._capacity;
            _reverse_flow = other._reverse_flow;
            _edge_is_reversed = other._edge_is_reversed;
            //cerr << "in OPEdge copy constructor:" << endl;
            //cerr << "  other._from:         " << memoryAddressAsString(other._from) << endl;
            //cerr << "  other._to:           " << memoryAddressAsString(other._to) << endl;
            //cerr << "  other._capacity:     " << other._capacity << endl;
            //cerr << "  other._flow:         " << other._flow << endl;
            //cerr << "  other._reverse_flow: " << other._reverse_flow << endl;
        }
        OPVertex * _from;
        OPVertex * _to;
        bool _open;
        double _capacity;
        double _flow;
        double _reverse_flow;
        bool _edge_is_reversed;
    };

}
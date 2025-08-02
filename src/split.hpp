#pragma once

#include <vector>
#include <memory>
#include <set>
#include <map>
#include <climits>
#include <cassert>
#include <boost/container/container_fwd.hpp>

using namespace std;
using namespace boost;

namespace op
    {

    class Split
        {
        public:
                                                                Split();
                                                                ~Split();

            bool                                                operator==(const Split & other) const;
            bool                                                operator<(const Split & other) const;

            void                                                clear();
            void                                                resize(unsigned nleaves);

            typedef unsigned long                               split_unit_t;
            typedef vector<split_unit_t>                        split_t;
            typedef set<Split>                                  treeid_t;
            typedef pair<treeid_t, treeid_t>                    treeid_pair_t;
            typedef pair<Split, Split>                          split_pair_t;
            typedef map< treeid_t, vector<unsigned> >           treemap_t;

            void                                                setEdgeLen(double v);
            double                                              getEdgeLen() const;

            const split_t &                                     getBits() const;

            unsigned                                            getSize() const;
            unsigned                                            getNumBitsSet() const;

            void                                                bitwiseAnd(const Split & other);
            unsigned                                            findFirstSetBit() const;
            void                                                invertBits();
            bool                                                subsumedIn(const Split & other) const;
            bool                                                compatibleWith(const Split & other) const;

            void                                                setBitAt(unsigned leaf_index);
            void                                                addSplit(const Split & other);

            string                                              createPatternRepresentation(bool show_edge_length = false) const;

        private:

            double                                              _edgelen;
            split_t                                             _bits;
            unsigned                                            _bits_per_unit;
            unsigned                                            _nleaves;

        public:

            typedef std::shared_ptr< Split >                    SharedPtr;
    };

inline Split::Split() {
    _edgelen = 0.0;
    _nleaves = 0;
    _bits_per_unit = (CHAR_BIT)*sizeof(Split::split_unit_t);
    clear();
    //cout << "Constructing a Split" << endl;
    }

inline Split::~Split() {
    //cout << "Destroying a Split" << endl;
    }

inline void Split::clear() {
    _edgelen = 0.0;
    for (auto & u : _bits) {
        u = 0L;
        }
    }

inline bool Split::operator==(const Split & other) const {
    return (_bits == other._bits);
    }

inline bool Split::operator<(const Split & other) const {
    assert(_bits.size() == other._bits.size());
    return (_bits < other._bits);
    }

inline const Split::split_t & Split::getBits() const {
    return _bits;
}

inline unsigned Split::getSize() const {
    return _nleaves;
}

inline void Split::resize(unsigned nleaves) {
    _nleaves = nleaves;
    unsigned nunits = 1 + ((nleaves - 1)/_bits_per_unit);
    _bits.resize(nunits);
    clear();
}

    inline void Split::setEdgeLen(double edgelen) {
    _edgelen = edgelen;
    }

inline double Split::getEdgeLen() const {
    return _edgelen;
}

inline void Split::setBitAt(unsigned leaf_index) {
    unsigned unit_index = leaf_index/_bits_per_unit;
    unsigned bit_index = leaf_index - unit_index*_bits_per_unit;
    split_unit_t bit_to_set = (split_unit_t)1 << bit_index;
    _bits[unit_index] |= bit_to_set;
    }

inline void Split::addSplit(const Split & other) {
    unsigned nunits = (unsigned)_bits.size();
    assert(nunits == other._bits.size());
    for (unsigned i = 0; i < nunits; ++i)
        {
        _bits[i] |= other._bits[i];
        }
    }

inline unsigned Split::getNumBitsSet() const {
    unsigned n = 0;
    unsigned ntax_added = 0;
    for (unsigned i = 0; i < _bits.size(); ++i) {
        for (unsigned j = 0; j < _bits_per_unit; ++j) {
            split_unit_t bitmask = ((split_unit_t)1 << j);
            bool bit_is_set = ((_bits[i] & bitmask) > (split_unit_t)0);
            if (bit_is_set)
                n++;
            if (++ntax_added == _nleaves)
                break;
        }
    }
    return n;
}

inline void Split::bitwiseAnd(const Split & other) {
    const split_t & other_bits = other.getBits();
    assert(_bits.size() == other_bits.size());
    for (unsigned i = 0; i < _bits.size(); ++i) {
        _bits[i] &= other_bits[i];
    }
}

inline unsigned Split::findFirstSetBit() const {
    unsigned bit_index = 0;
    for (auto & u : _bits) {
        if (u > 0) {
            while ((u & ((split_unit_t)1 << bit_index)) == 0)
                bit_index++;
            break;
        }
    }
    return bit_index;
}

inline void Split::invertBits() {
    for (auto & u : _bits) {
        u = ~u;
    }
}

inline bool Split::compatibleWith(const Split & other) const {
    const split_t & other_bits = other._bits;
    assert(_bits.size() == other_bits.size());
    bool is_compatible = true;
    for (unsigned i = 0; i < _bits.size(); ++i) {
        split_unit_t a       = _bits[i];
        split_unit_t b       = other_bits[i];
        split_unit_t a_and_b = (a & b);
        bool equals_a   = (a_and_b == a);
        bool equals_b   = (a_and_b == b);
        if (a_and_b && !(equals_a || equals_b)) {
            // A failure of any unit to be compatible makes the entire split incompatible
            is_compatible = false;
            break;
        }
    }
    return is_compatible;
}

inline bool Split::subsumedIn(const Split & other) const {
    // Example: this is NOT subsumed in other
    // 8421
    // **-* 13 this <----------+
    // **-- 12 other           | not equal
    // **-- 12 this & other <--+
    //
    // Example: this IS subsumed in other
    // 8421
    // **-- 12 this <----------+
    // **-* 13 other           | equal
    // **-- 12 this & other <--+
    //
    // Example: this is compatible with, but not subsumed in, other
    // 8421
    // **-- 12 this <----------+
    // --**  3 other           | not equal
    // ----  0 this & other <--+
    const split_t & other_bits = other.getBits();
    assert(_bits.size() == other_bits.size());
    for (unsigned i = 0; i < _bits.size(); ++i) {
        if ((_bits[i] & other_bits[i]) != _bits[i])
            return false;
    }
    return true;
}

inline string Split::createPatternRepresentation(bool show_edge_length) const {
    string s;
    if (show_edge_length)
        s += str(format("%.3f: ") % _edgelen);
    unsigned ntax_added = 0;
    for (unsigned i = 0; i < _bits.size(); ++i) {
        for (unsigned j = 0; j < _bits_per_unit; ++j) {
            split_unit_t bitmask = ((split_unit_t)1 << j);
            bool bit_is_set = ((_bits[i] & bitmask) > (split_unit_t)0);
            if (bit_is_set)
                s += '*';
            else
                s += '-';
            if (++ntax_added == _nleaves)
                break;
            }
        }
    return s;
    }

}

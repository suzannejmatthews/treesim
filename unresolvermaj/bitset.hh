// Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
//                             University of Aarhus, Denmark
//                             Contact: Thomas Mailund <mailund@birc.dk>

#ifndef BITSET_HH
#define BITSET_HH

#include <climits>
#include <cstdlib>
#include <iostream>


class BitSet {
    size_t  _size;
    size_t *_data;

    // sulsj
    int _ones;
    int _zeros;

    static const size_t WORD_BITS = sizeof(size_t)*CHAR_BIT;

    static size_t word_index(size_t i) {
        return i/WORD_BITS;
    }
    static size_t bit_index (size_t i) {
        return 1 << (WORD_BITS-1-(i%WORD_BITS));
    }

    size_t array_length() const {
        return _size/WORD_BITS + 1;
    }

public:
    BitSet() {};
    BitSet(size_t size);
    BitSet(const BitSet &b);
    ~BitSet();

    BitSet &operator=(const BitSet &b);

    size_t size() const {
        return _size;
    }

    class Bit {
        BitSet &_bs;
        size_t  _idx;
    public:
        Bit(BitSet &bs, size_t idx);
        ~Bit();

        operator bool() const;
        bool operator~() const;
        Bit& operator=(bool b);
        Bit& operator=(const Bit &b);
    }; // class Bit

    Bit  at(size_t idx) {
        return Bit(*this,idx);
    }
    bool at(size_t idx) const {
        return _data[word_index(idx)] & bit_index(idx);
    }

    Bit  operator[](size_t idx)       {
        return at(idx);
    }
    bool operator[](size_t idx) const {
        return at(idx);
    }

    bool operator==(const BitSet &b) const;
    bool operator!=(const BitSet &b) const {
        return !(*this == b);
    }
    bool operator< (const BitSet &b) const;

    BitSet &operator|=(const BitSet &b);
    BitSet &operator&=(const BitSet &b);
    BitSet &operator^=(const BitSet &b);

    void print(std::ostream &) const;

    // ssj
//  char* bs2string(BitSet &b);	// original
    char* bs2string();			// overloaded by ssj
    int num_ones()  const {
        return _ones;
    }
    int num_zeros() const {
        return _zeros;
    }
//	bool IsTrivial() const { return _ones == 1 or _zeros == 1; }

    void inc_ones() 	{
        ++_ones;
    }
//	bool IsTrivial();
};




inline BitSet::Bit::Bit(BitSet &bs, size_t idx)
    : _bs(bs), _idx(idx)
{}
inline BitSet::Bit::~Bit()
{}

inline BitSet::Bit::operator bool() const
{
    return (_bs._data[BitSet::word_index(_idx)] & BitSet::bit_index(_idx));
}

inline bool
BitSet::Bit::operator~() const
{
    return ~(_bs._data[BitSet::word_index(_idx)] & BitSet::bit_index(_idx));
}

inline BitSet::Bit&
BitSet::Bit::operator=(bool b)
{
    if (b) _bs._data[BitSet::word_index(_idx)] |= BitSet::bit_index(_idx);
    else   _bs._data[BitSet::word_index(_idx)] &= ~BitSet::bit_index(_idx);
    return *this;
}

inline BitSet::Bit&
BitSet::Bit::operator=(const Bit &b)
{
    if (bool(b)) _bs._data[BitSet::word_index(_idx)] |= BitSet::bit_index(_idx);
    else         _bs._data[BitSet::word_index(_idx)] &= ~BitSet::bit_index(_idx);
    return *this;
}




inline BitSet operator|(const BitSet &b1, const BitSet &b2)
{
    BitSet res(b1);
    res |= b2;
    return res;
}

inline BitSet operator&(const BitSet &b1, const BitSet &b2)
{
    BitSet res(b1);
    res &= b2;
    return res;
}




inline std::ostream &operator<<(std::ostream &os, const BitSet &bs)
{
    bs.print(os);
    return os;
}


#endif // BITSET_HH

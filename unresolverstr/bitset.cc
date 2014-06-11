// Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
//                             University of Aarhus, Denmark
//                             Contact: Thomas Mailund <mailund@birc.dk>

#include "bitset.hh"

#include <cassert>
#include <algorithm>
#include <functional>
using namespace std;

// -- CONSTRUCTION AND DESTRUCTION --------------------------------------
BitSet::BitSet(size_t size) : _size(size), _ones(0), _zeros(0)
{
    _data = new size_t[array_length()];
    fill(_data, _data+array_length(), 0);

////	for ( unsigned int i = 0; i < _size; ++i )
////		if ( at(i) ) 	++_ones;
////		else 	       ++_zeros;
}

BitSet::BitSet(const BitSet &b)
{
    _size = b._size;
    _data = new size_t[array_length()];
    copy(b._data,b._data+array_length(),_data);
}

BitSet::~BitSet()
{
    delete[] _data;
}

BitSet &
BitSet::operator=(const BitSet &b)
{
    if (&b != this) {
        delete[] _data;
        _size = b._size;
        _data = new size_t[array_length()];
        copy(b._data,b._data+array_length(),_data);
    }
    return *this;
}


// -- COMPARISONS -------------------------------------------------------
bool
BitSet::operator==(const BitSet &b) const
{
    if (_size != b._size) return false;
    return equal(_data, _data+array_length(), b._data);
}

bool
BitSet::operator< (const BitSet &b) const
{
    return lexicographical_compare(_data, _data+array_length(),
                                   b._data, b._data+b.array_length());
}


// -- SET OPERATIONS ----------------------------------------------------
namespace {
// stl only defines _logical_ && and ||, we need _bit_ & and |.
template<typename T> struct binary_or : public binary_function<T,T,T> {
    explicit binary_or() {};
    T operator()(T t1, T t2) const {
        return t1 | t2;
    }
};
template<typename T> struct binary_and : public binary_function<T,T,T> {
    explicit binary_and() {};
    T operator()(T t1, T t2) const {
        return t1 & t2;
    }
};
template<typename T> struct binary_xor : public binary_function<T,T,T> {
    explicit binary_xor() {};
    T operator()(T t1, T t2) const {
        return t1 ^ t2;
    }
};
};

BitSet &
BitSet::operator|=(const BitSet &b)
{
    assert(_size == b._size);
    transform(_data,_data+array_length(),b._data, _data, binary_or<size_t>());
    return *this;
}

BitSet &
BitSet::operator&=(const BitSet &b)
{
    assert(_size == b._size);
    transform(_data,_data+array_length(),b._data, _data, binary_and<size_t>());
    return *this;
}

BitSet &
BitSet::operator^=(const BitSet &b)
{
    assert(_size == b._size);
    transform(_data,_data+array_length(),b._data, _data, binary_xor<size_t>());
    return *this;
}


// --- for debugging, mainly ------------------------------------------
void
BitSet::print(std::ostream &os) const
{
    os << '[';
    for (size_t i = 0; i < _size; ++i) os << (*this)[i];
    os << ']';
}

//int
//BitSet::count_ones()
//{
//	int ones=0;
//	for ( int i = 0; i < size(); ++i )
//		if ( at(i) ) 	++ones;
//
//	return ones;
//}
//
//int
//BitSet::count_zeros()
//{
//	int zeros=0;
//	for ( int i = 0; i < size(); ++i )
//		if ( !at(i) ) 	++zeros;
//
//	return zeros;
//}

void
BitSet::count_ones_zeros()
{
    for ( unsigned int i = 0; i < _size; ++i )
        ( at(i) ) ? ++_ones : ++_zeros;
}

//bool
//BitSet::IsTrivial()
//{
//	int ones=0;
//	int zeros=0;
//
//	for ( unsigned int i = 0; i < _size; ++i ) {
//		if ( ones > 1 && zeros > 1 ) return false;	// bs is not trivial
//		( at(i) ) ? ++ones : ++zeros;
//	}
//
//	return true;	// trivial!
//}


// sulsj
char*
BitSet::bs2string()
{
    char* ch = new char[_size+1];
    {
        assert(ch != NULL);
    }
    for (size_t i = 0; i < _size; ++i) {
//		((*this)[i] == true) ? ch[i] = '1' : ch[i] = '0';
        if ((*this)[i] == true) {
            ch[i] = '1';
            ++_ones;
        } else {
            ch[i] = '0';
            ++_zeros;
        }
    }

    ch[_size] = '\0';
    return ch;
}





///*****************************************************/
//
// Copyright (C) 2006, 2007 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS DEFINITION
//		CHashFunc: Universal hash functions
//
///*****************************************************/

#ifndef HASHFUNC_HH
#define HASHFUNC_HH

#include <iostream>
#include "bitset.hh"


typedef struct {
    unsigned long hv1;
    unsigned long hv2;
} HV_STRUCT_T;


class CHashFunc {

    unsigned long   	_m1;  // prime number1 for hash function1
    unsigned long   	_m2;  // prime number1 for hash function2
    unsigned  		_t;  	// number of trees
    unsigned  		_n;  	// number of taxa
    unsigned long * 	_a1;  // random numbers for hash function1
    unsigned long * 	_a2;  // random numbers for hash function2
    float   		_r;  	// hash table factor
    unsigned long 	_c;  	// double collision factor: constant for c*t*n of hash function2;

public:
    CHashFunc() : _m1(0), _m2(0), _t(0), _n(0), _a1(NULL), _a2(NULL), _r(0.0), _c(0) {}
    CHashFunc(unsigned t, unsigned n, float r, unsigned long c);
    ~CHashFunc();

    void UHashfunc_init(unsigned t, unsigned n, float r, unsigned long c);

//	void UHashFunc(HV_STRUCT_T &hv, BitSet &bs, unsigned num_taxa);
//	void UHashFunc(HV_STRUCT_T &hv, const std::string &c, unsigned num_taxa); // with type3
    void UHashFunc(HV_STRUCT_T &hv, unsigned long long bs, unsigned num_taxa); // with type3

    // Implicit bp
    unsigned long getA1(unsigned idx) {
        return (_a1[idx]);
    }
    unsigned long getA2(unsigned idx) {
        return (_a2[idx]);
    }
    unsigned long getM1() {
        return _m1;
    }
    unsigned long getM2() {
        return _m2;
    }

    unsigned long GetPrime(unsigned long top_num, unsigned from);
};

#endif

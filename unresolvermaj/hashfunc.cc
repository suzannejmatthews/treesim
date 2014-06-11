///*****************************************************/
//
// Copyright (C) 2006, 2007 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS IMPLEMENTATION
//		CHashFunc: Class for Universal hash functions
//
///*****************************************************/

#include "hashfunc.hh"
#include "randomc.h" // MT random number generator

#include <iostream>
#include <cassert>

// constructor
CHashFunc::CHashFunc(
    unsigned t,
    unsigned n,
    float r,
    unsigned long c)
    : _m1(0), _m2(0), _t(t), _n(n), _a1(NULL), _a2(NULL), _r(0.0), _c(c)
{
    UHashfunc_init(t, n, r, c);
}

// destructor
CHashFunc::~CHashFunc()
{
    delete[] _a1;
    delete[] _a2;
}

void
CHashFunc::UHashfunc_init(
    unsigned t,
    unsigned n,
    float r,
    unsigned long c)
{
    // Init member variables
    _t = t;
    _n = n;
    _r = r;
    _c = c;

    // Get the prime number which is larger than t*n of the amount of
    // t*n*r. If t=100, n=100 and r=0.1 then the size of hash table
    // will be 10000+1000 = 11000
    unsigned long top = _t*_n + int(_t*_n*_r);
//	unsigned long p = GetPrime(top);
    unsigned long p =0;
    unsigned mul = 1;
    do {
        unsigned from = 100 * mul;
        p = GetPrime(top, from);
        ++mul;
    } while (p == 0);

    assert(top != 0);
    assert(p != 0);

    _m1 = p;   	// t*n ~~ m1

//	unsigned long top2 = _c*_t*_n;
//	unsigned long p2 = GetPrime(top2);
    unsigned long top2 = _c*top;
    unsigned long p2 = 0;
    mul = 1;
    do {
        unsigned from = 100 * mul;
        p2 = GetPrime(top2, from);
        ++mul;
    } while (p2 == 0);

    assert(top2 != 0);
    assert(p2 != 0);

    _m2 = p2; // m2 > c*t*n --> to avoid double collision ==> I just use _c*p for _m2

    int32 seed = time(0);
    TRandomMersenne rg(seed);

    // universal hash funciton = a * b mod m
    // where B=(b1, b2,...,bn) is bit-string
    _a1 = new unsigned long [_n];
    _a2 = new unsigned long [_n];


    // generate n random numbers between 0 and m1-1
    // for hash function1 and hash function2
    // rand() % 48
    // random number between 0~47
    for (unsigned int i=0; i<_n; i++) {
        _a1[i] = rg.IRandom(0,_m1-1);
        _a2[i] = rg.IRandom(0,_m2-1);
    }
}


void
CHashFunc::UHashFunc(
    HV_STRUCT_T &hv,
    unsigned long long bs,
    unsigned numBits)
{
//	unsigned len = numBits;
    hv.hv1=0;
    hv.hv2=0;
//	assert(num_taxa == len);

    for (unsigned i=0; i<numBits; i++) {
        if (bs&(1<<i)) {
//			hv.hv1 += _a1[i]*ii;
//			hv.hv2 += _a2[i]*ii;
            hv.hv1 += _a1[i];
            hv.hv2 += _a2[i];
        }
    }
//	assert(hv.hv1 != 0);
//	assert(hv.hv2 != 0);
    hv.hv1 %= _m1;
    hv.hv2 %= _m2;
}


//void
//CHashFunc::UHashFunc(
//	HV_STRUCT_T &hv,
//	BitSet &bs,
//	unsigned numTaxa)
//{
////	int len = bs.size();
//	hv.hv1=0;
//	hv.hv2=0;
////	assert(numTaxa == len);
//	for (unsigned i=0; i<numTaxa; ++i) {
//		if (bs[i]) {
////			bs.inc_ones(); // increase the number of '1's.
//			// compute the two hash values using a1, a2 and bs.
////			hv.hv1 += _a1[i]*bs[i];
////			hv.hv2 += _a2[i]*bs[i];
//			hv.hv1 += _a1[i];
//			hv.hv2 += _a2[i];
//		}
//	}
////	assert(hv.hv1 != 0);
////	assert(hv.hv2 != 0);
//	hv.hv1 %= _m1;
//	hv.hv2 %= _m2;
//}


//// with type3
//void
//CHashFunc::UHashFunc(
//	HV_STRUCT_T &hv,
//	const std::string &s,
//	unsigned numTaxa)
//{
////	int len = s.length();
////	assert(len == numTaxa);
////	const char* d = s.c_str();
////	unsigned long h=0;
//	hv.hv1=0;
//	hv.hv2=0;
//	for (unsigned i=0; i<numTaxa; ++i) {
//		if (s[i] != '0') {
////	for (int i=0; i<numTaxa; ++i, ++d) {
////		hv.hv1 += _a1[i]*(*d);
////		hv.hv2 += _a2[i]*(*d);
////			hv.hv1 += _a1[i]*(s[i]);
////			hv.hv2 += _a2[i]*(s[i]);
//			hv.hv1 += _a1[i];
//			hv.hv2 += _a2[i];
//		}
//	}
////	h = (h+_b2) % _m2;
//	assert(hv.hv1 != 0);
//	assert(hv.hv2 != 0);
//	hv.hv1 %= _m1;
//	hv.hv2 %= _m2;
//}

// Generate a prime number right after topNum
unsigned long
CHashFunc::GetPrime(unsigned long topNum, unsigned from)
{
    unsigned long primeNum=0;
    unsigned long candidate=0;

    if (topNum <= 100)
        candidate = 2;
    else
        candidate = topNum - from;
//    std::cout << "candi=" << candidate << " topNum=" << topNum << std:: endl;

    while(candidate <= topNum) {
        unsigned long trialDivisor = 2;
        int prime = 1;

        while(trialDivisor * trialDivisor <= candidate) {
            if(candidate % trialDivisor == 0) {
                prime = 0;
                break;
            }
            trialDivisor++;
        }
        if(prime) primeNum = candidate;

        candidate++;
    }

//	std::cout << "PRIME=" << primeNum << std::endl;
    return primeNum;
}




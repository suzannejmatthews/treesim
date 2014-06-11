///*****************************************************/
//
// Copyright (C) 2006, 2007 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS DEFINITION
//		HashMap: Class for hashing
//
///*****************************************************/

#ifndef HASHMAP_HH
#define HASHMAP_HH

#include <vector>
using namespace std;

#include "bitset.hh"
#include "hashfunc.hh"

// For HashConsense
typedef struct {
    unsigned long	long	_hv2;
    BitSet* 						_bs;
    unsigned long				_count;
} BUCKET_STRUCT_T_HASHCS;

//HashConsense
typedef vector<BUCKET_STRUCT_T_HASHCS> V_BUCKET_T_HASHCS;

typedef struct {
    V_BUCKET_T_HASHCS 	_v_bucket;
    unsigned long				_max_count;
} HASHTAB_T2;

typedef vector<V_BUCKET_T_HASHCS> HASHTAB_T_HASHCS;
typedef vector<HASHTAB_T2> HASHTAB_T_HASHCS2;


//typedef HASHTAB_T_HASHCS::const_iterator HASHTAB_HASHCS_ITER_T;

class HashMap {

public:
    HashMap() {}
    ~HashMap() {}

    HASHTAB_T_HASHCS 		_hashtab_hashCS;
//	HASHTAB_T_HASHCS2 	_hashtab_hashCS2;

    CHashFunc 					_HF;

//	void hashing_bs_nbits(BitSet& bs, unsigned treeIdx, unsigned numTaxa, unsigned numTrees, unsigned op_mode, unsigned long hv1,	unsigned long hv2); // fast hash-rf
//	void hashing_bs_nbits_wo_T2(BitSet& bs, unsigned treeIdx, unsigned numTaxa, unsigned numTrees, unsigned op_mode, unsigned long hv1,	unsigned long hv2, vector<BitSet*> & vec_bs); // fast hash-rf
    void hashing_bs_nbits_wo_T2(BitSet& bs, unsigned treeIdx, unsigned numTaxa, unsigned numTrees, unsigned long hv1,	unsigned long hv2, vector<BitSet*> & vec_bs); // fast hash-rf
//	void hashing_bs_64bits(uint64_t bs64, BitSet& bs, unsigned numBits, unsigned treeIdx, unsigned numTrees, unsigned numTaxa, unsigned op_mode);

    void uhashfunc_init(unsigned t, unsigned n, float r, unsigned long c);
    void HashMap_clear();
};


#endif



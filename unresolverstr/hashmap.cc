///*****************************************************/
//
// Copyright (C) 2006, 2007 Seung-Jin Sul
// 		Department of Computer Science
// 		Texas A&M University
// 		Contact: sulsj@cs.tamu.edu
//
// 		CLASS IMPLEMENTATION
//		HashMap: Class for hashing bit
//
///*****************************************************/

#include "hashmap.hh"
#include "hashfunc.hh"

//#include <cassert>


// Implicit BPs ///////////////
void
HashMap::hashing_bs_nbits_wo_T2(
    BitSet&  bs,
    unsigned treeIdx,
    unsigned numTaxa,
    unsigned numTrees,
//	unsigned op_mode,
    unsigned long hv1,
    unsigned long hv2,
    vector<BitSet*> & vec_bs)
{
    if (treeIdx == 0) { // The first tree: just populate the bipaartitions.
        BUCKET_STRUCT_T_HASHCS bk;
        bk._hv2 = hv2;
        BitSet* tempbs = new BitSet(numTaxa);
        *tempbs = bs;
        bk._bs = tempbs;
        bk._count = 1;
        _hashtab_hashCS[hv1].push_back(bk);
//		_hashtab_hashCS2[hv1]._v_bucket.push_back(bk);
//		_hashtab_hashCS2[hv1]._max_count = 1;
    } else if (treeIdx > 0 & treeIdx < numTrees-1) { // The other trees: just compare the biaprtitions in the hashtable
//	else {
        if (!_hashtab_hashCS[hv1].empty()) { // if the hash table location has already been occupied
//		if ((!_hashtab_hashCS2[hv1]._v_bucket.empty()) & (_hashtab_hashCS2[hv1]._max_count == treeIdx)) {

            for (unsigned i=0; i<_hashtab_hashCS[hv1].size(); ++i) { // get the number of items in the linked list
                if (_hashtab_hashCS[hv1][i]._hv2 == hv2) { // if the same bp is found
                    _hashtab_hashCS[hv1][i]._count++;
//					if (_hashtab_hashCS[hv1][i]._count == numTrees)
//						vec_bs.push_back(_hashtab_hashCS[hv1][i]._bs);
                    break;
                }
//				else {
//					if (*(_hashtab_hashCS[hv1][i]._bs) == bs) {
//						cout << "Data type overflow\n";
//						exit(0);
//					}
//				}
            }
        }
    }

////////			for (unsigned i=0; i<_hashtab_hashCS2[hv1]._v_bucket.size(); ++i) { // get the number of items in the linked list
////////				if (_hashtab_hashCS2[hv1]._v_bucket[i]._hv2 == hv2) { // if the same bp is found
////////					_hashtab_hashCS2[hv1]._v_bucket[i]._count++;
////////					_hashtab_hashCS2[hv1]._max_count++;
////////					if (_hashtab_hashCS2[hv1]._v_bucket[i]._count == numTrees)
////////						vec_bs.push_back(_hashtab_hashCS2[hv1]._v_bucket[i]._bs);
////////					break;
////////				}
////////			}
////////
////////		}
////////	}
    else { // if treeIdx == numTrees
        if (!_hashtab_hashCS[hv1].empty()) {
            for (unsigned i=0; i<_hashtab_hashCS[hv1].size(); ++i) {
                if (_hashtab_hashCS[hv1][i]._hv2 == hv2) {
//					_hashtab_hashCS[hv1][i]._count++;
                    if (_hashtab_hashCS[hv1][i]._count == numTrees-1)
                        vec_bs.push_back(_hashtab_hashCS[hv1][i]._bs);
                    break;
                }
            }
        }
    }
}


// Implicit BPs ///////////////
//void
//HashMap::hashing_bs_nbits(
//	BitSet&  bs,
//	unsigned treeIdx,
//	unsigned numTaxa,
//	unsigned numTrees,
//	unsigned op_mode,
//	unsigned long hv1,
//	unsigned long hv2)
//{
//	if (treeIdx == 0) { // The first tree: just populate the bipaartitions.
//		BUCKET_STRUCT_T_HASHCS bk;
//		bk._hv2 = hv2;
//		BitSet* tempbs = new BitSet(numTaxa);
//		*tempbs = bs;
//		bk._bs = tempbs;
//		bk._count = 1;
//		_hashtab_hashCS[hv1].push_back(bk);
//	}
//	else { // The other trees: just compare the biaprtitions in the hashtable
//		if (!_hashtab_hashCS[hv1].empty()) { // if the hash table location has already been occupied
//			for (unsigned i=0; i<_hashtab_hashCS[hv1].size(); ++i) { // get the number of items in the linked list
//				if (_hashtab_hashCS[hv1][i]._hv2 == hv2) { // if the same bp is found
//					if (op_mode == 0) { // if type2 checking
//						if (*(_hashtab_hashCS[hv1][i]._bs) != bs) {// Check TYPE-II error!
//							cout << "*** TYPE-II Collision! ***\n";
//							exit(0);
//						}
//					}
//					_hashtab_hashCS[hv1][i]._count++;
//					break;
//				}
//				else {
//					if (*(_hashtab_hashCS[hv1][i]._bs) == bs) {
//						cout << "Data type overflow\n";
//						exit(0);
//					}
//				}
//			}
//		}
//	}
//}



//void
//HashMap::hashing_bs_64bits(
//	uint64_t bs64,
//	BitSet 	 &bs,
//	unsigned numBits,
//	unsigned treeIdx,
//	unsigned numTaxa,
//	unsigned numTrees,
//	unsigned op_mode)
//{
//	HV_STRUCT_T hv; // structure for hv1 and hv2
//	_HF.UHashFunc(hv, bs64, numBits); // Generate hv1 and hv2 for the bs64 for pgm
//
//	if (treeIdx == 0) { // The first tree: just populate the bipaartitions.
//		BUCKET_STRUCT_T_HASHCS bk;
//		bk._hv2 = hv.hv2;
//		BitSet* tempbs = new BitSet(numTaxa);
//		*tempbs = bs;
//		bk._bs = tempbs;
//		bk._count = 1;
//		_hashtab_hashCS[hv.hv1].push_back(bk);
//	}
//	else { // The other trees: just compare the biaprtitions in the hashtable
//		if (!_hashtab_hashCS[hv.hv1].empty()) { // if the hash table location has already been occupied
//			for (unsigned i=0; i<_hashtab_hashCS[hv.hv1].size(); ++i) { // get the number of items in the linked list
//				if (_hashtab_hashCS[hv.hv1][i]._hv2 == hv.hv2) { // if the same bp is found
//					if (op_mode == 2) { // if type2 checking
//						if (*(_hashtab_hashCS[hv.hv1][i]._bs) != bs) {// Check TYPE-II error!
//							cout << "*** TYPE-II Collision! ***\n";
//							exit(0);
//						}
//					}
//					_hashtab_hashCS[hv.hv1][i]._count++;
//					break; // Assumption: there is no possibility that the bipartition exists twice or more.
//				}
//				else {
//					if (*(_hashtab_hashCS[hv.hv1][i]._bs) == bs) {
//						cout << "Data type overflow\n";
//						exit(0);
//					}
//				}
//			}
//		}
//	}
//}





void
HashMap::HashMap_clear()
{
//	for (unsigned i=0; i<_hashtab_hashCS.size(); ++i) {
//		for (unsigned j=0; j<_hashtab_hashCS[i].size(); ++j) {
//			if (_hashtab_hashCS[i][j]._bs) {
//				delete _hashtab_hashCS[i][j]._bs;
//				_hashtab_hashCS[i][j]._bs = NULL;
//			}
//			_hashtab_hashCS[i].clear();
//		}
//	}
//
//	_hashtab_hashCS.clear();

//	for (unsigned i=0; i<_hashtab_hashCS.size(); ++i) {
//		for (unsigned j=0; j<_hashtab_hashCS[i].size(); ++j) {
//			if (_hashtab_hashCS[i][j]._bs) {
//				delete _hashtab_hashCS[i][j]._bs;
//				_hashtab_hashCS[i][j]._bs = NULL;
//			}
//			_hashtab_hashCS[i].clear();
//		}
//	}
//
    _hashtab_hashCS.clear();
}


void
HashMap::uhashfunc_init(
    unsigned t,	// number of trees
    unsigned n, 	// number of taxa
    float r, 					// HASHTABLE_FACTOR
    unsigned long c)						// c value in m2 > c*t*n
{
    _HF.UHashfunc_init(t, n, r, c);
}



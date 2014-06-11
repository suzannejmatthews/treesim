/**********************************************************************
 unresolvermaj
**********************************************************************

    FILE: unresolvermaj.cc


    AUTHOR:
        Copyright (C) 2006, 2007 Seung-Jin Sul
        Dept. of Computer Science
        Texas A&M University
        College Station, Texas, U.S.A.
        (contact: sulsj@cs.tamu.edu)

    HISTORY:

    06.13.2007 (v.1.0.0)

    06.14.2007 (v.1.0.1)
        - Add Consense::ConsenseTrees()

    06.21.2007 (v.1.0.3)
        - Add "delete bs_str" in HashMap::hashing_bs_without_type3_nbits_hashconsense()
          to fix memory leak
        - Use assign() in Consense::ConsenseTrees() to improve

    06.22.2007 (v.1.0.4)
        - Modify conversion routine in HashMap::hashing_bs_without_type3_nbits_hashconsense()
          to optimize

    08.06.2007 (v.1.0.6)
        - For new SC constructing algorithm

    08.16.2007
        - Consense.cpp optimization

    08.23.2007
        - New consensus tree construction routine completed.

    09.07.2007
        - Integration with 64-bit bitstring
        - vec_bucket = iter_hashcs->second; and vec_bucket[i]
          ===> (iter_hashcs->second).size(); // remove copy
        - hashmap.cc ==> hashing_bs_without_type3_64bits()

    09.11.2007 Hash function
        - hv.hv1 += _a1[i]*ii;
            hv.hv2 += _a2[i]*ii;
            ==>
            hv.hv1 += _a1[i];
            hv.hv2 += _a2[i];

    09.19.2007
        - New tree constructing routine
        - map<string, SCNode*> ==> map<string, unsigned>

    10.23.2007 Implicit bipartitions
        - Each node stores its hash code
        - If leaf, just get the a1[i] and a2[i] values
        - If internal node, get the childrens hash code and modular by m1 and m2.

    10.24.2007 Implicit BP is only applied with n-bit implementation

    10.27.2007 Optimization
        - Hashmap => vector
        - string bitstring in bucket structure => BitSet*
        - uhashfunc_init ==> two verisons (64-bit and n-bit)
        - vvec_hashcs.resize()

    11.12.2007 DEBUG
        In hashfunc.cc --> top2 = c*t*n;

    11.13.2007 data type
        int --> unsigned
        unsigned long long

    11.16.2007 DEBUG
        Prime number generator: fix fot the bigger number (from)

*/
/**********************************************************************/

//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


// My classes
#include "hashfunc.hh"
#include "hashmap.hh"
#include "SCTree.hh"
#include "SCNode.hh"
#include "bitset.hh"

// From Split-Dist
#include "label-map.hh"
#include "bitset.hh"

// Etc
#include <cassert>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdint.h>

// For newick parser
extern "C" {
#include <newick.h>
}

#include <./tclap/CmdLine.h>

#include "randomc.h"                   // define classes for random number generators
//#include "mersenne.cc"                // members of class TRandomMersenne


#define BITS                            64 // the length of bit string using in PGM
#define LEFT                            0
#define RIGHT                           1
#define ROOT                            2

// Set a random number for m1 (= Initial size of hash table)
// m1 is the closest value to (t*n)+(t*n*HASHTABLE_FACTOR)
#define HASHTABLE_FACTOR                0.2

//#define TIMECHECK
//#define NDEBUG

// if double collision is reported, increase this!
// the c value of m2 > c*t*n in the paper
static unsigned long C                   = 1000;

//static const double MAX_COLLISION_CHANCE  = 1e-6; // chance of collision
static unsigned         NUM_TREES        = 0;    // Number of trees
static unsigned         NUM_TAXA         = 0;    // Number of taxa

void
GetTaxaLabels2(
    NEWICKNODE *node,
    LabelMap &lm)
{
    if (node->Nchildren == 0) {
        string temp(node->label);
        lm.push(temp);
    } else
        for(int i=0; i<node->Nchildren; ++i)
            GetTaxaLabels2(node->child[i], lm);
}

void
dfs_resolve_one(
    SCNode* node,
    TRandomMersenne& rg,
    SCTree* sctree)
//  unsigned seeding)
{
    if (node == NULL) return;

    unsigned numChildren = node->NumChildren();

    if (numChildren == 0) {
//      cout << "leaf=" << node->name << endl;
//      b_complete = false;
    } else {
        for (unsigned i=0; i<numChildren; ++i) {
            dfs_resolve_one(node->children[i], rg, sctree);
        }

        // if numChildren > 2
        // then make the children in to a binary tree
        //
        // 1. randomly select 2 children nodes
        // 2. add two internal nodes, intNodeA and intNodeB
        // 3. dangle two choosed nodes as children of intNodeA
        // 4. dangle the other nodes as children of intNodeB

//      if (numChildren > 2 and node->name != "root")

        if (numChildren == 3) {
//          cout << "intNode=" << node->name << " " << "node->NumChildren()=" << numChildren << endl;
//          for (unsigned i=0; i<numChildren; ++i)
//              cout << node->children[i]->name << " ";
//          cout << endl;

//          int32 seed = time(0);
//          seed += seeding;
//          TRandomMersenne rg(seed);
//          rg.RandomInit(seed);

            vector<int32> vec_r;
            int32 ir;
            for (unsigned i=0; i<2; ++i) {
                ir = rg.IRandom(0, numChildren-1);
                if (find(vec_r.begin(), vec_r.end(), ir) != vec_r.end()) {
                    --i;
                    continue;
                }
                vec_r.push_back(ir);
            }

//          cout << "random=" << vec_r[0] << " " << vec_r[1] << endl;
            assert(vec_r[0] >= 0);
            assert(vec_r[1] >= 0);

            // Make a new interenal node
            SCNode* intNodeA = new SCNode();
//          SCNode* intNodeB = new SCNode();
            intNodeA->name = "intX";
//          intNodeB->name = "intNodeB";

            // Add the new internal node in nodelist of the node
            sctree->nodelist.push_back(intNodeA);
//          sctree->nodelist.push_back(intNodeB);


            // update the children
            // dangle 2 picked children as intNodeA's children
//          cout << "selected1: " << node->children[vec_r[0]]->name << endl;
//          cout << "selected2: " << node->children[vec_r[1]]->name << endl;

            // dangle two randomlly chosen children as children of the new
            // internal node
            node->children[vec_r[0]]->parent = intNodeA;
            node->children[vec_r[1]]->parent = intNodeA;
            intNodeA->children.push_back(node->children[vec_r[0]]);
            intNodeA->children.push_back(node->children[vec_r[1]]);

            // remove the moved children from the nodelist in the node
//          cout << node->NumChildren() << endl;
            node->children[vec_r[0]] = NULL;
            node->children[vec_r[1]] = NULL;
//          cout << node->NumChildren() << endl;

            assert(intNodeA->NumChildren() == 2);
//          cout << "c=" << intNodeA->children[0]->name << endl;
//          cout << "c=" << intNodeA->children[1]->name << endl;

            assert(intNodeA->children[0]->parent == intNodeA);
            assert(intNodeA->children[1]->parent == intNodeA);
//          cout << "p=" << intNodeA->children[0]->parent->name << endl;
//          cout << "p=" << intNodeA->children[1]->parent->name << endl;


            // dangle the other nodes as intNodeB's children
//          for (unsigned i=0; i<numChildren; ++i)
//          {
//              if (i != vec_r[0] && i != vec_r[1])
//              {
////                    cout << "others=" << node->children[i]->name << endl;
//                  node->children[i]->parent = intNodeB;
//                  intNodeB->children.push_back(node->children[i]);
//              }
//          }

//          assert(intNodeB->NumChildren() == numChildren-2);

//          for (unsigned i=0; i<intNodeB->NumChildren(); ++i)
//          {
//              intNodeB->children[i]->parent = intNodeB;
//          }

            // update the original node
            // clear the original children info
//          assert(node->NumChildren() == numChildren);
//          node->ClearChildren();
            assert(node->NumChildren() == numChildren-2);
            assert(intNodeA != NULL);
//          assert(intNodeB != NULL);

//          cout << intNodeA->name << endl;
//          node->children[0] = intNodeA;
            vector<SCNode*> vec_temp = node->children;
//          cout << "vec_temp.size()=" << vec_temp.size() << endl;
            vec_temp.push_back(intNodeA);
//          cout << "vec_temp.size()=" << vec_temp.size() << endl;

//          cout << node->NumChildren() << endl;
//          node->children.push_back(intNodeA);
            node->children.clear();
//          cout << node->NumChildren() << endl;

            for (unsigned j=0; j<vec_temp.size(); ++j) {
                if (vec_temp[j] != NULL) node->children.push_back(vec_temp[j]);
            }
//          cout << node->NumChildren() << endl;

//          node->children[1] = intNodeB;
//          assert(node->NumChildren() == 2);
//          assert(node->children[0] != NULL);
//          assert(node->children[1] != NULL);

            // update the new internal nodes
            intNodeA->parent = node;
//          intNodeB->parent = node;
//          cout << "intNodeA->parent=" << intNodeA->parent->name << endl;
//          cout << "intNodeB->parent=" << intNodeB->parent->name << endl;

//          b_complete = false;
        }

        else if (numChildren > 3) {
//          cout << "intNode=" << node->name << endl;
//          for (unsigned i=0; i<numChildren; ++i)
//              cout << node->children[i]->name << " ";
//          cout << endl;

//          int32 seed = time(0);
//          seed += seeding;
//          TRandomMersenne rg(seed);
//          rg.RandomInit(seed);

            vector<int32> vec_r;
            int32 ir;
            for (unsigned i=0; i<2; ++i) {
                ir = rg.IRandom(0, numChildren-1);
                if (find(vec_r.begin(), vec_r.end(), ir) != vec_r.end()) {
                    --i;
                    continue;
                }
                vec_r.push_back(ir);
            }

//          cout << "random=" << vec_r[0] << " " << vec_r[1] << endl;
            assert(vec_r[0] >= 0);
            assert(vec_r[1] >= 0);

            SCNode* intNodeA = new SCNode();
            SCNode* intNodeB = new SCNode();
            intNodeA->name = "intNodeA";
            intNodeB->name = "intNodeB";
            sctree->nodelist.push_back(intNodeA);
            sctree->nodelist.push_back(intNodeB);


            // update the children
            // dangle 2 picked children as intNodeA's children
//          cout << " " << node->children[vec_r[0]]->name << endl;
//          cout << " " << node->children[vec_r[1]]->name << endl;

            node->children[vec_r[0]]->parent = intNodeA;
            node->children[vec_r[1]]->parent = intNodeA;

            intNodeA->children.push_back(node->children[vec_r[0]]);
            intNodeA->children.push_back(node->children[vec_r[1]]);

            assert(intNodeA->NumChildren() == 2);
//          cout << "c=" << intNodeA->children[0]->name << endl;
//          cout << "c=" << intNodeA->children[1]->name << endl;

//          intNodeA->children[0]->parent = intNodeA;
//          intNodeA->children[1]->parent = intNodeA;
//          cout << "p=" << intNodeA->children[0]->parent->name << endl;
//          cout << "p=" << intNodeA->children[1]->parent->name << endl;


            // dangle the other nodes as intNodeB's children
            for (unsigned i=0; i<numChildren; ++i) {
                if (i != vec_r[0] && i != vec_r[1]) {
//                  cout << "others=" << node->children[i]->name << endl;
                    node->children[i]->parent = intNodeB;
                    intNodeB->children.push_back(node->children[i]);
                }
            }

            assert(intNodeB->NumChildren() == numChildren-2);

//          for (unsigned i=0; i<intNodeB->NumChildren(); ++i)
//          {
//              intNodeB->children[i]->parent = intNodeB;
//          }

            // update the original node
            // clear the original children info
            assert(node->NumChildren() == numChildren);
            node->ClearChildren();
            assert(node->NumChildren() == 0);
            assert(intNodeA != NULL);
            assert(intNodeB != NULL);

//          cout << intNodeA->name << endl;
            node->children[0] = intNodeA;
            node->children[1] = intNodeB;
            assert(node->NumChildren() == 2);
            assert(node->children[0] != NULL);
            assert(node->children[1] != NULL);

            // update the new internal nodes
            intNodeA->parent = node;
            intNodeB->parent = node;
//          cout << "intNodeA->parent=" << intNodeA->parent->name << endl;
//          cout << "intNodeB->parent=" << intNodeB->parent->name << endl;

//          b_complete = false;
        }
//      else
//      {
//          b_complete = true;
//      }
    }
}


void
dfs_check_multifur(
    SCNode* node,
    unsigned &intMultifur)
{
    if (node == NULL) return;

    unsigned numChildren = node->NumChildren();
    if (numChildren > 2) {
        intMultifur++;
    }

    if (numChildren == 0) {
//      cout << "leaf=" << node->name << endl;
//      b_complete = false;
    } else {
        for (unsigned i=0; i<numChildren; ++i) {
            dfs_check_multifur(node->children[i], intMultifur);
        }
//      cout << "intNode=" << node->name << endl;
//      cout << "node->NumChildren()=" << numChildren << endl;

    }
}




// IMPLICIT BP
BitSet*
dfs_hashcs_SC_nbit_wo_T2_NEWICK(
    NEWICKNODE* startNode,
    LabelMap &lm,
    HashMap &vvec_hashcs,
    unsigned treeIdx,
    unsigned long m1,
    unsigned long m2,
    vector<BitSet> & vec_bs)
{
//////  if (treeIdx == 0)
//////  { // for the first tree
    // If the node is leaf node, just set the place of the taxon name in the bit string to '1'
    // and push the bit string into stack
//      if (startNode.IsLeaf()) {
    if (startNode->Nchildren == 0) {
        // leaf node
        BitSet* bs = new BitSet(NUM_TAXA);
        string temp(startNode->label);
        unsigned idx = lm[temp];
//          unsigned idx = lm[startNode.name];
        (*bs)[idx] = true;

        // Implicit BPs /////////////////////////
        // Set the hash values for each leaf node.
        startNode->hv1 = vvec_hashcs._HF.getA1(idx);
        startNode->hv2 = vvec_hashcs._HF.getA2(idx);

        return bs;
    } else {
////            BitSet* ebs1 = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[0], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs);
////          BitSet* ebs2 = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[1], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs);


        // At this point, we find a bipartition.
        // Thus, OR the bitstrings and make a bit string for the bipartition
        BitSet* bs = new BitSet(NUM_TAXA);

        for (int i=0; i<startNode->Nchildren; ++i) {
            BitSet* ebs = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[i], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs);
            *bs |= *ebs;
            if (ebs) {
                delete ebs;
                ebs = NULL;
            }
        }

        // If binary tree then 2 edges for every inner node
////          *bs |= *ebs1;
////            if (ebs1) {
////                delete ebs1;
////                ebs1 = NULL;
////            }
////
////          *bs |= *ebs2;
////          if (ebs2) {
////                delete ebs2;
////                ebs2 = NULL;
////            }

        // Implicit BPs ////////////
        // After an internal node is found, compute the hv1 and hv2
//          startNode->hv1  = (startNode.children[LEFT]->hv1 + startNode.children[RIGHT]->hv1) % m1;
//          startNode->hv2  = (startNode.children[LEFT]->hv2 + startNode.children[RIGHT]->hv2) % m2;
        unsigned long long temp1=0;
        unsigned long long temp2=0;
        for(int i=0; i<startNode->Nchildren; ++i) {
            temp1 += startNode->child[i]->hv1;
            temp2 += startNode->child[i]->hv2;
        }
        startNode->hv1 = temp1 % m1;
        startNode->hv2 = temp2 % m2;

        // Store bit strings in hash table
//          vvec_hashcs.hashing_bs_nbits_wo_T2(*bs, treeIdx, NUM_TAXA, NUM_TREES, OP_MODE, startNode->hv1, startNode->hv2, vec_bs);
//          vvec_hashcs.hashing_bs_nbits_wo_T2(*bs, treeIdx, NUM_TAXA, NUM_TREES, startNode->hv1, startNode->hv2, vec_bs);
//          cout << *bs << endl;


        vec_bs.push_back(*bs);

        return bs;
    }
//////  }
//////  else { // for the other trees
////////        if (startNode.IsLeaf()) {
//////      if (startNode->Nchildren == 0)
//////      { // leaf node
//////          string temp(startNode->label);
//////          unsigned idx = lm[temp];
////////            unsigned idx = lm[startNode.name];
//////          BitSet* bs = NULL;
//////
//////          // Implicit BPs /////////////////////////
//////          // Set the hash values for each leaf node.
//////          startNode->hv1 = vvec_hashcs._HF.getA1(idx);
//////          startNode->hv2 = vvec_hashcs._HF.getA2(idx);
//////
//////          return bs;
//////      }
//////      else
//////      {
////////            BitSet* ebs1 = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[0], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs);
////////          BitSet* ebs2 = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[1],lm, vvec_hashcs, treeIdx, m1, m2, vec_bs);
//////        for (int i=0; i<startNode->Nchildren; ++i)
//////        {
//////              dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[0], lm, vvec_hashcs, treeIdx, m1, m2, vec_bs);
//////          }
//////          // At this point, we find a bipartition.
//////          // Thus, OR the bitstrings and make a bit string for the bipartition
////////            BitSet* bs = new BitSet(NUM_TAXA);
//////          BitSet* bs = NULL;
//////
//////          // If binary tree then 2 edges for every inner node
////////          *bs |= *ebs1;
////////            if (ebs1) {
////////                delete ebs1;
////////                ebs1 = NULL;
////////            }
//////
////////          *bs |= *ebs2;
////////          if (ebs2) {
////////                delete ebs2;
////////                ebs2 = NULL;
////////            }
//////
//////          // Implicit BPs ////////////
//////          // After an internal node is found, compute the hv1 and hv2
////////            startNode->hv1  = (startNode.children[LEFT]->hv1 + startNode.children[RIGHT]->hv1) % m1;
////////            startNode->hv2  = (startNode.children[LEFT]->hv2 + startNode.children[RIGHT]->hv2) % m2;
//////              unsigned long long temp1=0;
//////              unsigned long long temp2=0;
//////              for(int i=0; i<startNode->Nchildren; ++i)
//////              {
//////                  temp1 += startNode->child[i]->hv1;
//////                  temp2 += startNode->child[i]->hv2;
//////              }
//////              startNode->hv1 = temp1 % m1;
//////              startNode->hv2 = temp2 % m2;
//////
//////          // Store bit strings in hash table
////////            vvec_hashcs.hashing_bs_nbits_wo_T2(*bs, treeIdx, NUM_TAXA, NUM_TREES, OP_MODE, startNode->hv1, startNode->hv2, vec_bs);
//////          vvec_hashcs.hashing_bs_nbits_wo_T2(*bs, treeIdx, NUM_TAXA, NUM_TREES, startNode->hv1, startNode->hv2, vec_bs);
//////
//////          return bs;
//////      }
//////  }
}



// int to string
string
itostr(int value, int base)
{
    enum { kMaxDigits = 5 };
    std::string buf;
    buf.reserve(kMaxDigits); // Pre-allocate enough space.

    // check that the base if valid
    if (base < 2 || base > 16) return buf;
    int quotient = value;

    // Translating number to string with base:
    do {
        buf += "0123456789abcdef"[std::abs(quotient % base)];
        quotient /= base;
    } while (quotient);

    // Append the negative sign for base 10
    if (value < 0 && base == 10) buf += '-';

    std::reverse(buf.begin(), buf.end());

    return buf;
}



int main(int argc, char** argv)
{
    string outfilename;
    float ResolutionRate = 0.0;
    unsigned numOutputTrees = 0;

    // TCLAP
    try {
        // Define the command line object.
        string    helpMsg  = "HashCS-sc\n";

        helpMsg += "Input file: \n";
        helpMsg += "   The current version of HashCS only supports the Newick format.\n";

        helpMsg += "Example of Newick tree: \n";
        helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
        helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
        helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
        helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";

        helpMsg += "File option: (default = outtree)\n";
        helpMsg += "   -o <export-file-name>, specify a file name to save the result tree.\n";

        helpMsg += "Examples: \n";
        helpMsg += "  hashcs-sc foo.tre 1000\n";
        helpMsg += "  hashcs-sc foo.tre 1000 -o out.tre\n";

        TCLAP::CmdLine cmd(helpMsg, ' ', "0.9");

        TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
        cmd.add( fnameArg );

        TCLAP::UnlabeledValueArg<int>  numtreeArg( "numtree", "number of trees", true, 2, "Number of trees"  );
        cmd.add( numtreeArg );

        TCLAP::UnlabeledValueArg<float>  rateArg( "rate", "resolution rate", true, 0.0, "Resolution rate"  );
        cmd.add( rateArg );

        TCLAP::UnlabeledValueArg<unsigned>  numoutTreeArg( "numTree", "number of output trees", true, 1, "Number of output trees"  );
        cmd.add( numoutTreeArg );

        TCLAP::ValueArg<int> cArg("c", "cvalue", "c value", false, 1000, "c value");
        cmd.add( cArg );

        TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "outtree", "Output file name");
        cmd.add( outfileArg );

        cmd.parse( argc, argv );


        //  NUM_TREES = atoi(argv[2]);
        NUM_TREES = numtreeArg.getValue();
        ResolutionRate = rateArg.getValue();
        numOutputTrees = numoutTreeArg.getValue();

        if (NUM_TREES == 0) {
            string strFileLine;
            unsigned long ulLineCount;
            ulLineCount = 0;

            ifstream infile(argv[1]);

            if (infile) {
                while (getline(infile, strFileLine)) {
                    ulLineCount++;
                }
            }
            cout << "*** Number of trees in the input file: " << ulLineCount << endl;
            NUM_TREES = ulLineCount;

            infile.close();
        }

//    if (NUM_TREES < 2)
//    {
//      cerr << "Fatal error: at least two trees expected.\n";
//      exit(2);
//    }

        if (cArg.getValue()) C = cArg.getValue();

        outfilename = outfileArg.getValue();

    } catch (TCLAP::ArgException &e) { // catch any exceptions
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }

    /*********************************************************************************/
    cout << "*** Reading a tree file and parsing the tree for taxon label collection ***\n";
    /*********************************************************************************/
    NEWICKTREE *newickTree;
    int err;
    FILE *fp;
    fp = fopen(argv[1], "r");
    if(!fp) {
        cout << "ERROR: file open error:" << argv[1] << endl;
        exit(0);
    }

    newickTree = loadnewicktree2(fp, &err);
    if(!newickTree) {
        switch(err) {
        case -1:
            printf("Out of memory\n");
            break;
        case -2:
            printf("parse error\n");
            break;
        case -3:
            printf("Can't load file\n");
            break;
        default:
            printf("Error %d\n", err);
        }
    }

    /*********************************************************************************/
    cout << "\n*** Collecting the taxon labels ***\n";
    /*********************************************************************************/
    LabelMap lm;

    try {
        GetTaxaLabels2(newickTree->root, lm);
    } catch (LabelMap::AlreadyPushedEx ex) {
        cerr << "ERROR: The label '" << ex.label << "' appeard twice in " << endl;
        exit(2);
    }
    NUM_TAXA = lm.size();
    cout << "    Number of taxa = " << NUM_TAXA << endl;
//  for (unsigned i=0; i<lm2.size(); ++i)
//      cout << lm2.name(i) << " ";
//  cout << endl;
//  printnewicknode(newickTree->root, 0);
    killnewicktree(newickTree);
    fclose(fp);

    // To print out the labels
//  for (unsigned i=0; i<lm.size(); ++i)
//      cout << lm.name(i) << " ";
//   cout << endl;



    /*********************************************************************************/
    cout << "\n*** Reading tree file and collecting bipartitions ***\n";
    /*********************************************************************************/

//  treeFile.open(argv[1]);
//  stack<BitSet*> stack_bs; // For collecting bipartitions for a tree
    HashMap vvec_hashcs; // Class HashRFMap for HashConsense

    // Init hash function class
    unsigned long M1=0;
    unsigned long M2=0;

    vvec_hashcs.uhashfunc_init(NUM_TREES, NUM_TAXA, HASHTABLE_FACTOR, C);
    M1 = vvec_hashcs._HF.getM1();
    M2 = vvec_hashcs._HF.getM2();
    vvec_hashcs._hashtab_hashCS.resize(M1); // Increase the size of hash table upto m1

    vector<BitSet> vec_bs; // to collect sc bipartitions

    fp = fopen(argv[1], "r");
    if(!fp) {
        cout << "ERROR: file open error\n";
        exit(0);
    }

    for (unsigned treeIdx=0; treeIdx<NUM_TREES; ++treeIdx) {
        newickTree = loadnewicktree2(fp, &err);
        if(!newickTree) {
            switch(err) {
            case -1:
                printf("Out of memory\n");
                break;
            case -2:
                printf("parse error\n");
                break;
            case -3:
                printf("Can't load file\n");
                break;
            default:
                printf("Error %d\n", err);
            }
        } else {
            dfs_hashcs_SC_nbit_wo_T2_NEWICK(newickTree->root, lm, vvec_hashcs, treeIdx, M1, M2, vec_bs);
            killnewicktree(newickTree);
        }
    }

    unsigned total_BPs = vec_bs.size()-2;
    cout << "    vec_bs.size() = " << vec_bs.size()-2 << endl;
    cout << "    Number of trees = " << NUM_TREES << endl;
    fclose(fp);

//  cout << "    Number of BPs = " << vec_bs.size() << endl;

//  vector<BitSet> vec_bs2;

    vec_bs.erase(vec_bs.end());
    vec_bs.erase(vec_bs.end());

//  for (unsigned i=0; i<vec_bs.size(); ++i)
//  {
//      cout << vec_bs[i] << endl;
//  }

    int32 seed = time(0);                // random seed
    TRandomMersenne rg(seed);            // make instance of random number generator

    vector<int32> vec_random;

    for ( unsigned i = 0; i < total_BPs; ++i ) {
        int32 ir = rg.IRandom(0, total_BPs-1);
        if (find(vec_random.begin(), vec_random.end(), ir) != vec_random.end()) {
            --i;
            continue;
        }
        vec_random.push_back(ir);
    }

    cout << "num of random numbers=" << vec_random.size() << endl;
//  for (unsigned i=0; i<vec_random.size(); ++i)
//  {
//      cout << vec_random[i] << endl;
//  }

    cout << "Rate=" << ResolutionRate << endl;
//  cout << "Before round=" << total_BPs * ResolutionRate << endl;
//  cout << "After  round=" << round(total_BPs * ResolutionRate) << endl;


    vector<int> vec_random_maj1;

    for (unsigned i=0; i<total_BPs; ++i) {
        int ir;
        ir = rg.IRandom(numOutputTrees/2+1, numOutputTrees);
        vec_random_maj1.push_back(ir);
    }
////    for (unsigned i=0; i<vec_random_maj1.size(); ++i)
////    {
////        cout << "vec_random_maj1[i]=" << vec_random_maj1[i] << endl;
////    }

    vector<int> vec_random_maj2;
    map<int, vector<int> > map_random_maj2;

    for (unsigned i=0; i<vec_random_maj1.size(); ++i) {

        if (vec_random_maj1[i] != numOutputTrees) {
            vector<int> vec_random_maj2;
            for (unsigned j=0; j<numOutputTrees; ++j) {
                int32 ir;
                ir = rg.IRandom(1, numOutputTrees);
                if (find(vec_random_maj2.begin(), vec_random_maj2.end(), ir) != vec_random_maj2.end()) {
                    --j;
                    continue;
                }
                vec_random_maj2.push_back(ir);
            }
            map_random_maj2.insert( pair<int, vector<int> >(i, vec_random_maj2) );
        }
    }
////    cout <<  "map_random_maj2.size()=" << map_random_maj2.size() << endl;

    map<int, vector<int> >::iterator iter;
////    for (iter=map_random_maj2.begin(); iter!=map_random_maj2.end(); ++iter)
////    {
////        cout << iter->first << " " << iter->second.size() << " : ";
////        for (unsigned i=0; i<iter->second.size(); ++i)
////            cout << iter->second[i] << " ";
////        cout << endl;
////    }




    multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster;
    vector<vector<SCNode*> > vvec_distinctClusters2;
    ofstream fout;
    if (outfilename != "outtree")
        fout.open(outfilename.c_str());
    else
        fout.open("outtree");


    int32 seed2 = time(0);
    TRandomMersenne rg2(seed2);

    ////////////////////////////////////////////////////////////
    for (unsigned numOut=0; numOut<numOutputTrees; ++numOut) {
        cout << numOut << endl;

//  unsigned numOut=1;

        // select BPs for consensus tree
        vector<BitSet> vec_bs_selected;
        for (unsigned i=0; i<vec_bs.size(); ++i) {
            if (vec_random[i] < round(total_BPs * ResolutionRate)) {
                if (vec_random_maj1[i] == numOutputTrees)
                    // it's strict consensus bipartition, just insert.
                    vec_bs_selected.push_back(vec_bs[i]);
                else if (vec_random_maj1[i] < numOutputTrees) {
                    // if it's majority consensus bipartition,
                    // based on map_random_maj2, decide whether a bipartition to be inserted in the
                    // vec_bs_selected or not.
                    if (map_random_maj2[i][numOut] <= vec_random_maj1[i]) {
////                    cout << "map_random_maj2[" << i << "][" << numOut << "]=" << map_random_maj2[i][numOut] << endl;
                        vec_bs_selected.push_back(vec_bs[i]);
                    }
                }
            }
        }

        BitSet* bs = new BitSet(NUM_TAXA);
        for (unsigned i=0; i<NUM_TAXA; ++i)
            (*bs)[i]=true;
//  cout << *bs << endl;
        vec_bs_selected.push_back(*bs);
//  vec_random_maj1.push_back(numOutputTrees);

//  for (unsigned i=0; i<vec_bs_selected.size(); ++i)
//  {
//      cout << vec_bs_selected[i] << endl;
//  }

//  vector<vector<SCNode*> > vvec_distinctClusters2;
        for (unsigned i=0; i<vec_bs_selected.size(); ++i) {
            vector<SCNode*> vec_nodes2;
            for (unsigned j=0; j<NUM_TAXA; ++j) {
//          cout << (vec_bs[i])[j] << endl;
                if ((vec_bs_selected[i])[j]) {
                    SCNode* aNode = new SCNode();
                    aNode->name = lm.name(j);
                    vec_nodes2.push_back(aNode);
                }
            }
            vvec_distinctClusters2.push_back(vec_nodes2);
        }

        /////////////////////////////////
        vvec_hashcs.HashMap_clear();
//  vec_bs.clear();
        vec_bs_selected.clear();
//  vec_random.clear();
//  vec_random_maj1.clear();

        // Create multimap with size of the cluster as key.
//  multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster;

        // Insert the size of distict vector and the index
        // To sort the distinct vectors by the size of clusters
        for (unsigned i=0; i<vvec_distinctClusters2.size(); ++i)
            mmap_cluster.insert(multimap<unsigned,unsigned>::value_type(vvec_distinctClusters2[i].size(), i));

        /*********************************************************************************/
//  cout << "\n*** Construct consensus tree ***\n";
        /*********************************************************************************/

        //////////////////////////////////////////////////////////////////////////////
        // 09.19.2007
        // Construct SC tree
        //////////////////////////////////////////////////////////////////////////////
        multimap<unsigned,unsigned>::iterator itr;
        SCTree *scTree = new SCTree();
        bool addedRoot=false;
        unsigned intNodeNum = 0;

        for (itr=mmap_cluster.begin(); itr!=mmap_cluster.end(); ++itr) {
            if (!addedRoot) {
                // The first cluster has all the taxa.
                // This constructs a star tree with all the taxa.
                // 1. Dangle all the taxa as root's children by adjusting parent link.
                // 2. Push all the node* in the root's children
                // 3. Push all the nodes in the tree's nodelist.
                // 4. Push all the nodes' parent in the tree's parentlist.

                for (unsigned i=0; i<vvec_distinctClusters2[itr->second].size(); ++i) {
                    vvec_distinctClusters2[itr->second][i]->parent = scTree->root;
                    scTree->root->children.push_back(vvec_distinctClusters2[itr->second][i]);
                    scTree->nodelist.push_back(vvec_distinctClusters2[itr->second][i]);
                    assert(scTree->nodelist[0]->name == "root");
                    scTree->parentlist2.insert(map<string,int>::value_type(vvec_distinctClusters2[itr->second][i]->name, 0));
                }
                addedRoot = true;
            } else {
                // For the next biggest cluster,
                // 1. Find node list to move (= vvec_distinctClusters2[itr->second]) and
                //    Get the parent node of the to-go nodes.
                // 2. Make an internal node.
                // 3. Insert the node in the nodelist of the tree and update the parentlist accordingly
                // 4. Adjust the to-go nodes' parent link.
                // 5. Adjust the parent node's link to children (delete the moved nodes from children).


                // 1. --------------------------------------------------------------------------
                SCNode* theParent = NULL;
                theParent = scTree->nodelist[scTree->parentlist2[vvec_distinctClusters2[itr->second][0]->name]];

                assert(theParent != NULL);
                assert(theParent->name != "");

                // 2. --------------------------------------------------------------------------
                string newIntNodeName = "int" + itostr(intNodeNum, 10);
                SCNode* newIntNode = new SCNode();
                newIntNode->name = newIntNodeName;
                newIntNode->parent = theParent;

                // 3. --------------------------------------------------------------------------
                assert(newIntNodeName.size() != 0);
                scTree->nodelist.push_back(newIntNode);
                assert(scTree->nodelist[scTree->nodelist.size()-1]->name == newIntNode->name);

                scTree->parentlist2.insert(map<string, unsigned>::value_type(newIntNodeName, scTree->nodelist.size()-1));

                for (unsigned i=0; i<vvec_distinctClusters2[itr->second].size(); ++i) {
                    // 4. --------------------------------------------------------------------------
                    vvec_distinctClusters2[itr->second][i]->parent = newIntNode;

                    // We have to update parentlist in the tree.
                    assert(vvec_distinctClusters2[itr->second][i]->parent->name == scTree->nodelist[scTree->nodelist.size()-1]->name);

                    scTree->parentlist2[vvec_distinctClusters2[itr->second][i]->name] = scTree->nodelist.size()-1;
                    newIntNode->children.push_back(vvec_distinctClusters2[itr->second][i]);

                    // 5. --------------------------------------------------------------------------
                    // Delete the moved nodes from parent's children.
                    vector<SCNode*>::iterator itr2;

                    for (itr2 = theParent->children.begin(); itr2 != theParent->children.end(); ++itr2) {
                        if (vvec_distinctClusters2[itr->second][i]->name == (*itr2)->name) {
                            theParent->children.erase(itr2);
                            break;
                        }
                    }
                }
                theParent->children.push_back(newIntNode);
                intNodeNum++;
            }
//      scTree->DrawOnTerminal();
        }

//  cout << scTree->GetTreeString() << endl;
//  fout << scTree->GetTreeString() << endl;


        unsigned numMultifur=0;




//  int intNodeNum = 0;
//  scTree->root = dfs_new(newickTree->root, scTree, intNodeNum);
//  scTree->root->name = "root";

        do {
            numMultifur=0;
//          dfs_resolve_one(scTree->root, numOutputTrees);
            dfs_resolve_one(scTree->root, rg2, scTree);
            dfs_check_multifur(scTree->root, numMultifur);
        } while (numMultifur != 0);

        fout << scTree->GetTreeString();
        fout << endl;

//      delete scTree;
//      scTree=NULL;
        scTree->DeleteAllNodes();
//      dfs_sctree_delete_nodes(scTree->root);

//      int intNodeNum = 0;
//      scTree->root = dfs_new(newickTree->root, scTree, intNodeNum);
//      scTree->root->name = "root";


//  scTree->DeleteAllNodes();
//  delete scTree;
//  killnewicktree(newickTree);



















//  scTree->DeleteAllNodes();
//  scTree=NULL;
        mmap_cluster.clear();
        for (unsigned i=0; i<vvec_distinctClusters2.size(); ++i)
            vvec_distinctClusters2[i].clear();
        vvec_distinctClusters2.clear();


    }
    ////////////////////////////////////////////////////////////

    /*********************************************************************************/
    cout << "\n*** Write the resulting multifurcating tree to file ***\n";
    /*********************************************************************************/
////    ofstream fout;
//////  string fname(argv[1]);
//////  char* fnamechar;
//////  string f_prefix = "SC-" + itostr(NUM_TAXA, 10) + "-" + itostr(NUM_TREES, 10) + ".tre";
//////  cout << f_prefix << endl;
//////  strcpy(fnamechar, f_prefix.c_str());
//////  if (b_outfname)
////    if (outfilename != "outtree")
//////      fout.open(argv[7]);
////        fout.open(outfilename.c_str());
////    else
////        fout.open("outtree");
//////  fout.open(fnamechar);
//////  assert(scTree->GetTreeString().size() > 10);  // Just in case
////    fout << scTree->GetTreeString();
    fout.close();

//  cout << "SC-" + itostr(NUM_TAXA, 10) + "-" + itostr(NUM_TREES, 10) << "\t" << scTree->GetTreeString() << endl;


    /*********************************************************************************/
//  cout << "\n*** Clear memory allocated ***\n";
    /*********************************************************************************/
    mmap_cluster.clear();
    map_random_maj2.clear();
//  vvec_hashcs.HashMap_clear();
//  vec_bs.clear();
//  vec_bs_selected.clear();
//  vec_random.clear();
//  vec_random_maj1.clear();
//
    for (unsigned i=0; i<vvec_distinctClusters2.size(); ++i)
        vvec_distinctClusters2[i].clear();
//      for (unsigned j=0; i<vvec_distinctClusters2[i].size(); ++j)
//          if (vvec_distinctClusters2[i][j]) {
//              delete vvec_distinctClusters2[i][j];
//              vvec_distinctClusters2[i][j] = NULL;
//          }
    vvec_distinctClusters2.clear();

//  assert(scTree != NULL);
//  free_SCTree(scTree);
//  if (scTree) {
//      delete scTree;
//      scTree = NULL;
//  }


//---------------------------------------------------------------------------
//printf("\n\nCollect bipartitions: %lds and %ldus\n", echodelay4.tv_sec, echodelay4.tv_usec);
//printf("Distinct BP: %lds and %ldus\n", echodelay2.tv_sec, echodelay2.tv_usec);
//printf("Construct SC: %lds and %ldus\n", echodelay3.tv_sec, echodelay3.tv_usec);
////printf("Delete trees: %lds and %ldus\n", echodelay6.tv_sec, echodelay6.tv_usec);
//printf("GetTreeString(): %lds and %ldus\n", echodelay6.tv_sec, echodelay6.tv_usec);
//printf("Total CPU time2: %lds and %ldus\n", echodelay.tv_sec, echodelay.tv_usec);
////printf("Collect bipartition and compute distance: %lds and %ldus\n", echodelay5.tv_sec+echodelay7.tv_sec, echodelay5.tv_usec+echodelay7.tv_usec);
////---------------------------------------------------------------------------

    // CPU time comsumed
    struct rusage a;
    if (getrusage(RUSAGE_SELF,&a) == -1) {
        cerr << "ERROR: getrusage failed.\n";
        exit(2);
    }
    cout << "\n    Total CPU time: " << a.ru_utime.tv_sec+a.ru_stime.tv_sec << " sec and ";
    cout << a.ru_utime.tv_usec+a.ru_stime.tv_usec << " usec.\n";


    return 1;
}



/**********************************************************************
 ** unresolverstr
 ***********************************************************************
 
    FILE: unresolverstr.cc from hashcs.cc


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
        - In hashfunc.cc --> top2 = c*t*n;

    11.13.2007 data type
        - int --> unsigned
        - unsigned long long

    11.16.2007 DEBUG
        - prime number generator: fix fot the bigger number (from)

      07.11.2008 New ditribution method for TCBB'08 paper
        - Serial random number gegnerator is replaced with random shuffling

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

// From Split-Dist
#include "label-map.hh"

// Etc
#include <cassert>
#include <sys/time.h>
#include <sys/resource.h>

// For newick parser
extern "C" {
#include <newick.h>
}

#include <./tclap/CmdLine.h>
#include "randomc.h"                   // define classes for random number generators

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
static unsigned         NUM_TAXA         = 0;    // Number of taxa

void GetTaxaLabels2(NEWICKNODE *node, LabelMap &lm) { 
  //recursive procedure to get all the taxa in a particular tree. Initially, node points to root.
  if (node->Nchildren == 0) {
    string temp(node->label);
    lm.push(temp);
  } 
  else {
    for(int i=0; i<node->Nchildren; ++i)
      GetTaxaLabels2(node->child[i], lm);
  }
}

void dfs_resolve_one(SCNode* node, TRandomMersenne& rg, SCTree* sctree, vector<SCNode*> &vec_trashcan_SCNODEp) {

  if (node == NULL) return; 
  //newly added shortcut: if the child is named intNodeA or intX, we know that it is binary already.
  if (node->name == "intX" || node->name=="intNodeA") return;   

  unsigned numChildren = node->NumChildren();
  //cout << "hello! my name is: " << node->name << endl;

  if (numChildren != 0) { //if this node is not a leaf
    //cout << "i am not a leaf... recursing" << endl;
    for (unsigned i=0; i<numChildren; ++i) {
      dfs_resolve_one(node->children[i], rg, sctree, vec_trashcan_SCNODEp); //recursively call the procedure until we hit a leaf node
    }

    // if numChildren > 2
    // then make the children in to a binary tree
    //
    // 1. randomly select 2 children nodes
    // 2. add two internal nodes, intNodeA and intNodeB
    // 3. dangle two choosed nodes as children of intNodeA
    // 4. dangle the other nodes as children of intNodeB
    //cout << "out of recusion! my name is:" << node->name << endl;
    if (numChildren == 3) {

      //cout << "I have exactly three children." << endl;
      //for (unsigned int i = 0; i < node->children.size(); i++){
      //if (node->children[i] == NULL)
      ///  continue;
      //cout << node->children[i]->name << endl;
      //}

      vector<int> vec_r;
      int32 ir;
      for (unsigned i=0; i<2; ++i) {
	ir = rg.IRandom(0, numChildren-1);
	if (find(vec_r.begin(), vec_r.end(), ir) != vec_r.end()) {
	  --i;
	  continue;
	}
	vec_r.push_back(ir);
      }
           
      assert(vec_r[0] >= 0);
      assert(vec_r[1] >= 0);
      
      // Make a new internal node
      SCNode* intNodeA = new SCNode();
      vec_trashcan_SCNODEp.push_back(intNodeA);
      
      intNodeA->name = "intX";
      
      // Add the new internal node in nodelist of the node
      sctree->nodelist.push_back(intNodeA);

      // dangle two randomly chosen children as children of the new
      // internal node
      node->children[vec_r[0]]->parent = intNodeA;
      node->children[vec_r[1]]->parent = intNodeA;
      intNodeA->children.push_back(node->children[vec_r[0]]);
      intNodeA->children.push_back(node->children[vec_r[1]]);

      // remove the moved children from the nodelist in the node
      node->children[vec_r[0]] = NULL;
      node->children[vec_r[1]] = NULL;
      
      assert(intNodeA->NumChildren() == 2);
      
      assert(intNodeA->children[0]->parent == intNodeA);
      assert(intNodeA->children[1]->parent == intNodeA);
      
      // update the original node
      // clear the original children info
      assert(node->NumChildren() == numChildren-2);
      assert(intNodeA != NULL);
      
      vector<SCNode*> vec_temp = node->children;
      vec_temp.push_back(intNodeA);
      
      node->children.clear();
      
      for (unsigned j=0; j<vec_temp.size(); ++j) {
	if (vec_temp[j] != NULL) node->children.push_back(vec_temp[j]);
      }
      
      // update the new internal nodes
      intNodeA->parent = node;
      //cout << "my children are now:" << endl;
      //for (unsigned int i = 0; i < node->children.size(); i++){
      //if (node->children[i] == NULL)
      //continue;
      //cout << node->children[i]->name << endl;
      //}
      //cout << "children of new node, intX:" << endl;
      //for (unsigned int i = 0; i < intNodeA->children.size(); i++){
      //if (node->children[i] == NULL)
      //continue;
      //cout << intNodeA->children[i]->name << endl;
      //}
    } //end if numChildren == 3
    else if (numChildren > 3) {
      //cout << "I have more than three children!" << endl;
      //for (unsigned int i = 0; i < node->children.size(); i++){
      //if (node->children[i] == NULL)
      // continue;
      //cout << node->children[i]->name << endl;
      //}
      vector<int> vec_r;
      int32 ir;

      for (unsigned i=0; i<2; ++i) {
	ir = rg.IRandom(0, numChildren-1);
	if (find(vec_r.begin(), vec_r.end(), ir) != vec_r.end()) {
	  --i;
	  continue;
	}
	vec_r.push_back(ir);
      }
      
      assert(vec_r[0] >= 0);
      assert(vec_r[1] >= 0);
      //cout << "the two children that are going to node a are:" << endl;
      //cout << node->children[vec_r[0]]->name << endl;
      //cout << node->children[vec_r[1]]->name << endl;

      SCNode* intNodeA = new SCNode();
      SCNode* intNodeB = new SCNode();
      
      vec_trashcan_SCNODEp.push_back(intNodeA);
      vec_trashcan_SCNODEp.push_back(intNodeB);
      
      intNodeA->name = "intNodeA";
      intNodeB->name = "intNodeB";
      sctree->nodelist.push_back(intNodeA);
      sctree->nodelist.push_back(intNodeB);
      
      // update the children
      // dangle 2 picked children as intNodeA's children
      node->children[vec_r[0]]->parent = intNodeA;
      node->children[vec_r[1]]->parent = intNodeA;
      
      intNodeA->children.push_back(node->children[vec_r[0]]);
      intNodeA->children.push_back(node->children[vec_r[1]]);
      
      assert(intNodeA->NumChildren() == 2);
      
      // dangle the other nodes as intNodeB's children
      for (int i=0; i<node->children.size(); ++i) {
	if (node->children[i] == NULL)
	  continue;
	//cout << "i is: " << i << endl;
	//cout << "name is: " << node->children[i]->name << endl;
	if (i != vec_r[0] && i != vec_r[1]) {
	  //cout << "this node is being added as B's child.." << endl;
	  node->children[i]->parent = intNodeB;
	  intNodeB->children.push_back(node->children[i]);
	}
      }

      assert(intNodeB->NumChildren() == numChildren-2);

       // update the original node
      // clear the original children info
      assert(node->NumChildren() == numChildren);


      node->ClearChildren();
   
      assert(node->NumChildren() == 0);
      assert(intNodeA != NULL);
      assert(intNodeB != NULL);


      
      node->children[0] = intNodeA;
      node->children[1] = intNodeB;
      assert(node->NumChildren() == 2);
      assert(node->children[0] != NULL);
      assert(node->children[1] != NULL);
      
      // update the new internal nodes
      intNodeA->parent = node;
      intNodeB->parent = node;
      //cout << "my children are now:" << endl;
      //for (unsigned int i = 0; i < node->children.size(); i++){
      //if (node->children[i] == NULL){
      // cout << "(null)" << endl;
      //continue;
      //}
      //cout << node->children[i]->name << endl;
      //}
      //cout << "children of first new node, intNodeA:" << endl;
      //for (unsigned int i = 0; i < intNodeA->children.size(); i++){
      //if (intNodeA->children[i] == NULL){
      //cout << "(null)" << endl;
      //continue;
      //}
      //cout << intNodeA->children[i]->name << endl;
      //}
      //cout << "children of second new node, intNodeB:" << endl;
      //for (unsigned int i = 0; i < intNodeB->children.size(); i++){
      //if (intNodeB->children[i] == NULL){
      //cout << "(null)" << endl;
      //continue;
      //}
      //cout << intNodeB->children[i]->name << endl;
      //}
      if (intNodeB->NumChildren() > 2){ //add this check to recursively resolve anytime B contains more than 2 nodes
	dfs_resolve_one(intNodeB, rg, sctree, vec_trashcan_SCNODEp); //recursively call the procedure until we hit a leaf node
      }
    } // end numchildren greater than 3
    //cout << "I have " << node->NumChildren()  << " children. My name is: " << node->name << ". Exiting recusion!" << endl;
  } //end num children greater than 0
}


void dfs_check_multifur(SCNode* node, unsigned &intMultifur) {
  //this function determines if a particular tree is multifurcating or not
  //if at any point there is a multifurcation (children > 3), it returns with
  //initMultifur incremented by one.

  if (node == NULL) return;
  
  unsigned numChildren = node->NumChildren();
  if (numChildren > 2) {
    intMultifur++;
    return;
  }

  if (numChildren != 0) {
    for (unsigned i=0; i<numChildren; ++i) {
      dfs_check_multifur(node->children[i], intMultifur);
    }
  }
}



bool * dfs_collect_bp(NEWICKNODE* startNode, LabelMap &lm, vector<bool *> & vec_bs_resolved) {
  //collects bipartitions in a tree (seung's note: "implicit bp")
  if (startNode->Nchildren == 0) { //leaf node
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    string temp(startNode->label);
    unsigned idx = lm[temp];
    bs[idx] = 1;
    //return a bitstring with the taxon associated with that position set to 1
    return bs; 
  } else {
    // At this point, we find a bipartition.
    // Thus, OR the bitstrings and make a bit string for the bipartition
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    for (int i=0; i<startNode->Nchildren; ++i) {
      bool * ebs = dfs_collect_bp(startNode->child[i], lm, vec_bs_resolved);
      for (unsigned int j = 0; j < NUM_TAXA; j++){
	  bs[j] |= ebs[j];
      }
    }

    if (find(vec_bs_resolved.begin(), vec_bs_resolved.end(), bs) == vec_bs_resolved.end()) {
      vec_bs_resolved.push_back(bs);
    }
    return bs;
  }
}

// IMPLICIT BP
bool * dfs_hashcs_SC_nbit_wo_T2_NEWICK(NEWICKNODE* startNode, LabelMap &lm, vector<bool *> & vec_bs) {
  if (startNode->Nchildren == 0) { //leaf node
    //cout << "I am a leaf node!" << endl;
    //cout << "my name is: " << startNode->label << endl;
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    string temp(startNode->label);
    unsigned idx = lm[temp];
    //cout << "idx is: " << idx << endl;
    //cout << "bs is first: ";
    //for (unsigned int i = 0; i < NUM_TAXA; i++)
    //  cout << bs[i];
    //cout << endl;
    bs[idx] = 1;
    //cout << "i am now: "; 
    //for (unsigned int i = 0; i < NUM_TAXA; i++)
    //cout << bs[i];
    //cout << endl;
    return bs;
  } else { //internal node
        
    //cout << "I am an internal node! I have: " << startNode->Nchildren << " children."<< endl;
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    //assert(startNode->Nchildren < 3);
    
    for (int i=0; i<startNode->Nchildren; ++i) {
      //cout <<"recursing on first child.." << endl;
      bool * ebs = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[i], lm, vec_bs);
      //cout << "out of recusion! my child is: " << startNode->child[i]->label << endl;
      //cout << "ebs is: ";
      //for (unsigned int j = 0; j < NUM_TAXA; j++)
      //cout << ebs[j];
      //cout << endl;

      for (unsigned int j = 0; j < NUM_TAXA; j++)
	bs[j] |= ebs[j];
      //if (ebs) {
      //delete ebs;
      //ebs = NULL;
      //}
      //cout << "bs is now: " << endl;
      //for (unsigned int j = 0; j < NUM_TAXA; j++)
      //cout << bs[j];
      //cout << endl;
      //cout << "{";
      //for (unsigned int j = 0; j < NUM_TAXA; j++){
      //if (bs[j] == 1)
      //  cout << lm.name(j)<< ",";
      //}
      //cout << "}" << endl;
    }
    vec_bs.push_back(bs);
    
    return bs;
  }
}


string itos(int i)  // convert int to string
{
    stringstream s;
    s << i;
    return s.str();
}


int main(int argc, char** argv) {
  string outfilename;
  float ResolutionRate = 0.0;
  unsigned numOutputTrees = 0;
  
  // TCLAP
  try {
    // Define the command line object.
    string  helpMsg  = "unresolverstr\n";
    
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
    helpMsg += "  unresolverstr foo.tre 0.75 1000\n";
    helpMsg += "  unresolverstr foo.tre 0.75 1000\n";

    TCLAP::CmdLine cmd(helpMsg, ' ', "0.9");
    
    TCLAP::UnlabeledValueArg<string>  fnameArg( "name", "file name", true, "intree", "Input tree file name"  );
    cmd.add( fnameArg );
    
    TCLAP::UnlabeledValueArg<float>  rateArg( "rate", "resolution rate", true, 0.0, "Resolution rate"  );
    cmd.add( rateArg );
    
    TCLAP::UnlabeledValueArg<unsigned>  numoutTreeArg( "numTree", "number of output trees", true, 1, "Number of output trees"  );
    cmd.add( numoutTreeArg );
    
    TCLAP::ValueArg<int> cArg("c", "cvalue", "c value", false, 1000, "c value");
    cmd.add( cArg );
    
    TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "outtree", "Output file name");
    cmd.add( outfileArg );
    
    cmd.parse( argc, argv );
    
    ResolutionRate = rateArg.getValue();
    numOutputTrees = numoutTreeArg.getValue();
    
    if (cArg.getValue()) C = cArg.getValue();
    
    outfilename = outfileArg.getValue();
    
    } catch (TCLAP::ArgException &e) { // catch any exceptions
    
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
  
  /*********************************************************************************/
  cout << "*** Reading a single tree from input file and parsing the tree for taxon label collection ***\n";
  /*********************************************************************************/
   
  NEWICKTREE *newickTree;
  int err;
  FILE *fp;
  fp = fopen(argv[1], "r");
  if(!fp) {
    cout << "ERROR: file open error\n";
    exit(0);
  }

  newickTree = loadnewicktree2(fp, &err); //reads in a single tree
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
  
  /**************************************************/
  cout << "\n*** Collecting the taxon labels ***\n";
  /**************************************************/
  LabelMap lm;

  try {
    GetTaxaLabels2(newickTree->root, lm);
  } catch (LabelMap::AlreadyPushedEx ex) {
    cerr << "ERROR: The label '" << ex.label << "' appeard twice in " << endl;
    exit(2);
  }
  NUM_TAXA = lm.size();
  cout << "    Number of taxa = " << NUM_TAXA << endl;
  killnewicktree(newickTree);
  
  fclose(fp);



  /*******************************************************************/
  cout << "\n*** Reading tree file and collecting bipartitions ***\n";
  /*******************************************************************/

  vector<bool *> vec_bs; // to collect sc bipartitions
  
  fp = fopen(argv[1], "r"); //reopen file
  if(!fp) {
    cout << "ERROR: file open error\n";
    exit(0);
  }

  newickTree = loadnewicktree2(fp, &err); //re-read first tree in file
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
    dfs_hashcs_SC_nbit_wo_T2_NEWICK(newickTree->root, lm, vec_bs); //get the strict consensus bipartitions (in this case, all of them)
    killnewicktree(newickTree);
  }

  
  unsigned total_BPs = vec_bs.size()-1;
  //vec_bs.erase(vec_bs.end()); //remove star bipartition

  cout << "    vec_bs.size() = " << vec_bs.size() << endl;
  fclose(fp);

  /*for (unsigned int i = 0; i < lm.size(); i++)
    cout << lm.name(i);
  cout << endl;
  cout << "printing out the bipartitions:" << endl;
  for (unsigned int i = 0; i < vec_bs.size(); i++){
    cout << "{";
    for (unsigned int j = 0; j < NUM_TAXA; j++){
      if (vec_bs[i][j] == 1)
	cout << lm.name(j) << ",";
    }
    cout << "}" << endl;
    } */

  //return 2;
  
  // Remove two biaprtitions unnecessary (WHY DO WE DO THIS?)

  //vec_bs.erase(vec_bs.end());
  

  int32 seed = time(0);                // random seed
  TRandomMersenne rg(seed);            // make instance of random number generator


  vector<int> vec_random;

  for (unsigned i=0; i<total_BPs; ++i) { //this is the total number of bipartitions
    vec_random.push_back(i);
  }
  
  
  random_shuffle(vec_random.begin(), vec_random.end()); //shuffle the order of numbers
  
  cout << "Rate=" << ResolutionRate << endl;
  cout << "Before round=" << total_BPs * ResolutionRate << endl;
  cout << "After  round=" << round(total_BPs * ResolutionRate) << endl;
  
  unsigned numBPLimit = int(round(total_BPs * ResolutionRate));



  // Select r% (numBPLimit) bipartitoins from vec_bs

  vector<bool *> vec_bs_in;

  for (unsigned i=0; i<total_BPs; ++i) {
    if (vec_random[i] < numBPLimit)
      vec_bs_in.push_back(vec_bs[i]);
  }
  
  unsigned remaining_BPs = vec_bs.size()-1 - vec_bs_in.size();
  cout << "remaining BPs=" << remaining_BPs << endl;
  
  
  ///////////////
  //for (unsigned int i = 0; i < vec_bs.size(); i++){
  //  bool * garbage = vec_bs[i];
  //  vec_bs[i] = NULL;
  //  delete[] garbage;
  //}
  //vec_bs[i] = NULL;
  //}
  vec_bs.clear(); //we don't care about this array anymore -- now all we care about is vec_bs_in
  ///////////////
  

  cout << "# of bipartitions considered (in) = " << vec_bs_in.size() << endl;

  multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster;
  vector<vector<SCNode*> > vvec_distinctClusters2;
  ofstream fout;
  
  if (outfilename != "output.tre")
        fout.open(outfilename.c_str());
    else
      fout.open("output.tre");
  
  int32 seed2 = time(0);
  TRandomMersenne rg2(seed2);

  vector<SCNode*> vec_trashcan_SCNODEp;
  vector<string> vec_trashcan_STRING;
  vector<bool *> vec_bs_selected;
  vector<bool *> vec_bs_newly_resolved;
  
  vector< vector<int> > vvec_assigned_BID(numOutputTrees);
  map<unsigned, vector<unsigned> > map_assigned_BID;
  map<unsigned, vector<unsigned> >::iterator mItr;
  
  ////////////////// 1 ///////////////////////////////////
  //
  // Collect r% strict biaprtition only
  //
  ////////////////// 1 ///////////////////////////////////
  //to each of the trees being generated, assign the strict consensus bipartitions.
  cout << "Assigning Strict Consensus Bipartitions to each tree..." << endl;
  for (unsigned numOut=0; numOut<numOutputTrees; ++numOut) {
    //cout << numOut << endl;
    
    for (unsigned i=0; i<vec_bs_in.size(); ++i) {
      map_assigned_BID[numOut].push_back(i);
    }
  }
  cout << "Done." << endl;

  ////////////////// 2 ///////////////////////////////////
  //
  // Construct and resolve each tree to find newly
  // resolved bipartitions
  // If there is no conflict between the newly resolved
  // bipartition and the already collected BPs in vec_bs_in
  // then insert the bipartition in the BP list of the tree
  //
  ////////////////// 2 ///////////////////////////////////



  // Collect actual bipartitoins and make cluster data structure one tree by one.
  // For each cluster information, a multifurcating tree (scTree) is constructed.
  // And them the tree will be resolved one by one and the newly found (resolved)
  // bipartitions are distributed to the other trees.


  //create a star bipartition
  bool * bs = new bool[NUM_TAXA];
  for (unsigned i=0; i<NUM_TAXA; ++i) {
    bs[i]=1;
  }
    
  //beging building trees
  unsigned printMod = numOutputTrees/10;
  short percent = 0;
  cout << "Building trees..." << endl;
  for (unsigned numOut=0; numOut<numOutputTrees; ++numOut) {
    //cout << numOut << endl;
    if (numOut % printMod == 0){
      cout << "Building tree " << numOut << "(" << percent << "% done)" << endl;
      percent+=10;
    }
    //add star bipartition 
    vec_bs_selected.push_back(bs);
    
    // add strict bipartitions corresponding to tree to vec_bs_selected
    for (unsigned i=0; i<map_assigned_BID[numOut].size(); ++i) { 
      vec_bs_selected.push_back(vec_bs_in[map_assigned_BID[numOut][i]]);
    }
    
    // Clear new resolved bipartition
    vec_bs_newly_resolved.clear();
    
    for (unsigned i=0; i<vec_bs_selected.size(); ++i) {
      vector<SCNode*> vec_nodes2;
      for (unsigned j=0; j<NUM_TAXA; ++j) {
	if ((vec_bs_selected[i])[j]) {
	  SCNode* aNode = new SCNode();
	  aNode->name = lm.name(j);
	  vec_nodes2.push_back(aNode);
	  vec_trashcan_SCNODEp.push_back(aNode);
	}
      }
      vvec_distinctClusters2.push_back(vec_nodes2);
    }
    

    // Insert the cluster into a multi-map for sorting in descending order of
    // the number of '1's.
    // This is to construct tree using the bipartitions

    for (unsigned i=0; i<vvec_distinctClusters2.size(); ++i) {
      mmap_cluster.insert(multimap<unsigned,unsigned>::value_type(vvec_distinctClusters2[i].size(), i));
    }
    
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
	string newIntNodeName = "int" + itos(intNodeNum);
	vec_trashcan_STRING.push_back(newIntNodeName);
	
	SCNode* newIntNode = new SCNode();
	vec_trashcan_SCNODEp.push_back(newIntNode);
	
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
    }

    //      cout << scTree->GetTreeString() << endl;
    //      fout << scTree->GetTreeString() << endl;
    


    // The constructed trees are multifurcating trees. Thus, resolve it.

    unsigned numMultifur=0;
    
    do {
      numMultifur=0;
      dfs_resolve_one(scTree->root, rg2, scTree, vec_trashcan_SCNODEp);
      dfs_check_multifur(scTree->root, numMultifur);
    } while (numMultifur != 0);

    //=========unfortunately, this tree is not fully resolved! Let's fix this first!!
    ofstream fout_resolved;
    fout_resolved.open("resolved_tree.tre");
    string temptree = scTree->GetTreeString();
    
    fout_resolved << temptree << endl;
    fout_resolved.close();

    //return 2; //remember to comment this out!


    // Read the fully resolved tree and collect bipartitions

    fp = fopen("resolved_tree.tre", "r");
    if(!fp) {
      cout << "ERROR: file open error\n";
      exit(0);
    }

    NEWICKTREE *newickTree = loadnewicktree2(fp, &err);
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
      dfs_collect_bp(newickTree->root, lm, vec_bs_newly_resolved);
      killnewicktree(newickTree);
    }
    fclose(fp);
    
    //AGAIN, WHY THE HELL ARE WE DOING THIS?
    //vec_bs_newly_resolved.erase(vec_bs_newly_resolved.end());
    //vec_bs_newly_resolved.erase(vec_bs_newly_resolved.end());
    
    vector<int> temp_rand;
    
    for (unsigned j=0; j<numOutputTrees; ++j) {
      temp_rand.push_back(j+1);
    }
    
    random_shuffle(temp_rand.begin(), temp_rand.end());
    
    int newlyInserted=0;
    for (unsigned i=0; i<vec_bs_newly_resolved.size(); ++i) {
      if (find(vec_bs_in.begin(), vec_bs_in.end(), vec_bs_newly_resolved[i]) == vec_bs_in.end()) {
	vec_bs_in.push_back(vec_bs_newly_resolved[i]);
	
	for (unsigned k=0; k<numOutputTrees; ++k) {
	  unsigned numTreeToDitribute = 0;
	  numTreeToDitribute = int(round(numOutputTrees*ResolutionRate));
	  
	  if (temp_rand[k] < numTreeToDitribute) {
	    if (map_assigned_BID[k].size() != NUM_TAXA-3) {
	      assert(find(map_assigned_BID[k].begin(), map_assigned_BID[k].end(), vec_bs_in.size()-1) == map_assigned_BID[k].end());

	      bool bCheck = true;

	      for (unsigned m=0; m<map_assigned_BID[k].size(); ++m) {
		// Check conflict in bipartitions
		bool * temp_bs = new bool[NUM_TAXA];
		for (unsigned int z = 0; z < NUM_TAXA; z++)
		  temp_bs[z] = vec_bs_in[map_assigned_BID[k][m]][z];
		for (unsigned int z = 0; z < NUM_TAXA; z++)
		  temp_bs[z] &= vec_bs_newly_resolved[i][z];
		unsigned int countOnes = 0;
		for (unsigned int z = 0; z < NUM_TAXA; z++){
		  if (temp_bs[z] == 1)
		    countOnes++;
		}

		if (countOnes < 2) {
		  bCheck = false;
		  break;
		} else {
		  if (map_assigned_BID[k].size() != NUM_TAXA-3 and find(map_assigned_BID[k].begin(), map_assigned_BID[k].end(), vec_bs_in.size()-1) == map_assigned_BID[k].end()) {
		    newlyInserted++;
		  } else break;
		}
		delete[] temp_bs;
	      } //end for
	      
	    }
	  }
	}
	
      }
    }

    vec_bs_selected.clear();
    mmap_cluster.clear();
    vvec_distinctClusters2.clear();
    
    for (unsigned ii=0; ii<vec_trashcan_SCNODEp.size(); ++ii) {
      if (vec_trashcan_SCNODEp[ii]) {
	delete vec_trashcan_SCNODEp[ii];
      }
      vec_trashcan_SCNODEp[ii] = NULL;
    }
    vec_trashcan_SCNODEp.clear();
    
    if (vec_trashcan_STRING.size()) vec_trashcan_STRING.clear();
    
  }

  ////////////////// 3 ///////////////////////////////////
  //
  // Now we have complete set of BP BIDs for each tree
  // in map_assigned_BID[][]
  //
  ////////////////// 3 ///////////////////////////////////

  cout << "Building final trees.." << endl;
  percent = 0;
  for (unsigned numOut=0; numOut<numOutputTrees; ++numOut) {
    if (numOut % printMod == 0){
      cout << "Building tree " << numOut << "(" << percent << "% done)" << endl;
      percent += 10;
    }
    vec_bs_selected.clear();
    vec_bs_selected.push_back(bs);

    // strict bipartition
    for (unsigned i=0; i<map_assigned_BID[numOut].size(); ++i) {
      vec_bs_selected.push_back(vec_bs_in[map_assigned_BID[numOut][i]]);
    }

    for (unsigned i=0; i<vec_bs_selected.size(); ++i) {
      vector<SCNode*> vec_nodes2;
      for (unsigned j=0; j<NUM_TAXA; ++j) {
	if ((vec_bs_selected[i])[j]) {
	  SCNode* aNode = new SCNode();
	  aNode->name = lm.name(j);
	  vec_nodes2.push_back(aNode);
	  vec_trashcan_SCNODEp.push_back(aNode);
	}
      }
      vvec_distinctClusters2.push_back(vec_nodes2);
    }
    
    for (unsigned i=0; i<vvec_distinctClusters2.size(); ++i)
      mmap_cluster.insert(multimap<unsigned,unsigned>::value_type(vvec_distinctClusters2[i].size(), i));
    
    
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
	string newIntNodeName = "int" + itos(intNodeNum);
	vec_trashcan_STRING.push_back(newIntNodeName);
	
	SCNode* newIntNode = new SCNode();
	vec_trashcan_SCNODEp.push_back(newIntNode);
	
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
    
    }

    //  cout << scTree->GetTreeString() << endl;
    //      fout << scTree->GetTreeString() << endl;
    
    unsigned numMultifur=0;
    
    do {
      numMultifur=0;
      dfs_resolve_one(scTree->root, rg2, scTree, vec_trashcan_SCNODEp);
      dfs_check_multifur(scTree->root, numMultifur);
    } while (numMultifur != 0);
    

    // The completed tree !!!

    //be sure to change false back to true!!
    //string temptree = scTree->GetTreeString();
    string temptree = scTree->GetTreeString(true, 100.0);

    fout << temptree;
    fout << endl;
    
    vec_bs_selected.clear();
    mmap_cluster.clear();
    vvec_distinctClusters2.clear();
    
    for (unsigned ii=0; ii<vec_trashcan_SCNODEp.size(); ++ii) {
      if (vec_trashcan_SCNODEp[ii]) delete vec_trashcan_SCNODEp[ii];
      vec_trashcan_SCNODEp[ii] = NULL;
    }
    vec_trashcan_SCNODEp.clear();
    
    if (vec_trashcan_STRING.size()) {
      vec_trashcan_STRING.clear();
    }
  }
  fout.close();
  

  // CPU time comsumed
  struct rusage a;
  if (getrusage(RUSAGE_SELF,&a) == -1) {
    cerr << "ERROR: getrusage failed.\n";
    exit(2);
  }
  cout << "\n    Total CPU time: " << a.ru_utime.tv_sec+a.ru_stime.tv_sec << " sec and ";
  cout << a.ru_utime.tv_usec+a.ru_stime.tv_usec << " usec.\n";


 
  return 0;
}
// eof




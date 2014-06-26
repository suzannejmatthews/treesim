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
       
        Copyright (c) 2014, Suzanne Matthews
        Dept. of Electrical Engineering & Computer Science
        United States Military Academy
        West Point, NY 10996, U.S.A.
        (contact: suzanne.matthews@usma.edu)

    HISTORY:

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

// if double collision is reported, increase this!
// the c value of m2 > c*t*n in the paper
static unsigned long C                   = 1000;

//static const double MAX_COLLISION_CHANCE  = 1e-6; // chance of collision
static unsigned         NUM_TREES        = 0;    // Number of trees
static unsigned         NUM_TAXA         = 0;    // Number of taxa

void GetTaxaLabels2(NEWICKNODE *node, LabelMap &lm) {
  if (node->Nchildren == 0) {
    string temp(node->label);
    lm.push(temp);
  } else
    for(int i=0; i<node->Nchildren; ++i)
      GetTaxaLabels2(node->child[i], lm);
}

vector<int> pickChildren(unsigned int pick, unsigned int numChildren) {
  vector<int> vec_r;
  int32 ir;
  for (unsigned int i = 0; i < pick; ++i) {
    ir = rand() % numChildren;
    if (find(vec_r.begin(), vec_r.end(), ir) != vec_r.end()) {
      --i;
      continue;
    }
    vec_r.push_back(ir);
  }
  for (unsigned int i = 0; i < vec_r.size(); i++)
    assert(vec_r[i] >= 0);
  return vec_r;
}

void add_internal_node(SCNode * node, SCTree *sctree, string label, vector<int> vec_r, unsigned int numChildren){
  SCNode* intNodeA = new SCNode();
  //vec_trashcan_SCNODEp.push_back(intNodeA);
  intNodeA->name = label;
      
  // Add the new internal node in nodelist of the node
  sctree->nodelist.push_back(intNodeA);

  // dangle selected children as children of the new
  // internal node
  if (vec_r.size() != 0){
    random_shuffle(vec_r.begin(), vec_r.end());
    for (unsigned int i = 0; i < vec_r.size(); i++){
      node->children[vec_r[i]]->parent = intNodeA;
      intNodeA->children.push_back(node->children[vec_r[i]]);
      // remove the moved children from the nodelist in the node
      node->children[vec_r[i]] = NULL;
      assert(intNodeA->children[i]->parent == intNodeA);
    }        
    
    assert(intNodeA->NumChildren() == vec_r.size());
      
    // update the original node
    assert(node->NumChildren() == numChildren - vec_r.size());
  }
  else{
    unsigned int count = 0;
    for (unsigned int i=0; i<node->children.size(); ++i) {
      if (node->children[i] == NULL || node->children[i]->name == "intA")
	continue;
      node->children[i]->parent = intNodeA;
      intNodeA->children.push_back(node->children[i]);
      node->children[i] = NULL;
      assert(intNodeA->children[count]->parent == intNodeA);
      count++;
    }
  }
  assert(intNodeA != NULL);
      
  vector<SCNode*> vec_temp = node->children;
  vec_temp.push_back(intNodeA);
  node->children.clear();
      
  for (unsigned j=0; j<vec_temp.size(); ++j) {
    if (vec_temp[j] != NULL) node->children.push_back(vec_temp[j]);
  }
      
  // update the new internal nodes
  intNodeA->parent = node;
  /*cout << "children of new node " << label << ":" << endl;
  for (unsigned int i = 0; i < intNodeA->children.size(); i++){
    if (node->children[i] == NULL)
      cout << "(null)" << endl;
    cout << intNodeA->children[i]->name << endl;
    }*/
}

void dfs_resolve_one(SCNode* node, SCTree* sctree) {

  if (node == NULL) return; 
  //newly added shortcut: if the child is named intX, we know that it is binary already.
  if (node->name == "intX") return;   

  unsigned numChildren = node->NumChildren();
  //cout << "hello! my name is: " << node->name << endl; //printing

  // if numChildren > 2, make subtree into a binary tree
  if (numChildren != 0) { //if this node is not a leaf
    //cout << "i am not a leaf... recursing" << endl; //printing
    for (unsigned i=0; i<numChildren; ++i) {
      dfs_resolve_one(node->children[i], sctree); //recursively call the procedure until we hit a leaf node
    }

    //cout << "out of recusion! my name is:" << node->name << endl;
    if (numChildren == 3) {

      //printing
      /*cout << "I have exactly three children." << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
      if (node->children[i] == NULL)
        continue;
      cout << node->children[i]->name << endl;
      }*/
      //printing

      vector<int> vec_r = pickChildren(2, numChildren); //pick two random children
      add_internal_node(node, sctree, "intX", vec_r, numChildren); //make a new internal node, dangling children off it
      
      //printing      
      /*cout << "my children are now:" << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
      if (node->children[i] == NULL)
	continue;
      cout << node->children[i]->name << endl;
      }*/
      //printing

    } //end if numChildren == 3
    else if (numChildren > 3) {
      //printing
      /*cout << "I have more than three children!" << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
	if (node->children[i] == NULL)
	  continue;
	cout << node->children[i]->name << endl;
	}*/
      //printing

      unsigned int first = rand() % (numChildren-1) + 1; 
      if (first ==1) 
	first = numChildren-1; //either way, it will result in 1 child on one side, and n-1 children on the other
      vector<int> vec_r = pickChildren(first, numChildren);
      if (first != numChildren-1){ //that means the number is between 2 .. numChildren -2; so create two internal nodes
	add_internal_node(node, sctree, "intA", vec_r, numChildren); //make a new internal node, dangling children off it
	assert(node->NumChildren() == numChildren-first+1);
	vec_r.clear(); // we are going to add the remainder of the nodes to intB
	add_internal_node(node, sctree, "intB", vec_r, numChildren); //make a new internal node, dangling remainder of children off it
	
      }
      else{//create one new internal node
	add_internal_node(node, sctree, "intY", vec_r, numChildren); //make a new internal node, dangling n-1 children off it
      }

      assert(node->NumChildren() == 2);
      assert(node->children[0] != NULL);
      assert(node->children[1] != NULL);
      if (first != numChildren - 1){
	assert(node->children[1]->NumChildren() == numChildren-first);  
      }
      else{
	assert(node->children[1]->NumChildren() == first);  
      }
      // update the new internal nodes


      //printing      
      /*cout << "my children are now:" << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
	if (node->children[i] == NULL){
	  cout << "(null)" << endl;
	  continue;
	  }
	cout << node->children[i]->name << endl;
	}*/
      //printing      
      
      if (node->children[0]->NumChildren() > 2){ //add this check to recursively resolve anytime B contains more than 2 nodes
	dfs_resolve_one(node->children[0], sctree); //recursively call the procedure until we hit a leaf node
      }
      if (node->children[1]->NumChildren() > 2){
	dfs_resolve_one(node->children[1], sctree); //recursively call the procedure until we hit a leaf node
      }
    } // end numchildren greater than 3
    //cout << "I have " << node->NumChildren()  << " children. My name is: " << node->name << ". Exiting recusion!" << endl; //printing
  } //end num children greater than 0
}


// IMPLICIT BP
bool * dfs_hashcs_SC_nbit_wo_T2_NEWICK(NEWICKNODE* startNode, LabelMap &lm, unsigned treeIdx, vector<bool *> & vec_bs)
{
  if (startNode->Nchildren == 0) {
    // leaf node
    bool* bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    string temp(startNode->label);
    unsigned idx = lm[temp];
    bs[idx] = 1;

    //startNode->hv1 = vvec_hashcs._HF.getA1(idx);
    //startNode->hv2 = vvec_hashcs._HF.getA2(idx);
    return bs;
  } else {
    bool* bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    for (int i=0; i<startNode->Nchildren; ++i) {
      bool * ebs = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[i], lm, treeIdx, vec_bs);
      for (unsigned int j = 0; j < NUM_TAXA; j++)
	bs[j] |= ebs[j];

      //if (ebs) {
      //delete [] ebs;
      //ebs = NULL;
      //}
    }

    // Implicit BPs ////////////
    // After an internal node is found, compute the hv1 and hv2
    //unsigned long long temp1=0;
    //unsigned long long temp2=0;
    //for(int i=0; i<startNode->Nchildren; ++i) {
    //  temp1 += startNode->child[i]->hv1;
    //  temp2 += startNode->child[i]->hv2;
    // }
    //startNode->hv1 = temp1 % m1;
    //startNode->hv2 = temp2 % m2;

    vec_bs.push_back(bs);
    
    return bs;
  }
}



// int to string
string itostr(int value, int base){
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

void print_error(int err){
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

vector<unsigned int> genRandomNums(unsigned int howmany, unsigned int max){
  vector<unsigned int> vec_rand;
  for ( unsigned i = 0; i < howmany; ++i ) {
    unsigned int ir = rand()%max;
    if (find(vec_rand.begin(), vec_rand.end(), ir) != vec_rand.end()) {
      --i;
      continue;
    }
    vec_rand.push_back(ir);
  }
  return vec_rand;
}

string buildtree(vector<bool *> vec_bs_selected, LabelMap lm){
  multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster;
  vector<vector<SCNode*> > vvec_distinctClusters2;
  //cout << "selected bipartitions:" << endl;
  for (unsigned i=0; i<vec_bs_selected.size(); ++i) {
    vector<SCNode*> vec_nodes2;
    //cout << "{";
    for (unsigned j=0; j<NUM_TAXA; ++j) {
      if ((vec_bs_selected[i])[j]) {
	SCNode* aNode = new SCNode();
	aNode->name = lm.name(j);
	//cout << aNode->name << ",";
	vec_nodes2.push_back(aNode);
      }
    }
    //cout << "}" << endl;
    vvec_distinctClusters2.push_back(vec_nodes2);
  }




  /////////////////////////////////
  //vvec_hashcs.HashMap_clear();
  vec_bs_selected.clear();
  
  // Insert the size of distict vector and the index
  // To sort the distinct vectors by the size of clusters
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
    
  }
  //cout << "tree before resolving:" << scTree->GetTreeString() << endl;
  dfs_resolve_one(scTree->root, scTree);
  string tree = scTree-> GetTreeString(); 
  //cout << "tree after resolving:" << tree << endl;
  //clean up
  scTree->DeleteAllNodes();
  
  mmap_cluster.clear();
  for (unsigned i=0; i<vvec_distinctClusters2.size(); ++i)
    vvec_distinctClusters2[i].clear();
  vvec_distinctClusters2.clear();

  //return the tree
  return tree;
}

int main(int argc, char** argv)
{
    string outfilename;
    float ResolutionRate = 0.0;
    unsigned int numOutputTrees = 0;

    // TCLAP
    try {
        // Define the command line object.
        string    helpMsg  = "unresolvermaj\n";

        helpMsg += "Input file: \n";
        helpMsg += "   The current version of unresolvermaj only supports the Newick format.\n";

        helpMsg += "Example of Newick tree: \n";
        helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
        helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
        helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
        helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";

        helpMsg += "File option: (default = output.tre)\n";
        helpMsg += "   -o <export-file-name>, specify a file name to save the result tree.\n";

        helpMsg += "Examples: \n";
        helpMsg += "  unresolvermaj foo.tre 1 0.75 1000\n";
        helpMsg += "  unresolvermaj foo.tre 1 0.75 1000 -o out.tre\n";

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

        TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "output.tre", "Output file name");
        cmd.add( outfileArg );

        cmd.parse( argc, argv );

        NUM_TREES = numtreeArg.getValue();
        ResolutionRate = rateArg.getValue();
	if (ResolutionRate > 1 || ResolutionRate < 0){
	  cerr << "Resolution rate must be between 0 and 1!" << endl;
	  return 1;
	}
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
      print_error(err);
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
    killnewicktree(newickTree);
    fclose(fp);

    /*********************************************************************************/
    cout << "\n*** Reading tree file and collecting bipartitions ***\n";
    /*********************************************************************************/

    vector<bool *> vec_bs; // to collect sc bipartitions (why do we want the strict consensus bipartitions?)

    fp = fopen(argv[1], "r");
    if(!fp) {
        cout << "ERROR: file open error\n";
        exit(0);
    }

    //collect the bipartitions from only ONE tree.
    for (unsigned int treeIdx=0; treeIdx<1; ++treeIdx) {
        newickTree = loadnewicktree2(fp, &err);
        if(!newickTree) {
	  print_error(err);
        } else {
	  dfs_hashcs_SC_nbit_wo_T2_NEWICK(newickTree->root, lm, treeIdx, vec_bs);
	  killnewicktree(newickTree);
        }
    }
    fclose(fp);

    unsigned int total_BPs = vec_bs.size()-1;
    unsigned int resolution_rate = int(total_BPs*ResolutionRate);
    cout << "resolution_rate is:" << resolution_rate << endl;
    cout << "    vec_bs.size() = " << vec_bs.size()-1 << endl;
    cout << "    Number of Output trees = " << numOutputTrees << endl;
    cout << "    Number of bipartitions that will be in 50% or more of the trees = " << resolution_rate << endl;
        
    srand(time(NULL)); //seed random number generator

    //select r% of the bipartitions at random
    vector<unsigned int> vec_random = genRandomNums(resolution_rate, total_BPs);
    cout << "Requested Rate=" << ResolutionRate << endl;
    cout << "Attempted Rate=" << float(vec_random.size())/total_BPs << endl;
    
    /*cout << "printing out all the bipartitions:" << endl;
    for (unsigned int i = 0; i < vec_bs.size(); i++){
      //cout << "{";
      for (unsigned int j = 0; j < NUM_TAXA; j++){
	cout << vec_bs[i][j];
	//if (vec_bs[j])
	//cout << lm.name(j) << ",";
      }
      cout << endl;
      //cout << "}" << endl;
      }*/
    //generate t vectors of vectors (n x t)
    vector< vector<bool*> > tree_matrix;
    tree_matrix.resize(numOutputTrees);
    //for each of the selected bipartitions:
    for (unsigned int i = 0; i < vec_random.size(); i++){
      unsigned int perc = rand() % 50 + 50; // generate a random number between 50 .. 100
      //cout << "Adding bipartition: " << vec_random[i] << " to " << perc << "% of the trees" << endl;
      perc = int((float(perc)/100) * numOutputTrees);
      //cout << "perc will actually be: " << perc << endl;
      vector<unsigned int> selected_trees = genRandomNums(perc, numOutputTrees);
      for (unsigned int j = 0; j < selected_trees.size(); j++){
	unsigned int tree_id = selected_trees[j];
	tree_matrix[tree_id].push_back(vec_bs[i]);
	//cout << "added bipartition " << vec_random[i] << "to tree " << tree_id << endl;
      }
    }
    //add the star bipartition to each node
    for (unsigned int i = 0; i < numOutputTrees; i++)
      tree_matrix[i].push_back(vec_bs[total_BPs]);
    

    //commence building
    ofstream fout;
    if (outfilename != "output.tre")
      fout.open(outfilename.c_str());
    else
      fout.open("output.tre");
    
    cout << "building trees!" << endl;
    for (unsigned numOut=0; numOut<numOutputTrees; ++numOut) {
      if (numOut % 10 == 0)
	cout << numOut << endl;
      string tree = buildtree(tree_matrix[numOut], lm);
      fout << tree << endl;
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

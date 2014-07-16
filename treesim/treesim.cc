/**********************************************************************
treesim
**********************************************************************

    FILE: treesim.cc
   
    DESCRIPTION: Creates "accurate" simluated datasets based on a given strict
    and majority consensus rate.

    Based on the software unresolverstr and unresolvermaj by SeungJin Sul
    and TreeZip by Suzanne J. Matthews

    AUTHOR:
        Copyright (C) 2014 Suzanne J. Matthews
        Dept. of Electrical Engineering & Computer Science
        United States Military Academy
        West Point, NY, 10996, USA
        (contact: suzanne.matthews@usma.edu)
       
    HISTORY:
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


#include "SCTree.h"
#include "SCNode.h"

// From Split-Dist
#include "label-map.hh"

#include <cassert>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdint.h>

// For newick parser
extern "C" {
#include <newick.h>
}

#include <iostream>
#include <./tclap/CmdLine.h>

using namespace std;

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
  unsigned int ir;
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

void add_internal_node(SCNode * node, SCTree *sctree, string label, vector<int> vec_r, unsigned int numChildren, vector<SCNode *> &vec_trashcan_SCNODEp){
  SCNode* intNodeA = new SCNode();
  vec_trashcan_SCNODEp.push_back(intNodeA);
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

}

void resolve_tree(SCNode* node, SCTree* sctree, vector<SCNode*> & trashcan) {

  if (node == NULL) return; 
  //newly added shortcut: if the child is named intX, we know that it is binary already.
  if (node->name == "intX") return;   

  unsigned numChildren = node->NumChildren();
  //cout << "hello! my name is: " << node->name << endl; //printing

  // if numChildren > 2, make subtree into a binary tree
  if (numChildren != 0) { //if this node is not a leaf
    //cout << "i am not a leaf... recursing" << endl; //printing
    for (unsigned i=0; i<numChildren; ++i) {
      resolve_tree(node->children[i], sctree, trashcan); //recursively call the procedure until we hit a leaf node
    }

    //cout << "out of recusion! my name is:" << node->name << endl;
    if (numChildren == 3) {

      //printing

      //printing

      vector<int> vec_r = pickChildren(2, numChildren); //pick two random children
      add_internal_node(node, sctree, "intX", vec_r, numChildren, trashcan); //make a new internal node, dangling children off it
      
      //printing      
      //printing

    } //end if numChildren == 3
    else if (numChildren > 3) {
      //printing

      //printing

      unsigned int first = rand() % (numChildren-1) + 1; 
      if (first ==1) 
	first = numChildren-1; //either way, it will result in 1 child on one side, and n-1 children on the other
      vector<int> vec_r = pickChildren(first, numChildren);
      if (first != numChildren-1){ //that means the number is between 2 .. numChildren -2; so create two internal nodes
	add_internal_node(node, sctree, "intA", vec_r, numChildren, trashcan); //make a new internal node, dangling children off it
	assert(node->NumChildren() == numChildren-first+1);
	vec_r.clear(); // we are going to add the remainder of the nodes to intB
	add_internal_node(node, sctree, "intB", vec_r, numChildren, trashcan); //make a new internal node, dangling remainder of children off it
	
      }
      else{//create one new internal node
	add_internal_node(node, sctree, "intY", vec_r, numChildren, trashcan); //make a new internal node, dangling n-1 children off it
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

      //printing      
      
      if (node->children[0]->NumChildren() > 2){ //add this check to recursively resolve anytime B contains more than 2 nodes
	resolve_tree(node->children[0], sctree, trashcan); //recursively call the procedure until we hit a leaf node
      }
      if (node->children[1]->NumChildren() > 2){
	resolve_tree(node->children[1], sctree, trashcan); //recursively call the procedure until we hit a leaf node
      }
    } // end numchildren greater than 3
    //cout << "I have " << node->NumChildren()  << " children. My name is: " << node->name << ". Exiting recusion!" << endl; //printing
  } //end num children greater than 0
}


bool * collect_biparts(NEWICKNODE* startNode, LabelMap &lm, unsigned treeIdx, vector<bool *> & vec_bs, unsigned int NUM_TAXA)
{
  if (startNode->Nchildren == 0) {
    // leaf node
    bool* bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    string temp(startNode->label);
    unsigned idx = lm[temp];
    bs[idx] = 1;
    return bs;
  } 
  else { //internal node
    bool* bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    for (int i=0; i<startNode->Nchildren; ++i) {
      bool * ebs = collect_biparts(startNode->child[i], lm, treeIdx, vec_bs, NUM_TAXA);
      for (unsigned int j = 0; j < NUM_TAXA; j++)
	bs[j] |= ebs[j];

      //if (ebs) {
      //delete [] ebs;
      //ebs = NULL;
      //}
    }

    vec_bs.push_back(bs);
    
    return bs;
  }
}

//itostr: from unresolvermaj
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

string buildtree(vector<bool *> vec_bs_selected, LabelMap lm, unsigned int NUM_TAXA){
  multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster;
  vector<vector<SCNode*> > vvec_distinctClusters2;
  vector<SCNode*> vec_trashcan_SCNODEp;
  vector<string> vec_trashcan_STRING;

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
	vec_trashcan_SCNODEp.push_back(aNode);
      }
    }
    //cout << "}" << endl;
    vvec_distinctClusters2.push_back(vec_nodes2);
  }
  //cout << endl;
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
  //cout << "tree before resolving:" << scTree->GetTreeString() << endl;
  resolve_tree(scTree->root, scTree, vec_trashcan_SCNODEp);
  string tree = scTree-> GetTreeString(true); 
  //cout << "tree after resolving:" << tree << endl;

  //clean up    
  for (unsigned i=0; i<vec_trashcan_SCNODEp.size(); ++i) { //empty node trashcan
    if (vec_trashcan_SCNODEp[i]){ 
      delete vec_trashcan_SCNODEp[i];
    }
  }
  if (scTree->root){
    SCNode *tmp = scTree -> root;
    delete tmp;
    scTree->root = NULL;
  }
  scTree->nodelist.clear();
  scTree->parentlist2.clear();
  delete scTree;
  //return the tree
  return tree;
}

int main(int argc, char** argv)
{
  string outfilename, startingFile;
  float majRate = 0.0, strictRate = 0.0;
  unsigned int NUM_TAXA = 0, NUM_TREES = 0, unique_trees = 0;

    // TCLAP
    try {
        // Define the command line object.
        string    helpMsg  = "treesim\n";

        helpMsg += "Input file: \n";
        helpMsg += "   The current version of treesim only supports the Newick format.\n";
        helpMsg += "Sample Newick file: \n";
        helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
        helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
        helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";
        helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";

        helpMsg += "File option: (default = output.tre)\n";
        helpMsg += "   -o <export-file-name>, specify a file name to save the result tree.\n";

        helpMsg += "Example: \n";
        helpMsg += "  treesim 500 1000 0.5 .75 -o out.tre\n";
	helpMsg += "  generates a random 500-taxa, 1000 tree dataset with a strict\n";
	helpMsg += "  consensus rate of at least 50% and a majority consensus rate\n";
	helpMsg += "  of 75%. The result is stored in out.tre\n";

        TCLAP::CmdLine cmd(helpMsg, ' ', "0.1");

        TCLAP::UnlabeledValueArg<int>  numtaxaArg( "numtaxa", "number of taxa", true, 5, "Number of taxa"  );
        cmd.add( numtaxaArg );

        TCLAP::UnlabeledValueArg<unsigned>  numtreeArg( "numtree", "number of output trees", true, 1, "Number of output trees"  );
        cmd.add( numtreeArg );

        TCLAP::UnlabeledValueArg<float>  sRateArg( "sRate", "strict consensus rate", true, 0.0, "strict rate"  );
        cmd.add( sRateArg );

        TCLAP::UnlabeledValueArg<float>  mRateArg( "mRate", "majority consensus rate", true, 0.0, "majority rate"  );
        cmd.add( mRateArg );

        TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "output.tre", "Output file name");
        cmd.add( outfileArg );

	TCLAP::ValueArg<string> sfileArg("t", "startingtree", "starting tree file name", false, "starting.tre", "input starting tree file name");
        cmd.add( sfileArg );

	TCLAP::ValueArg<int> uArg("u", "uniquetree", "number of unique trees", false, 0, "number of unique trees in file");
        cmd.add( uArg );


        cmd.parse( argc, argv );
	NUM_TAXA = numtaxaArg.getValue();
        NUM_TREES = numtreeArg.getValue();
        strictRate = sRateArg.getValue();
	majRate = mRateArg.getValue();
	unique_trees = uArg.getValue();

	if (strictRate > 1 || strictRate < 0){
	  cerr << "ERROR: strict consensus rate must be between 0 and 1!" << endl;
	  return 1;
	}
	if (majRate > 1 || majRate < 0){
	  cerr << "ERROR: strict consensus rate must be between 0 and 1!" << endl;
	  return 1;
	}
	if (strictRate > majRate){
	  cerr << "ERROR: strict consensus rate must be lesser than or equal to majority consensus rate!" << endl;
	  return 1;
	}
        outfilename = outfileArg.getValue();
	startingFile = sfileArg.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }


    //generate a random tree with the specified number of taxa
   LabelMap lm;
   ofstream fout;
   NEWICKTREE *newickTree;
   int err;
   FILE *fp;
    if (startingFile == "starting.tre"){
      cout << "We get here!" << endl;
      fprintf(stderr, "Generating a random tree with %u taxa...\n", NUM_TAXA);
      vector<bool *> random_tree_bs;
      bool * star = new bool[NUM_TAXA];
      for (unsigned int i = 0; i < NUM_TAXA; i++){
	star[i] = 1;
	string taxa;
	stringstream ss;
	ss << i+1;
	ss >> taxa;
	lm.push(taxa);
      }
      random_tree_bs.push_back(star);
      string random_tree = buildtree(random_tree_bs, lm, NUM_TAXA); 
      fout.open("starting.tre");
      fout << random_tree << endl;
      fout.close();
      fprintf(stderr, "Done. starting tree outputted to starting.tre.\n");
      random_tree_bs.clear();
      delete [] star;
    }
    else{
      //collect label map from file
      cout << "getting labels from file: " << startingFile << endl;
      fp = fopen(startingFile.c_str(), "r");
      if (!fp){
        cout << "ERROR: file open error:" << argv[1] << endl;
        exit(0);
      }
      newickTree = loadnewicktree2(fp, &err);
      if(!newickTree) {
	print_error(err);
      }
      else{
	GetTaxaLabels2(newickTree->root, lm);
	killnewicktree(newickTree);
	fclose(fp);
      }
      assert(NUM_TAXA == lm.size());
    }
    fprintf(stderr, "Building Collection...\n");
 

    vector<bool *> vec_bs;
    fp = fopen(startingFile.c_str(), "r");
    if(!fp) {
        cout << "ERROR: file open error:" << argv[1] << endl;
        exit(0);
    }

    newickTree = loadnewicktree2(fp, &err);
    if(!newickTree) {
      print_error(err);
    }
    else{
      collect_biparts(newickTree->root, lm, 0, vec_bs, NUM_TAXA);
      killnewicktree(newickTree);
    }

    fclose(fp);

    unsigned int total_BPs = vec_bs.size()-1;
    unsigned int majority_resolution_rate = int(total_BPs*majRate);
    unsigned int strict_resolution_rate = int(total_BPs*strictRate);
    unsigned int difference = majority_resolution_rate - strict_resolution_rate;
    cout << "    vec_bs.size() = " << vec_bs.size()-1 << endl;
    cout << "    Number of Output trees = " << NUM_TREES << endl;
    cout << "    Number of bipartitions that will be in all the trees = " << strict_resolution_rate << endl;
    cout << "    Number of bipartitions that will be additionally in 50% or more of trees  = " << difference << endl;
    
    srand(time(NULL)); //seed random number generator

    //select r% of the bipartitions at random
    vector<unsigned int> vec_random = genRandomNums(majority_resolution_rate, total_BPs);
    //cout << "printing out random numbers:" << endl;
    //for (unsigned int i = 0; i < vec_random.size(); i++){
    //  cout << vec_random[i] << " ";
    //}
    //cout << endl;
    cout << "Requested Majority Rate=" << majRate << endl;
    cout << "Attempted Majority Rate=" << float(vec_random.size())/total_BPs << endl;    
    cout << "Requested Strict Rate=" << strictRate << endl;
    cout << "Attempted Strict Rate=" << float(strict_resolution_rate)/total_BPs << endl;    
    //generate t vectors of vectors (n x t)
    vector< vector<bool*> > tree_matrix;
    vector<unsigned int> duplicates;
    unsigned int NUM_TO_BUILD = 0;
    if (unique_trees == 0){ //if this parameter is not specfied
      NUM_TO_BUILD = NUM_TREES;
    }
    else{
      assert(unique_trees < NUM_TREES);
      assert(unique_trees > 0);
      NUM_TO_BUILD = unique_trees;
      duplicates = genRandomNums(NUM_TREES-unique_trees, NUM_TREES);
      sort(duplicates.begin(), duplicates.end()); //sort the duplicates vector
    }
    //tree_matrix.resize(NUM_TREES);
    tree_matrix.resize(NUM_TO_BUILD);
    for (unsigned int i = 0; i < strict_resolution_rate; i++){
      //cout << "adding bipartition " << vec_random[i] << "to all the trees" << endl;
      unsigned int bipart = vec_random[i];
      for (unsigned int j = 0; j < NUM_TO_BUILD; j++){
	tree_matrix[j].push_back(vec_bs[bipart]);
      }
    }

    //for each of the remainder of the selected bipartitions:
    for (unsigned int i = strict_resolution_rate; i < vec_random.size(); i++){
      unsigned int bipart = vec_random[i];
      unsigned int perc = rand() % 51 + 50; // generate a random number between 50 .. 100
      //cout << "Adding bipartition: " << vec_random[i] << " to " << perc << "% of the trees" << endl;
      perc = int((float(perc)/100) * NUM_TO_BUILD);
      //cout << "perc will actually be: " << perc << endl;
      vector<unsigned int> selected_trees = genRandomNums(perc, NUM_TO_BUILD);
      for (unsigned int j = 0; j < selected_trees.size(); j++){
	unsigned int tree_id = selected_trees[j];
	tree_matrix[tree_id].push_back(vec_bs[bipart]);
	//cout << "added bipartition " << vec_random[i] << "to tree " << tree_id << endl;
      }
    }

    //add the star bipartition to all the trees
    for (unsigned int i = 0; i < NUM_TO_BUILD; i++)
      tree_matrix[i].push_back(vec_bs[total_BPs]);

    //commence building

    if (outfilename != "output.tre")
      fout.open(outfilename.c_str());
    else
      fout.open("output.tre");
    
    cout << "building trees!" << endl;
    unsigned int dupCount = 0;
    for (unsigned numOut=0; numOut<NUM_TO_BUILD; ++numOut) {
      if (numOut % 1000 == 0)
	cout << numOut << endl;
      string tree = buildtree(tree_matrix[numOut], lm, NUM_TAXA);
      fout << tree << endl;
      if (unique_trees != 0 && numOut == duplicates[dupCount]){
	fout << tree << endl;
	dupCount++;
      }
      //tree_matrix[numOut].clear(); //remove the bipartitions from this
    }
    fout.close();

    //clean up
    for (unsigned int i = 0; i < vec_bs.size(); i++){
      if (vec_bs[i]!=NULL){
	delete [] vec_bs[i];
	vec_bs[i] = NULL;
      }
    }
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

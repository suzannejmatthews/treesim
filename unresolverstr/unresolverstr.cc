/**********************************************************************
 ** unresolverstr
 ***********************************************************************
 
    FILE: unresolverstr.cc from hashcs.cc


    AUTHOR:
        Copyright (C) 2014 Suzanne J. Matthews
        Copyright (C) 2006, 2007 Seung-Jin Sul
        Dept. of Computer Science
        Texas A&M University
        College Station, Texas, U.S.A.
        (contact: sulsj@cs.tamu.edu)

    HISTORY:
    (for earlier history, please see old version of software)
    11.16.2007 DEBUG
        - prime number generator: fix fot the bigger number (from)

    07.11.2008 New ditribution method for TCBB'08 paper
        - Serial random number gegnerator is replaced with random shuffling

    06.14.2014 (SJM) I discovered that this software does not work anymore. 
               Proceeded to rewrite the entire damn thing.
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

static unsigned         NUM_TAXA         = 0;    // Number of taxa

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

LabelMap getLabels(char * treefile){
  NEWICKTREE *newickTree;
  int err;
  FILE *fp;
  fp = fopen(treefile, "r");
  if(!fp) {
    cout << "ERROR: file open error\n";
    exit(0);
  }
  
  newickTree = loadnewicktree2(fp, &err); //reads in a single tree
  if(!newickTree) {
    print_error(err);
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

  killnewicktree(newickTree);
  fclose(fp);

  return lm;
}

void dfs_resolve_one(SCNode* node, TRandomMersenne& rg, SCTree* sctree, vector<SCNode*> &vec_trashcan_SCNODEp) {

  if (node == NULL) return; 
  //newly added shortcut: if the child is named intX, we know that it is binary already.
  if (node->name == "intX") return;   

  unsigned numChildren = node->NumChildren();
  cout << "hello! my name is: " << node->name << endl; //printing

  if (numChildren != 0) { //if this node is not a leaf
    cout << "i am not a leaf... recursing" << endl; //printing
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

      //printing
      
      cout << "I have exactly three children." << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
      if (node->children[i] == NULL)
        continue;
      cout << node->children[i]->name << endl;
      }
      //printing

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
      //printing
      
      cout << "my children are now:" << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
      if (node->children[i] == NULL)
	continue;
      cout << node->children[i]->name << endl;
      }
      cout << "children of new node, intX:" << endl;
      for (unsigned int i = 0; i < intNodeA->children.size(); i++){
	if (node->children[i] == NULL)
	  continue;
	cout << intNodeA->children[i]->name << endl;
      }
      
      //printing

    } //end if numChildren == 3
    else if (numChildren > 3) {
      //printing
      
      cout << "I have more than three children!" << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
	if (node->children[i] == NULL)
	  continue;
	cout << node->children[i]->name << endl;
	}
      //printing
      unsigned int first = rg.IRandom(2, numChildren-2); //each node (A and B) should get some number of children (at least two per)
      vector<int> vec_r;
      int32 ir;


      for (unsigned i=0; i<first; ++i) {
	ir = rg.IRandom(0, numChildren-1); //each node (A and B) should get some number of children (at least two per)
	if (find(vec_r.begin(), vec_r.end(), ir) != vec_r.end()) {
	  --i;
	  continue; //ensure all random numbers chosen are indeed unique
	}
	vec_r.push_back(ir);
      }

      for (unsigned int x = 0; x < vec_r.size(); x++) //make sure a number was chosen for each      
	assert(vec_r[x] >= 0);

      //printing
      
      cout << "the children that are going to node a are:" << endl;
      for (unsigned int x = 0; x < vec_r.size(); x++)
	cout << node->children[vec_r[x]]->name<< endl;
            
      //printing

      SCNode* intNodeA = new SCNode();
      SCNode* intNodeB = new SCNode();
      
      vec_trashcan_SCNODEp.push_back(intNodeA);
      vec_trashcan_SCNODEp.push_back(intNodeB);
      
      intNodeA->name = "intNodeA";
      intNodeB->name = "intNodeB";
      sctree->nodelist.push_back(intNodeA);
      sctree->nodelist.push_back(intNodeB);
      
      // update the children
      // dangle first picked children as intNodeA's children
      for (unsigned int x = 0; x < vec_r.size(); x++){
	node->children[vec_r[x]]->parent = intNodeA;
	intNodeA->children.push_back(node->children[vec_r[x]]);
	node->children[vec_r[x]] = NULL;
      }
      
      assert(intNodeA->NumChildren() == first);
      
      // dangle the other nodes as intNodeB's children
      for (unsigned int i=0; i<node->children.size(); ++i) {
	if (node->children[i] == NULL)
	  continue;
	//debug
	//cout << "i is: " << i << endl;
	//cout << "name is: " << node->children[i]->name << endl;
	//cout << "this node is being added as B's child.." << endl;
	

	node->children[i]->parent = intNodeB;
	intNodeB->children.push_back(node->children[i]);
      }

      assert(intNodeB->NumChildren() == numChildren-first);
      
      // update the original node
      // clear the original children info
      //assert(node->NumChildren() == numChildren);
      
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

      //printing
      
      cout << "my children are now:" << endl;
      for (unsigned int i = 0; i < node->children.size(); i++){
	if (node->children[i] == NULL){
	  cout << "(null)" << endl;
	  continue;
	  }
	cout << node->children[i]->name << endl;
      }
      
      cout << "children of first new node, intNodeA:" << endl;
      for (unsigned int i = 0; i < intNodeA->children.size(); i++){
	if (intNodeA->children[i] == NULL){
	  cout << "(null)" << endl;
	  continue;
	}
	cout << intNodeA->children[i]->name << endl;
      }
      cout << "children of second new node, intNodeB:" << endl;
      for (unsigned int i = 0; i < intNodeB->children.size(); i++){
	if (intNodeB->children[i] == NULL){
	  cout << "(null)" << endl;
	  continue;
	}
	cout << intNodeB->children[i]->name << endl;
	}
      //printing

      if (intNodeB->NumChildren() > 2){ //add this check to recursively resolve anytime B contains more than 2 nodes
	dfs_resolve_one(intNodeB, rg, sctree, vec_trashcan_SCNODEp); //recursively call the procedure until we hit a leaf node
      }
      if (intNodeA->NumChildren() > 2){
	dfs_resolve_one(intNodeA, rg, sctree, vec_trashcan_SCNODEp); //recursively call the procedure until we hit a leaf node
      }
    } // end numchildren greater than 3
    cout << "I have " << node->NumChildren()  << " children. My name is: " << node->name << ". Exiting recusion!" << endl; //printing
  } //end num children greater than 0
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

    bs[idx] = 1;
    return bs;
  } else { //internal node
        
    //cout << "I am an internal node! I have: " << startNode->Nchildren << " children."<< endl;
    bool * bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    
    for (int i=0; i<startNode->Nchildren; ++i) {
      bool * ebs = dfs_hashcs_SC_nbit_wo_T2_NEWICK(startNode->child[i], lm, vec_bs);
      //cout << "out of recusion! my child is: " << startNode->child[i]->label << endl;

      for (unsigned int j = 0; j < NUM_TAXA; j++)
	bs[j] |= ebs[j];
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

string buildTree(vector<bool *> vec_bs_in, LabelMap lm, TRandomMersenne &rg2){ //if this is not passed by reference, it gets all messed up...
  // Collect actual bipartitions and make cluster data structure one tree by one.
  // For each cluster information, a multifurcating tree (scTree) is constructed.
  // And them the tree will be resolved one by one and the newly found (resolved)
  // bipartitions are distributed to the other trees.  

  multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster;
  vector<vector<SCNode*> > vvec_distinctClusters2;  
  vector<SCNode*> vec_trashcan_SCNODEp;
  vector<string> vec_trashcan_STRING;
  
  for (unsigned i=0; i<vec_bs_in.size(); ++i) {
    vector<SCNode*> vec_nodes2;
    for (unsigned j=0; j<NUM_TAXA; ++j) {
      if ((vec_bs_in[i])[j]) {
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
  // modified by sjm: 06.09.2014
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
  } //end creating multifurcating tree
  
  // The constructed tre is multifurcating. Thus, resolve it.
  
  dfs_resolve_one(scTree->root, rg2, scTree, vec_trashcan_SCNODEp);
  string myTree = scTree->GetTreeString();

  //clean up
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

  //return tree
  return myTree;
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
    helpMsg += "   The current version of unresolverstr only supports the Newick format.\n";
    
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
        
    TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "output.tre", "Output file name");
    cmd.add( outfileArg );
    
    cmd.parse( argc, argv );
    
    ResolutionRate = rateArg.getValue();
    numOutputTrees = numoutTreeArg.getValue();    
    outfilename = outfileArg.getValue();
    
    } catch (TCLAP::ArgException &e) { // catch any exceptions
    
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  if (ResolutionRate > 1){
    cerr << "Error! Resolution Rate must be betwen 0 and 1!" << endl;
    return 2;
  }
  
  /*********************************************************************************/
  cout << "*** Reading a single tree from input file and parsing the tree for taxon label collection ***\n";
  /*********************************************************************************/
  LabelMap lm = getLabels(argv[1]);
  NUM_TAXA = lm.size();
  cout << "    Number of taxa = " << NUM_TAXA << endl;
  /*******************************************************************/
  cout << "\n*** Reading tree file and collecting bipartitions ***\n";
  /*******************************************************************/

  vector<bool *> vec_bs; // to collect sc bipartitions
  NEWICKTREE *newickTree;
  int err;
  FILE * fp = fopen(argv[1], "r"); //reopen file
  if(!fp) {
    cout << "ERROR: file open error\n";
    exit(0);
  }

  newickTree = loadnewicktree2(fp, &err); //re-read first tree in file
  if(!newickTree) {
    print_error(err);
  } else {
    dfs_hashcs_SC_nbit_wo_T2_NEWICK(newickTree->root, lm, vec_bs); //get the strict consensus bipartitions (in this case, all of them)
    killnewicktree(newickTree);
  }
  fclose(fp);

  cout << "    vec_bs.size() = " << vec_bs.size() << endl;  
 

  /*cout << "printing out the bipartitions:" << endl;
  for (unsigned int i = 0; i < vec_bs.size(); i++){
    cout << "{";
    for (unsigned int j = 0; j < NUM_TAXA; j++){
      if (vec_bs[i][j] == 1)
	cout << lm.name(j) << ",";
    }
    cout << "}" << endl;
    }*/
  
  //generate the random set of bipartitions that will chosen for strict consensus set   
  vector<unsigned int> vec_random;
  unsigned total_BPs = vec_bs.size()-1; //we don't want to include the star bipartition in the random choice
  for (unsigned int i=0; i<total_BPs; ++i) { //this is the total number of bipartitions
    vec_random.push_back(i);
  }
  
  srand(time(NULL));
  random_shuffle(vec_random.begin(), vec_random.end()); //shuffle the order of numbers

  cout << "Rate=" << ResolutionRate << endl;
  cout << "Before round=" << total_BPs * ResolutionRate << endl;
  cout << "After  round=" << round(total_BPs * ResolutionRate) << endl;
  
  unsigned numBPLimit = int(round(total_BPs * ResolutionRate));
  cout << "numBPLimit is:" << numBPLimit << endl;
  cout << "Total BPs are:" << total_BPs << endl;

  ////////////////// 1 ///////////////////////////////////
  //
  // Collect r% strict biaprtition only
  //
  ////////////////// 1 ///////////////////////////////////
  // Select r% (numBPLimit) bipartitions from vec_bs
  vector<bool *> vec_bs_in;

  for (unsigned int i=0; i<numBPLimit; ++i) {
    int pos = vec_random[i];
    vec_bs_in.push_back(vec_bs[pos]);
  }
  
  unsigned remaining_BPs = vec_bs.size()-1 - vec_bs_in.size();
  cout << "remaining BPs=" << remaining_BPs << endl;
  
  cout << "the bipartitions that will be used:" << endl;
  for (unsigned int i = 0; i < vec_bs_in.size(); i++){
    cout << "{";
    for (unsigned int j = 0; j < NUM_TAXA; j++){
      if (vec_bs_in[i][j] == 1)
	cout << lm.name(j) << ",";
    }
    cout << "}" << endl;
  }
  
  vec_bs_in.push_back(vec_bs[total_BPs]); //add the star bipartition too
    
  ///////////////
  //deallocate the pointers we won't use (the following code block doesn't work: not sure why)
  /*
  for (unsigned int i = numBPLimit; i < vec_bs.size(); i++){
    int pos = vec_random[i];
    delete vec_bs[pos];
    }*/

  vec_bs.clear(); //we don't care about this array anymore -- now all we care about is vec_bs_in
  ///////////////
  
  cout << "# of bipartitions considered (in) = " << vec_bs_in.size() << endl;
 
  ofstream fout;
  
  if (outfilename != "output.tre")
        fout.open(outfilename.c_str());
    else
      fout.open("output.tre");
 

  int32 seed2 = time(NULL); //generate random seed (what happens if we generate it every time, instead of once?
  TRandomMersenne rg2(seed2);

  //beging building trees
  unsigned printMod = 1;
  if (numOutputTrees > 10)
    printMod = numOutputTrees/10;
  short percent = 0;
  cout << "Building trees..." << endl;

  ////////////////// 2 ///////////////////////////////////
  //
  // Construct and resolve each tree to find newly
  // resolved bipartitions
  // If there is no conflict between the newly resolved
  // bipartition and the already collected BPs in vec_bs_in
  // then insert the bipartition in the BP list of the tree
  //
  ////////////////// 2 ///////////////////////////////////

  for (unsigned numOut=0; numOut<numOutputTrees; ++numOut) {
    if (numOut % printMod == 0){
      cout << "Building tree " << numOut << "(" << percent << "% done)" << endl;
      percent+=10;
    }
    
    string temptree = buildTree(vec_bs_in, lm, rg2);
    fout << temptree;
    fout << endl;

  } //end loop
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




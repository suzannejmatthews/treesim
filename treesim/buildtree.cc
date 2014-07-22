//*****************************************************************/
/*
This is TreeZip 2.0, a compression software for phylogenetic trees. 
It is based on HashRF and HashCS, developed by SeungJin Sul

(c) 2011 TreeZip 2.0: Suzanne Matthews
(c) 2009 HashRF : SeungJin Sul 
(c) 2009 HashCS : SeungJin Sul

This file is part of TreeZip.

TreeZip is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TreeZip is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TreeZip.  If not, see <http://www.gnu.org/licenses/>.

*/
/*****************************************************/

#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <valarray>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits>
#include <bitset>

#include "label-map.hh"
#include "SCTree.h"


// For newick parser
extern "C" {
#include <newick.h>
}

using namespace std;

/*bool mmap_sorter(const pair<unsigned unsigned> & a, const pair<unsigned, unsigned> &b) { 
  return a.first > b.first;
  }*/

bool vvec_sorter(const vector<SCNode *> & a, const vector<SCNode*> &b) { 
  return a.size() > b.size();
}

string itos(int i)	
{
  stringstream s;
  s << i;
  return s.str();
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

void add_internal_node(SCNode * node, SCTree *sctree, string label, vector<int> vec_r, unsigned int numChildren, vector<SCNode *> &vec_garbageCan){
  SCNode* intNodeA = new SCNode();
  vec_garbageCan.push_back(intNodeA);
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

//modified from TreeZip 2.0 (modified code from HashCS (c) SeungJin Sul)
string compute_tree(
    LabelMap lm,
    vector< bool * > my_bs,
    unsigned int NUM_TAXA){

  vector<vector<SCNode*> > vvec_distinctClusters2;
  vector<SCNode *> vec_garbageCan;
  multimap<unsigned, unsigned, greater<unsigned> > mmap_cluster;

  //update distinct clusters
  for (unsigned int i = 0; i < my_bs.size(); ++i) {
    vector<SCNode*> vec_nodes2;
    unsigned int lmIndex = 0;
    for (unsigned int j = 0; j < NUM_TAXA; j++) {
      if (my_bs[i][j]) {
	SCNode* aNode = new SCNode();
	aNode->name = lm.name(lmIndex);
	vec_nodes2.push_back(aNode);
	vec_garbageCan.push_back(aNode);
      }
      lmIndex++;
    }
    vvec_distinctClusters2.push_back(vec_nodes2);
  }
    
  //sort (different from TreeZip, where bipartitions are pre-sorted)
  //from hashcs
  for (unsigned int i =0; i < vvec_distinctClusters2.size(); ++i)
    mmap_cluster.insert(multimap<unsigned, unsigned>::value_type(vvec_distinctClusters2[i].size(),i));


  multimap<unsigned, unsigned>::iterator itr;
  SCTree *scTree = new SCTree();
  bool addedRoot = false;
  unsigned int intNodeNum = 0;
  
  for (itr=mmap_cluster.begin(); itr!=mmap_cluster.end(); ++itr) {
    if (!addedRoot) {
      // The first cluster has all the taxa.
      // This constructs a star tree with all the taxa.
      // 1. Dangle all the taxa as root's children by adjusting parent link.
      // 2. Push all the node* in the root's children
      // 3. Push all the nodes in the tree's nodelist.
      // 4. Push all the nodes' parent in the tree's parentlist.
      
      for (unsigned int i=0; i < vvec_distinctClusters2[itr->second].size(); ++i) {
	vvec_distinctClusters2[itr->second][i]->parent = scTree->root;
	scTree->root->children.push_back(vvec_distinctClusters2[itr->second][i]);
	scTree->nodelist.push_back(vvec_distinctClusters2[itr->second][i]);
	assert(scTree->nodelist[0]->name == "root");
	scTree->parentlist2.insert(map<string,int>::value_type(vvec_distinctClusters2[itr->second][i]->name, 0));
      }
      addedRoot = true;
    }
    else {
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
      //			string newIntNodeName = "int" + itostr(intNodeNum, 10);
      string newIntNodeName = "int" + itos(intNodeNum);
      SCNode* newIntNode = new SCNode();
      vec_garbageCan.push_back(newIntNode);
      newIntNode->name = newIntNodeName;
      newIntNode->parent = theParent;
      
      // 3. --------------------------------------------------------------------------
      assert(newIntNodeName.size() != 0);
      scTree->nodelist.push_back(newIntNode);
      assert(scTree->nodelist[scTree->nodelist.size()-1]->name == newIntNode->name);
      
      scTree->parentlist2.insert(map<string, unsigned>::value_type(newIntNodeName, scTree->nodelist.size()-1));
      
      for (unsigned int i=0; i<vvec_distinctClusters2[itr->second].size(); ++i) {
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
  //cout << "tree before resolving:" << scTree->GetTreeString(true) << endl;  
  resolve_tree(scTree->root, scTree, vec_garbageCan);
  string mytree;
  mytree = scTree->GetTreeString(true);

 
 //clean up 
  for (unsigned i=0; i<vec_garbageCan.size(); ++i) {
    if (vec_garbageCan[i]) {
      SCNode* temp = vec_garbageCan[i];
      delete temp;
      vec_garbageCan[i] = NULL;
    }
  }
  
  if (scTree->root) {
    SCNode * tmp = scTree->root;
    delete tmp;
    scTree->root = NULL;
  }
  scTree->nodelist.clear();
  scTree->parentlist2.clear();
  delete scTree;

  return mytree;
}

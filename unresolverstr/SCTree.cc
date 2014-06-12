///*****************************************************/
//
// Copyright (C) 2006, 2007 Seung-Jin Sul
//      Department of Computer Science
//      Texas A&M University
//      Contact: sulsj@cs.tamu.edu
//
//      CLASS DEFINITION
//      SCTree
//
///*****************************************************/

#include "SCTree.hh"
#include <iostream>

SCTree::SCTree()
{
    // Create a root node for the tree
    root = new SCNode();
    root->name = "root";
    root->parent = NULL;
    nodelist.push_back(root);
}


SCTree::SCTree(string name)
{
    // Create a root node for the tree
    root = new SCNode();
    root->name = name;
    root->parent = NULL;
    nodelist.push_back(root);
}


SCTree::~SCTree()
{
    // Delete Nodes in the Tree
    for (unsigned i = 0; i < nodelist.size(); ++i) {
        if (nodelist[i] != NULL) {
            delete nodelist[i];
            nodelist[i] = NULL;
        }
    }
}


//! Rapid display on stdout of the tree structure
void
SCTree::DrawOnTerminal(bool distances)
{
    cout << "\n\nTerminal Representation of the Tree:\n";
    root->DrawOnTerminal(0, distances);
}


//void
//SCTree::PrintParentlist()
//{
//  map<string, SCNode*>::const_iterator itr;
//  for (itr=parentlist.begin(); itr!=parentlist.end(); ++itr)
//      cout << itr->first << " | " << itr->second->name << "(" << itr->second << ")" << endl;
//
//}

//SCNode*
//SCTree::FindParent(string name)
//{
////    map<string, SCNode*>::const_iterator itr;
////    itr=parentlist.find(name);
//
////    if (itr == parentlist.end())  {
////        cout << "Fatal error: taxon name does not exist in the parentlist\n";
////        exit(0);
////    }
////    return itr->second;
//  return parentlist[name];
//}

unsigned
SCTree::FindParent2(string name)
{
//  map<string, SCNode*>::const_iterator itr;
//  itr=parentlist.find(name);

//  if (itr == parentlist.end())  {
//      cout << "Fatal error: taxon name does not exist in the parentlist\n";
//      exit(0);
//  }
//  return itr->second;
    return parentlist2[name];
}


string
SCTree::GetTreeString(
    bool distances,
    double scaleFactor)
{
    string tree;

    GetTreeRecurse(tree, root, distances, scaleFactor);
//  GetTreeRecurse(tree, root);

    return tree + ";";
}


void
SCTree::GetTreeRecurse(
    string& ret,
    SCNode* node,
    bool distances,
    double scaleFactor)
{
    string distance = "";

    if (distances) {
      //double bl = node->GetDistance() / scaleFactor;
      distance="0.123456";
    }

    distance="0.123456";

    if (node->IsLeaf()) {
      ret.append(node->name + ":" + distance);
      //	ret.append(node->name);
        return;
    }

    ret += "(";

// CHECK THIS

    unsigned numChildren = node->NumChildren();

    if (numChildren == 2) {
        if (node->children[0]->name > node->children[1]->name) {
            SCNode *temp = node->children[0];
            node->children[0] = node->children[1];
            node->children[1] = temp;
        }

        if (!node->children[0]->IsLeaf() && !node->children[1]->IsLeaf()) {
            SCNode *temp, *temp2;
            temp = GetLeastSubtree(node);
            temp2 = temp->parent;

            while (temp2 != node) {
                temp = temp2;
                temp2 = temp2->parent;
            }

            if (temp->name == node->children[1]->name) {
                node->children[1] = node->children[0];
                node->children[0] = temp;
            }
        }
    }

    for (unsigned i = 0; i < numChildren; i++) {
        GetTreeRecurse(ret, node->children[i], distances, scaleFactor);
//      GetTreeRecurse(ret, node->children[i]);
        if (i != numChildren - 1)
            ret += ",";
    }

    if (node->IsRoot()) {
        ret += ")";
        return;
    }

    ret = ret + ")";

    // Output node support as
    // internal node label
    // Do it reguardless of distances, because some consense
    // trees might have same tree strings
    // as other trees, we want to distinguish
    // between them in covSEARCH
    if (node->support) {
        char buf[16];
        sprintf(buf, "%d", node->support);
        ret = ret + buf;
    }

    ret = ret + ":" + distance;
//  ret = ret;
}


SCNode*
SCTree::GetLeastSubtree(SCNode *node)
{
    if (node->IsLeaf())
        return node;

    SCNode* temp = node->children[0];
    bool allInts = true;
    unsigned numChildren = node->NumChildren();

    for (unsigned i = 0; i < numChildren; i++) {
        if (node->children[i]->IsLeaf()) {
            allInts = false;
            if (!temp->IsLeaf() || (temp->IsLeaf() && node->children[i]->name < temp->name))
                temp = node->children[i];
        }
    }

    if (allInts) {
        vector<SCNode*> kids;
        kids.resize(numChildren);

        for (unsigned i = 0; i < numChildren; i++)
            kids[i] = GetLeastSubtree(node->children[i]);

        temp = kids[0];
        for (unsigned i = 0 ;  i < numChildren; i++)
            if (kids[i]->name < temp->name)
                temp = kids[i];
    }

    return temp;
}

void
SCTree::DeleteAllNodes()
{
    for (unsigned i=0; i<nodelist.size(); ++i) {
        if (nodelist[i] != NULL) {
            delete nodelist[i];
            nodelist[i] = NULL;
        }
    }
}


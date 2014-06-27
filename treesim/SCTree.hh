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

#ifndef __SCTree_hh__
#define __SCTree_hh__

#include  "SCNode.hh"

class SCTree {

public:
    SCTree();
    SCTree(string name);
    ~SCTree();

    SCNode* root;
    vector<SCNode*> nodelist;
//  vector<SCNode*> parentlist;
//  map<string, SCNode*> parentlist;
    map<string, unsigned> parentlist2;

    void DrawOnTerminal(bool distances = false);
    void PrintParentlist();
//  SCNode* FindParent(string name);
    unsigned FindParent2(string name);

    string GetTreeString(bool distances = false, double scaleFactor = 100.0);
    void GetTreeRecurse(string& ret, SCNode* node, bool distances, double scaleFactor);
//  void GetTreeRecurse(string& ret, SCNode* node);
    SCNode* GetLeastSubtree(SCNode *node);
    void DeleteAllNodes();
};





#endif



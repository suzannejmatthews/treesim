///*****************************************************/
//
// Copyright (C) 2006, 2007 Seung-Jin Sul
//      Department of Computer Science
//      Texas A&M University
//      Contact: sulsj@cs.tamu.edu
//
//      CLASS DEFINITION
//      HashMap: Class for hashing
//
///*****************************************************/

#ifndef __SCNode_hh__
#define __SCNode_hh__


using namespace std;

#include <vector>
#include <map>
#include <fstream>

class SCNode {

public:

    SCNode() ;
    ~SCNode();

    string name;
    SCNode *parent;
    vector<SCNode*> children;


    // obsolete
    unsigned support;
    double bl;

    bool IsLeaf();
    bool IsRoot();
    void SetDistance(double distance);
    double GetDistance() const;
    void ClearChildren();

    unsigned NumChildren();
//  string name() { return name; };
    void DrawOnTerminal(int dp, bool distances = false);
//  void SCDeleteChild(unsigned i);
//  void SCAddChild(SCNode* node);
//  void SCClearChildren();

};


#endif



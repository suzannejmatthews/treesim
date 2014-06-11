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

#include "SCTree.hh"
#include "SCNode.hh"
#include <iostream>
#include <cassert>

SCNode::SCNode()
{
    parent = NULL;
    ClearChildren();

    support = 0;
}


SCNode::~SCNode()
{
    // Delete Nodes in the Tree
//  for (unsigned i = 0; i < children.size(); ++i) {
//      if (children[i] != NULL) {
//          delete children[i];
//          children[i] = NULL;
//      }
//  }
}


void
SCNode::ClearChildren()
{
    for (unsigned i=0 ; i < children.size(); ++i)
        children[i] = NULL;
}



unsigned
SCNode::NumChildren()
{

    unsigned cnt=0;
    for (unsigned i=0; i<children.size(); ++i)
        if (children[i] != NULL) cnt++;

    return cnt;

}


void
SCNode::SetDistance(double distance)
{
    // One of manu's datasets gives us intermediate negative bl
    if (distance < 0)
        bl = 0.0;
    else
        bl = distance;
}

double
SCNode::GetDistance() const
{
    return bl;
}


void SCNode::DrawOnTerminal(int dp, bool distances)
{
    int i = 0;

    while (i < dp) {
        cout << "    ";
        ++i;
    }

    cout << name;

//  if (support)
//      cout << " (" << support << ")";

//  if (distances)
//      cout << " " << bl;

    cout << endl;

    for (unsigned i = 0; i < NumChildren(); ++i)
        children[i]->DrawOnTerminal(dp + 1, distances);
}

bool
SCNode::IsRoot()
{
    return (parent == NULL);
}

bool
SCNode::IsLeaf()
{
//  assert(children[2] == NULL);
//  return (children[0] == NULL && children[1] == NULL);
    return (NumChildren() == 0);
}


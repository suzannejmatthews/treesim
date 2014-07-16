/*
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



#include "SCTree.h"
#include "SCNode.h"
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
//	for (unsigned int i = 0; i < children.size(); ++i) {
//		if (children[i] != NULL) {		
//			delete children[i];
//			children[i] = NULL;
//		}
//	}
}


void 
SCNode::ClearChildren()
{	
	for (unsigned int i=0 ; i < children.size(); ++i)
		children[i] = NULL;	
}



unsigned int
SCNode::NumChildren()
{
	
	unsigned int cnt=0;
	for (unsigned int i=0; i<children.size(); ++i)
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
	
//	if (support)
//		cout << " (" << support << ")";

	if (distances)
		cout << " " << bl;
	
	cout << endl;

	for (unsigned int i = 0; i < NumChildren(); ++i)
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
	return (NumChildren() == 0);
}


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

#ifndef __SCNode_h__
#define __SCNode_h__

#include <vector>
#include <map>
#include <fstream>

using namespace std;

class SCNode {

public:
		
	SCNode() ;
	~SCNode();
	
	string name;
	SCNode *parent;
	//float dist;
	vector<SCNode*> children;
	
	double bl;

	// obsolete
	unsigned int support;
	

	bool IsLeaf();
	bool IsRoot();
	void SetDistance(double distance);
	double GetDistance() const;
	void ClearChildren();
			
	unsigned int NumChildren();
	void DrawOnTerminal(int dp, bool distances = false);
	
};


#endif





#ifndef __SCTree_h__
#define __SCTree_h__

#include  "SCNode.h"

class SCTree {
	
public:
	SCTree();
	SCTree(string name);
	~SCTree();
	
	SCNode* root;
	vector<SCNode*> nodelist;
	map<string, unsigned int> parentlist2;
	
	void DrawOnTerminal(bool distances);
	void PrintParentlist();
	unsigned int FindParent2(string name);
	
	string GetTreeString(bool distances);
	void GetTreeRecurse(string& ret, SCNode* node, bool distances);
	SCNode* GetLeastSubtree(SCNode *node);
	void DeleteAllNodes();
};





#endif



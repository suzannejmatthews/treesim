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


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <limits>
#include <bitset>
#include <map>
#include <utility>

#include "label-map.hh"
#include "SCTree.h"
using namespace std;

#ifndef _BUILDTREE_H_
#define _BUILDTREE_H_

//treezip's compute tree
/*string compute_tree( 
    LabelMap lm,
    vector< bool* > my_bs,
    vector<float> my_branches,
    unsigned id,
    bool branch);
*/

//additional treesim functions
vector<int> pickChildren(unsigned int pick, unsigned int numChildren);

void add_internal_node(SCNode *node, SCTree *sctree, string label, vector<int> vec_r, unsigned int numChildren, vector<SCNode *> & vec_garbageCan);

void resolve_tree(SCNode * node, SCTree *sctree, vector<SCNode *> & trashcan);

//treesim's compute tree
string compute_tree( 
    LabelMap lm,
    vector< bool* > my_bs,
    unsigned int NUM_TAXA);



#endif 

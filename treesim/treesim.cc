/**********************************************************************
treesim
**********************************************************************

    FILE: treesim.cc
   
    DESCRIPTION: Creates "accurate" simluated datasets based on a given strict
    and majority consensus rate.

    Based on the software unresolverstr and unresolvermaj by SeungJin Sul
    and TreeZip by Suzanne J. Matthews

    AUTHOR:
        Copyright (C) 2014 Suzanne J. Matthews
        Dept. of Electrical Engineering & Computer Science
        United States Military Academy
        West Point, NY, 10996, USA
        (contact: suzanne.matthews@usma.edu)
       
    HISTORY:
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


#include "SCTree.h"
#include "SCNode.h"

// From Split-Dist
#include "label-map.hh"
#include "buildtree.h"

#include <cassert>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdint.h>

// For newick parser
extern "C" {
#include <newick.h>
}

#include <iostream>
#include <./tclap/CmdLine.h>

using namespace std;

void GetTaxaLabels2(NEWICKNODE *node, LabelMap &lm) {
  if (node->Nchildren == 0) {
    string temp(node->label);
    lm.push(temp);
  } else
    for(int i=0; i<node->Nchildren; ++i)
      GetTaxaLabels2(node->child[i], lm);
}


bool * collect_biparts(NEWICKNODE* startNode, LabelMap &lm, unsigned treeIdx, vector<bool *> & vec_bs, unsigned int NUM_TAXA)
{
  if (startNode->Nchildren == 0) {
    // leaf node
    bool* bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    string temp(startNode->label);
    unsigned idx = lm[temp];
    bs[idx] = 1;
    return bs;
  } 
  else { //internal node
    bool* bs = new bool[NUM_TAXA];
    for (unsigned int i = 0; i < NUM_TAXA; i++)
      bs[i] = 0;
    for (int i=0; i<startNode->Nchildren; ++i) {
      bool * ebs = collect_biparts(startNode->child[i], lm, treeIdx, vec_bs, NUM_TAXA);
      for (unsigned int j = 0; j < NUM_TAXA; j++)
	bs[j] |= ebs[j];

      //if (ebs) {
      //delete [] ebs;
      //ebs = NULL;
      //}
    }

    vec_bs.push_back(bs);
    
    return bs;
  }
}


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

vector<unsigned int> genRandomNums(unsigned int howmany, unsigned int max){
  vector<unsigned int> vec_rand;
  for ( unsigned i = 0; i < howmany; ++i ) {
    unsigned int ir = rand()%max;
    if (find(vec_rand.begin(), vec_rand.end(), ir) != vec_rand.end()) {
      --i;
      continue;
    }
    vec_rand.push_back(ir);
  }
  return vec_rand;
}

int main(int argc, char** argv)
{
  string outfilename, startingFile;
  float majRate = 0.0, strictRate = 0.0;
  unsigned int NUM_TAXA = 0, NUM_TREES = 0, unique_trees = 0;
  bool verbose = false, precise = false;
    // TCLAP
    try {
        // Define the command line object.
        string    helpMsg  = "treesim\n";

        helpMsg += "Input file: \n";
        helpMsg += "   The current version of treesim only supports the Newick format.\n";
        helpMsg += "Sample Newick file: \n";
        helpMsg += "   (('Chimp':0.052625,'Human':0.042375):0.007875,'Gorilla':0.060125,\n";
        helpMsg += "   ('Gibbon':0.124833,'Orangutan':0.0971667):0.038875);\n";
        helpMsg += "   ('Chimp':0.052625,('Human':0.042375,'Gorilla':0.060125):0.007875,\n";


        helpMsg += "File option: (default = output.tre)\n";
        helpMsg += "   -o <export-file-name>, specify a file name to save the result tree.\n";

        helpMsg += "Example: \n";
        helpMsg += "  treesim 500 1000 0.5 .75 -o out.tre\n";
	helpMsg += "  generates a random 500-taxa, 1000 tree dataset with a strict\n";
	helpMsg += "  consensus rate of at least 50% and a majority consensus rate\n";
	helpMsg += "  of 75%. The result is stored in out.tre\n";

        TCLAP::CmdLine cmd(helpMsg, ' ', "0.1");

        TCLAP::UnlabeledValueArg<int>  numtaxaArg( "numtaxa", "number of taxa", true, 5, "Number of taxa"  );
        cmd.add( numtaxaArg );

        TCLAP::UnlabeledValueArg<unsigned>  numtreeArg( "numtree", "number of output trees", true, 1, "Number of output trees"  );
        cmd.add( numtreeArg );

        TCLAP::UnlabeledValueArg<float>  sRateArg( "sRate", "strict consensus rate", true, 0.0, "strict rate"  );
        cmd.add( sRateArg );

        TCLAP::UnlabeledValueArg<float>  mRateArg( "mRate", "majority consensus rate", true, 0.0, "majority rate"  );
        cmd.add( mRateArg );

        TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "output.tre", "Output file name");
        cmd.add( outfileArg );

	TCLAP::ValueArg<string> sfileArg("t", "startingtree", "starting tree file name", false, "starting.tre", "input starting tree file name");
        cmd.add( sfileArg );

	TCLAP::ValueArg<int> uArg("u", "uniquetree", "number of unique trees", false, 0, "number of unique trees in file");
        cmd.add( uArg );

	TCLAP::SwitchArg vArg("v", "verbose", "print out progress", false);
        cmd.add( vArg );

	TCLAP::SwitchArg pArg("p", "precise", "turn on precise version of algorithm", false);
        cmd.add( pArg );


        cmd.parse( argc, argv );
	NUM_TAXA = numtaxaArg.getValue();
        NUM_TREES = numtreeArg.getValue();
        strictRate = sRateArg.getValue();
	majRate = mRateArg.getValue();
	unique_trees = uArg.getValue();
	verbose = vArg.getValue();
	precise = pArg.getValue();

	if (strictRate > 1 || strictRate < 0){
	  cerr << "ERROR: strict consensus rate must be between 0 and 1!" << endl;
	  return 1;
	}
	if (majRate > 1 || majRate < 0){
	  cerr << "ERROR: strict consensus rate must be between 0 and 1!" << endl;
	  return 1;
	}
	if (strictRate > majRate){
	  cerr << "ERROR: strict consensus rate must be lesser than or equal to majority consensus rate!" << endl;
	  return 1;
	}
        outfilename = outfileArg.getValue();
	startingFile = sfileArg.getValue();
    } catch (TCLAP::ArgException &e) { // catch any exceptions
        cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    }


    //generate a random tree with the specified number of taxa
   LabelMap lm;
   ofstream fout;
   NEWICKTREE *newickTree;
   int err;
   FILE *fp;
    if (startingFile == "starting.tre"){
      cout << "We get here!" << endl;
      fprintf(stderr, "Generating a random tree with %u taxa...\n", NUM_TAXA);
      vector<bool *> random_tree_bs;
      bool * star = new bool[NUM_TAXA];
      for (unsigned int i = 0; i < NUM_TAXA; i++){
	star[i] = 1;
	string taxa;
	stringstream ss;
	ss << i+1;
	ss >> taxa;
	lm.push(taxa);
      }
      random_tree_bs.push_back(star);
      string random_tree = compute_tree(lm, random_tree_bs, NUM_TAXA); 
      fout.open("starting.tre");
      fout << random_tree << endl;
      fout.close();
      fprintf(stderr, "Done. starting tree outputted to starting.tre.\n");
      random_tree_bs.clear();
      delete [] star;
    }
    else{
      //collect label map from file
      cout << "getting labels from file: " << startingFile << endl;
      fp = fopen(startingFile.c_str(), "r");
      if (!fp){
        cout << "ERROR: file open error:" << argv[1] << endl;
        exit(0);
      }
      newickTree = loadnewicktree2(fp, &err);
      if(!newickTree) {
	print_error(err);
      }
      else{
	GetTaxaLabels2(newickTree->root, lm);
	killnewicktree(newickTree);
	fclose(fp);
      }
      assert(NUM_TAXA == lm.size());
    }
    fprintf(stderr, "Building Collection...\n");
 

    vector<bool *> vec_bs;
    fp = fopen(startingFile.c_str(), "r");
    if(!fp) {
        cout << "ERROR: file open error:" << argv[1] << endl;
        exit(0);
    }

    newickTree = loadnewicktree2(fp, &err);
    if(!newickTree) {
      print_error(err);
    }
    else{
      collect_biparts(newickTree->root, lm, 0, vec_bs, NUM_TAXA);
      killnewicktree(newickTree);
    }

    fclose(fp);

    unsigned int total_BPs = vec_bs.size()-1;
    unsigned int majority_resolution_rate = int(total_BPs*majRate);
    unsigned int strict_resolution_rate = int(total_BPs*strictRate);
    unsigned int difference = majority_resolution_rate - strict_resolution_rate;
    cout << "    vec_bs.size() = " << vec_bs.size()-1 << endl;
    cout << "    Number of Output trees = " << NUM_TREES << endl;
    cout << "    Number of bipartitions that will be in all the trees = " << strict_resolution_rate << endl;
    cout << "    Number of bipartitions that will be additionally in 50% or more of trees  = " << difference << endl;
    
    srand(time(NULL)); //seed random number generator

    //select r% of the bipartitions at random
    vector<unsigned int> vec_random = genRandomNums(majority_resolution_rate, total_BPs);
    bool * select;
    vector<unsigned int> other_biparts;
    unsigned int place;
    if (precise){
      select = (bool*)calloc(total_BPs,sizeof(bool));
      for (unsigned int i = 0; i < vec_random.size(); i++){
	place = vec_random[i];
	select[place]=1;
      }
      for (unsigned int i = 0; i < total_BPs; i++){
	if (select[i] == 0)
	  other_biparts.push_back(i);
      }
      free(select);
    }
    cout << "Requested Majority Rate=" << majRate << endl;
    cout << "Attempted Majority Rate=" << float(vec_random.size())/total_BPs << endl;    
    cout << "Requested Strict Rate=" << strictRate << endl;
    cout << "Attempted Strict Rate=" << float(strict_resolution_rate)/total_BPs << endl;    
    //generate t vectors of vectors (n x t)
    vector< vector<bool*> > tree_matrix;
    vector<unsigned int> duplicates;
    unsigned int NUM_TO_BUILD = 0;
    if (unique_trees == 0){ //if this parameter is not specfied
      NUM_TO_BUILD = NUM_TREES;
    }
    else{
      assert(unique_trees < NUM_TREES);
      assert(unique_trees > 0);
      NUM_TO_BUILD = unique_trees;
      duplicates = genRandomNums(NUM_TREES-unique_trees, NUM_TREES);
      sort(duplicates.begin(), duplicates.end()); //sort the duplicates vector
    }
    //tree_matrix.resize(NUM_TREES);
    tree_matrix.resize(NUM_TO_BUILD);
    for (unsigned int i = 0; i < strict_resolution_rate; i++){
      //cout << "adding bipartition " << vec_random[i] << "to all the trees" << endl;
      unsigned int bipart = vec_random[i];
      for (unsigned int j = 0; j < NUM_TO_BUILD; j++){
	tree_matrix[j].push_back(vec_bs[bipart]);
      }
    }

    //for each of the remainder of the selected bipartitions:
    for (unsigned int i = strict_resolution_rate; i < vec_random.size(); i++){
      unsigned int bipart = vec_random[i];
      unsigned int perc = rand() % 51 + 50; // generate a random number between 50 .. 100
      //cout << "Adding bipartition: " << vec_random[i] << " to " << perc << "% of the trees" << endl;
      perc = int((float(perc)/100) * NUM_TO_BUILD);
      //cout << "perc will actually be: " << perc << endl;
      vector<unsigned int> selected_trees = genRandomNums(perc, NUM_TO_BUILD);
      for (unsigned int j = 0; j < selected_trees.size(); j++){
	unsigned int tree_id = selected_trees[j];
	tree_matrix[tree_id].push_back(vec_bs[bipart]);
	//cout << "added bipartition " << vec_random[i] << "to tree " << tree_id << endl;
      }
    }

    //add the star bipartition to all the trees
    for (unsigned int i = 0; i < NUM_TO_BUILD; i++)
      tree_matrix[i].push_back(vec_bs[total_BPs]);

    //commence building

    if (outfilename != "output.tre")
      fout.open(outfilename.c_str());
    else
      fout.open("output.tre");
    
    cout << "building trees!" << endl;
    unsigned int dupCount = 0;
    for (unsigned numOut=0; numOut<NUM_TO_BUILD; ++numOut) {
      if (verbose && numOut % 1000 == 0)
	cout << numOut << endl;
      string tree = compute_tree(lm, tree_matrix[numOut], NUM_TAXA);
      fout << tree << endl;
      if (unique_trees != 0 && numOut == duplicates[dupCount]){
	fout << tree << endl;
	dupCount++;
      }
      //tree_matrix[numOut].clear(); //remove the bipartitions from this
    }
    fout.close();

    //clean up
    for (unsigned int i = 0; i < vec_bs.size(); i++){
      if (vec_bs[i]!=NULL){
	delete [] vec_bs[i];
	vec_bs[i] = NULL;
      }
    }
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

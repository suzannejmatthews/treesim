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

void PopulateTaxaLabels( NEWICKNODE *node, LabelMap &lm) {
  if (node->Nchildren == 0) { //if leaf node
    string temp(node->label);
    lm.add(temp);
  }
  else
    for (int i=0;i<node->Nchildren;++i)
      PopulateTaxaLabels(node->child[i], lm);
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
  if (howmany < max){
    for ( unsigned i = 0; i < howmany; ++i ) {
      unsigned int ir = rand()%max;
      if (find(vec_rand.begin(), vec_rand.end(), ir) != vec_rand.end()) {
	--i;
	continue;
      }
      vec_rand.push_back(ir);
    }
  }
  else{
    for (unsigned i = 0; i < howmany; i++){
      unsigned int ir= rand()% max;
      vec_rand.push_back(ir);
    }
  }
  return vec_rand;
}

void simulate_from_random(string startingFile, unsigned int NUM_TAXA, unsigned int NUM_TREES, LabelMap lm, float strictRate, float majRate, bool precise, unsigned int unique_trees, unsigned int DUPLICATES, bool weighted, string outfilename, bool verbose){
  //this code contains the original treesim code, where we were simulating trees 
  //from purely the strict and majority consensus rate. The starting tree is always 
  //random. 
  NEWICKTREE *newickTree;
  ofstream fout;
  int err;
  FILE *fp;
  fprintf(stderr, "Building Collection...\n");
  vector<bool *> vec_bs;
  fp = fopen(startingFile.c_str(), "r");
  if(!fp) {
    cout << "ERROR: file open error:" << startingFile << endl;
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
  //unsigned int algType = 0;
  if (majRate == strictRate && majRate == 0){
    cout << "Algorithm: Random Collection" << endl;
    //algType = 0;
  }
  if (majRate == strictRate && majRate > 0){
    cout << "Algorithm: Strict Consensus Collection" << endl;
    //algType = 1;
  }
  if (majRate > 0 && strictRate == 0){
    cout << "Algorithm: Majority Consensus Collection" << endl;
    if (precise)
      cout << "Precise option: ON" << endl;
    //algType = 2;
  }
  if (majRate > 0 && strictRate > 0 && majRate != strictRate){
    cout << "Algorithm: Combined Consensus Collection" << endl;
    if (precise)
      cout << "Precise option: ON" << endl;
    //algType = 3;
  }
  
  unsigned int total_BPs = vec_bs.size()-1;
  unsigned int majority_resolution_rate = int(total_BPs*majRate);
  unsigned int strict_resolution_rate = int(total_BPs*strictRate);
  unsigned int difference = majority_resolution_rate - strict_resolution_rate;
  vector<unsigned int> duplicates, vec_random, other_biparts;
  unsigned int NUM_TO_BUILD = 0;
  vector< vector<bool*> > tree_matrix;
  cout << "    vec_bs.size() = " << vec_bs.size()-1 << endl;
  cout << "    Number of Output trees = " << NUM_TREES << endl;
  cout << "    Number of bipartitions that will be in all the trees = " << strict_resolution_rate << endl;
  cout << "    Number of bipartitions that will be additionally in 50% or more of trees  = " << difference << endl;
  
  srand(time(NULL)); //seed random number generator
  
  //select r% of the bipartitions at random
  vec_random = genRandomNums(majority_resolution_rate, total_BPs);
  bool * select;
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
  
  //generate duplicates, if needed
  
  if (unique_trees == 0){ //if this parameter is not specfied
    NUM_TO_BUILD = NUM_TREES;
  }
  else{
    assert(unique_trees < NUM_TREES);
    assert(unique_trees > 0);
    NUM_TO_BUILD = unique_trees;
    DUPLICATES = NUM_TREES - unique_trees;
    duplicates = genRandomNums(DUPLICATES, NUM_TO_BUILD);
    sort(duplicates.begin(), duplicates.end()); //sort the duplicates vector
  }
  
  //generate t vectors of vectors (n x t)
  tree_matrix.resize(NUM_TO_BUILD);
  for (unsigned int i = 0; i < strict_resolution_rate; i++){
    //cout << "adding bipartition " << vec_random[i] << "to all the trees" << endl;
    unsigned int bipart = vec_random[i];
    for (unsigned int j = 0; j < NUM_TO_BUILD; j++){
      tree_matrix[j].push_back(vec_bs[bipart]);
    }
  }
  
  //for each of the remainder of the selected bipartitions: (this code doesn't execute when strict_resolution rate = majority_resolution_rate
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
  if (precise){
    for (unsigned int i = 0; i < other_biparts.size(); i++){
      unsigned int bipart = other_biparts[i];
      unsigned int perc = rand() % 50; //generate a random number between 0 .. 49
      perc = int(float(perc)/100)* NUM_TO_BUILD;
      vector<unsigned int> selected_trees = genRandomNums(perc, NUM_TO_BUILD);
      for (unsigned int j = 0; j < selected_trees.size(); j++){
	unsigned int tree_id = selected_trees[j];
	tree_matrix[tree_id].push_back(vec_bs[bipart]);
      }
    }
  }
  
  //add the star bipartition to all the trees
  for (unsigned int i = 0; i < NUM_TO_BUILD; i++)
    tree_matrix[i].push_back(vec_bs[total_BPs]);
  
  //commence building
  fout.open(outfilename.c_str());

  cout << "building trees!" << endl;
  unsigned int dupCount = 0;
  assert(NUM_TO_BUILD+DUPLICATES == NUM_TREES);
  for (unsigned numOut=0; numOut<NUM_TO_BUILD; numOut++) {
    if (verbose && numOut % 1000 == 0)
      cout << numOut << endl;
    string tree = compute_tree(lm, tree_matrix[numOut], NUM_TAXA, weighted);
    fout << tree << endl;
    if (unique_trees != 0 && dupCount < DUPLICATES){
      while (numOut == duplicates[dupCount]){
	fout << tree << endl;
	dupCount++;
      }
    }
    //tree_matrix[numOut].clear(); //remove the bipartitions from this
  }
  fout.close();
  assert(dupCount == DUPLICATES);
  //clean up
  for (unsigned int i = 0; i < vec_bs.size(); i++){
    if (vec_bs[i]!=NULL){
      delete [] vec_bs[i];
      vec_bs[i] = NULL;
    }
  }
  
}

void simulate_from_strict(string strictFile, unsigned int NUM_TAXA, unsigned int NUM_TREES, bool weighted, string outfilename, bool verbose){
  NEWICKTREE *newickTree;
  ofstream fout;
  int err;
  FILE *fp;
  vector<bool *> vec_bs;
  LabelMap lm;
  //collect labels
  cout << "Collecting labels from " << strictFile << endl;
  fp = fopen(strictFile.c_str(), "r");
  if(!fp) {
    cerr << "ERROR: file open error:" << strictFile << endl;
    exit(0);
  }  
  newickTree = loadnewicktree2(fp, &err);
  if(!newickTree) {
    print_error(err);
  }
  else{
    PopulateTaxaLabels(newickTree->root, lm);
    killnewicktree(newickTree);
  }
  fclose(fp);
  cout << "Found: " << lm.size() << " taxa" << endl;
  assert(lm.size() == NUM_TAXA);

  //collect strict bipartitions from the strict consensus tree
  cout << "Collecting strict consensus bipartitions from " << strictFile << endl;
  fp = fopen(strictFile.c_str(), "r");
  if(!fp) {
    cerr << "ERROR: file open error:" << strictFile << endl;
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
  cout << "done." << endl;

  cerr << "**** PRINTING OUT BIPARTITIONS***" << endl;
  unsigned int count = 0;
  for (unsigned int i = vec_bs.size()-2; i < vec_bs.size(); i++){
    for (unsigned int j = 0; j < NUM_TAXA; j++)
      count+=vec_bs[i][j];
  }
  if (count == 2*NUM_TAXA)
    vec_bs.pop_back();

  unsigned int total_BPs = vec_bs.size()-1;
  float strict_rate = (float)total_BPs/(NUM_TAXA-3);
  cout << "    vec_bs.size() = " << vec_bs.size()-1 << endl;
  cout << "    Number of Output trees = " << NUM_TREES << endl;
  cout << "    Number of bipartitions that will be in all the trees = " << total_BPs << endl;
  cout << "Detected Strict Rate=" << strict_rate << endl;

  //we currently don't handle duplicates -- right now, we're just trying to implement the simple building 
  fout.open(outfilename.c_str());
  for (unsigned numOut=0; numOut<NUM_TREES; numOut++) {
    if (verbose && numOut % 1000 == 0)
      cout << numOut << endl;
    string tree = compute_tree(lm, vec_bs, NUM_TAXA, weighted);
    fout << tree << endl;
    //tree_matrix[numOut].clear(); //remove the bipartitions from this
  }
  fout.close();
  //delete [] star;
  //clean up
  for (unsigned int i = 0; i < vec_bs.size(); i++){
    if (vec_bs[i]!=NULL){
      delete [] vec_bs[i];
      vec_bs[i] = NULL;
    }
  }
}

void simulate_from_majority(string majFile, unsigned int NUM_TAXA, unsigned int NUM_TREES, bool weighted, string outfilename, bool verbose){
  NEWICKTREE *newickTree;
  ofstream fout;
  int err;
  FILE *fp;
  vector<bool *> vec_bs;
  LabelMap lm;
  //collect labels
  cout << "Collecting labels from " << majFile << endl;
  fp = fopen(majFile.c_str(), "r");
  if(!fp) {
    cerr << "ERROR: file open error:" << majFile << endl;
    exit(0);
  }  
  newickTree = loadnewicktree2(fp, &err);
  if(!newickTree) {
    print_error(err);
  }
  else{
    PopulateTaxaLabels(newickTree->root, lm);
    killnewicktree(newickTree);
  }
  fclose(fp);

  //collect strict bipartitions from the strict consensus tree
  cout << "Collecting majority consensus bipartitions from tree located in " << majFile << endl;
  fp = fopen(majFile.c_str(), "r");
  if(!fp) {
    cerr << "ERROR: file open error:" << majFile << endl;
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
  float maj_rate = (float)total_BPs/(NUM_TAXA-3);
  cout << "    vec_bs.size() = " << vec_bs.size()-1 << endl;
  cout << "    Number of Output trees = " << NUM_TREES << endl;
  cout << "    Number of bipartitions that will be in a majority of all the trees = " << total_BPs << endl;
  cout << "Detected Majority Rate=" << maj_rate << endl;

}

int main(int argc, char** argv)
{
  string outfilename, strictFile, majFile;
  float majRate = 0.0, strictRate = 0.0;
  unsigned int NUM_TAXA = 0, NUM_TREES = 0, unique_trees = 0, DUPLICATES=0;
  bool verbose = false, precise = false, weighted=false;
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

    TCLAP::UnlabeledValueArg<float>  sRateArg( "sRate", "strict consensus rate", true, 0.0, "strict rate");
    cmd.add( sRateArg );
    
    TCLAP::UnlabeledValueArg<float>  mRateArg( "mRate", "majority consensus rate", true, 0.0, "majority rate");
    cmd.add( mRateArg );
    
    TCLAP::ValueArg<string> outfileArg("o", "outfile", "Output file name", false, "output.tre", "Output file name");
    cmd.add( outfileArg );
    
    TCLAP::ValueArg<string> mfileArg("m", "majFile", "majority consensus tree", false, "maj.tre", "majority consensus tree file");
    cmd.add( mfileArg );
    
    TCLAP::ValueArg<string> sfileArg("s", "strictFile", "strict consensus tree", false, "strict.tre", "strict consensus tree file");
    cmd.add( sfileArg );
    
    TCLAP::ValueArg<int> uArg("u", "uniquetree", "number of unique trees", false, 0, "number of unique trees in file");
    cmd.add( uArg );
    
    TCLAP::SwitchArg vArg("v", "verbose", "print out progress", false);
    cmd.add( vArg );
    
    TCLAP::SwitchArg wArg("w", "weighted", "create weighted trees", false);
    cmd.add( wArg );
    
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
    weighted = wArg.getValue(); 
    outfilename = outfileArg.getValue();
    strictFile = sfileArg.getValue();
    majFile = mfileArg.getValue();   
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

  } catch (TCLAP::ArgException &e) { // catch any exceptions
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  
  //generate a random tree with the specified number of taxa
  LabelMap lm;
  ofstream fout;
  //NEWICKTREE *newickTree;
  //int err;
  //FILE *fp;
  if (strictFile == "strict.tre" && majFile == "maj.tre"){
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
    string random_tree = compute_tree(lm, random_tree_bs, NUM_TAXA, weighted); 
    fout.open("starting.tre");
    fout << random_tree << endl;
    fout.close();
    fprintf(stderr, "Done. starting tree outputted to starting.tre.\n");
    random_tree_bs.clear();
    delete [] star;
    simulate_from_random("starting.tre", NUM_TAXA, NUM_TREES, lm, strictRate, majRate, precise, unique_trees, DUPLICATES, weighted, outfilename, verbose);
  }
  else{
    if (strictFile != "strict.tre" && majFile == "maj.tre"){
      cout << "Algorithm: Strict Consensus Resolution (From Consensus Input)" << strictFile << "!" << endl;
      simulate_from_strict(strictFile, NUM_TAXA, NUM_TREES, weighted, outfilename, verbose);
      
    }
    else if (strictFile == "strict.tre" && majFile != "maj.tre"){
      cout << "Algorithm: Majority Consensus Resolution (From Consensus Input)" << majFile << "!" << endl;
    }
    else{
      cout << "we will construct the combined consensus resolution tree using the trees located in " << strictFile << " and " << majFile << "!" << endl;
    }
  }
  
  // CPU time consumed
  struct rusage a;
  if (getrusage(RUSAGE_SELF,&a) == -1) {
    cerr << "ERROR: getrusage failed.\n";
    exit(2);
  }
  cout << "\n    Total CPU time: " << a.ru_utime.tv_sec+a.ru_stime.tv_sec << " sec and ";
  cout << a.ru_utime.tv_usec+a.ru_stime.tv_usec << " usec.\n";
  
  return 0;
}

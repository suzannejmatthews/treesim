// Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
//                             University of Aarhus, Denmark
//                             Contact: Thomas Mailund <mailund@birc.dk>

#include "label-map.hh"
#include <iostream>
#include <algorithm>
using namespace std;

size_t
LabelMap::push(string label) throw(AlreadyPushedEx)
{
    if (_map.find(label) != _map.end()) throw AlreadyPushedEx(label);
    _names.push_back(label);  
    return _map[label] = _count++;
}

size_t
LabelMap::add(string label) {
  if (_map.find(label) == _map.end()){
    _names.push_back(label);  
    return _map[label] = _count++;
  }
  else 
    return _map[label];
}

size_t
LabelMap::operator[](string label) throw(UnkownLabelEx)
{
    map<string, size_t>::const_iterator i = _map.find(label);
    if (i == _map.end()) {
      std::cout << "cannot find label : " << label << std::endl; 
      throw UnkownLabelEx(label);
    }
    return i->second;
}

std::string
LabelMap::name(unsigned int idx) const throw(std::out_of_range)
{
    return _names.at(idx);
}

void 
LabelMap::rename(string currentl, string newl, unsigned int loc)
{
  _names[loc] = newl;
  _map[newl] = _map[currentl];
}

void
LabelMap::sortTaxa(){
  _names.clear();
  size_t pos = 0;
  string label;
  map<string, size_t>::iterator i = _map.begin();
  while (i  != _map.end()){
    label = i->first;
    i->second = pos;
    _names.push_back(label);
    pos++;
    i++;
  }
}

void 
LabelMap::printMap(){
  map<string, size_t>::iterator i = _map.begin();
  while (i  != _map.end()){
    cout << i->first << ": " << i->second << endl;
    i++;
  }
  cout << "labels are:" << endl;
  for (unsigned int i =0; i < _names.size(); i++)
    cout << "_names[" << i << "]: " << _names[i] << endl;
}

void 
LabelMap::clear(){
  _names.clear();
  _map.clear();
  _count = 0;
}



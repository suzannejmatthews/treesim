// Copyright (C) 2003, 2004 by BiRC -- Bioinformatics Research Center
//                             University of Aarhus, Denmark
//                             Contact: Thomas Mailund <mailund@birc.dk>

#ifndef LABEL_MAP_HH
#define LABEL_MAP_HH

#include <map>
#include <vector>
#include <stdexcept>
#include <string>

class LabelMap {
  size_t _count;
  std::map<std::string, size_t> _map;
  std::vector<std::string>  _names;

public:
  LabelMap() : _count(0) {};

  struct AlreadyPushedEx {
    std::string label;
    AlreadyPushedEx(std::string l) : label(l) {}
  };
    
  struct UnkownLabelEx {
    std::string label;
    UnkownLabelEx(std::string l) : label(l) {}
  };

  size_t push(std::string label) throw(AlreadyPushedEx);
size_t push2(std::string label);
  size_t add(std::string label);
  void rename(std::string oldl, std::string newl, unsigned int loc);
  size_t size() const { return _count; }
  void sortTaxa();
  void clear();
  void printMap();
  size_t operator[](std::string label) throw(UnkownLabelEx);
  std::string name(unsigned int idx) const throw(std::out_of_range);
};

#endif // LABEL_MAP_HH

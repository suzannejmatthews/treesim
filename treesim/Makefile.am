
AM_CXXFLAGS = -static -Wall
#CXXFLAGS += -O3 -pg
AM_CXXFLAGS += -O3 -fexpensive-optimizations

#AM_CXXFLAGS += -DNDEBUG

# misc.

# real program stuff
bin_PROGRAMS = treesim

treesim_SOURCES =



# From Split-Dist
treesim_SOURCES += label-map.hh label-map.cc
treesim_SOURCES += newick.h newick.c

# For Consensus tree
treesim_SOURCES += SNode.h SCNode.cc
treesim_SOURCES += SCTree.h SCTree.cc
treesim_SOURCES += buildtree.h buildtree.cc


# Main
treesim_SOURCES += treesim.cc






#eof

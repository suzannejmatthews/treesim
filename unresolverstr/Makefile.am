
#AM_CXXFLAGS = -static
AM_CXXFLAGS = -Wall
#CXXFLAGS += -O3 -pg
#CXXFLAGS += -g
AM_CXXFLAGS += -O3 -pg -fexpensive-optimizations

#AM_CXXFLAGS += -DNDEBUG

# misc.

# real program stuff
bin_PROGRAMS = unresolverstr

unresolverstr_SOURCES =



# From Split-Dist
unresolverstr_SOURCES += label-map.hh label-map.cc
unresolverstr_SOURCES += bitset.hh bitset.cc

# My classes
unresolverstr_SOURCES += hashfunc.hh hashfunc.cc
unresolverstr_SOURCES += hashmap.hh hashmap.cc

# From Mersenne Twist 
unresolverstr_SOURCES += mersenne.cc

# For Consensus tree
unresolverstr_SOURCES += SCNode.cc SCTree.cc


unresolverstr_SOURCES += newick.h newick.c

#unresolverstr_SOURCES += stocc.h stoc1.cpp


# Main
unresolverstr_SOURCES += unresolverstr.cc






#eof

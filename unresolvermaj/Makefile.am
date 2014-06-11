
#AM_CXXFLAGS = -static
AM_CXXFLAGS = -Wall
#CXXFLAGS += -O3 -pg
#CXXFLAGS += -g
AM_CXXFLAGS += -O3 -fexpensive-optimizations

#AM_CXXFLAGS += -DNDEBUG

# misc.

# real program stuff
bin_PROGRAMS = unresolvermaj

unresolvermaj_SOURCES =



# From Split-Dist
unresolvermaj_SOURCES += label-map.hh label-map.cc
unresolvermaj_SOURCES += bitset.hh bitset.cc

# My classes
unresolvermaj_SOURCES += hashfunc.hh hashfunc.cc
unresolvermaj_SOURCES += hashmap.hh hashmap.cc

# From Mersenne Twist 
unresolvermaj_SOURCES += mersenne.cc

# For Consensus tree
unresolvermaj_SOURCES += SCNode.cc SCTree.cc


unresolvermaj_SOURCES += newick.h newick.c


# Main
unresolvermaj_SOURCES += unresolvermaj.cc






#eof

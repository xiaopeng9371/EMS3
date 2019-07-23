#CC=gcc
CXX=g++
RM=rm -f
#CXXFLAGS=-Wall -O0 -g -fopenmp -pthread -std=c++11 -Wl,--no-as-needed
CXXFLAGS=-Wall -O3 -fopenmp -pthread -lpthread -std=c++11 
#-Wl,--no-as-needed
LDFLAGS=-fopenmp -pthread -std=c++11
#-Wl,--no-as-needed
LDLIBS=-lpthread

SRCS=ems.hpp \
     utils.h \
     motif.hpp \
     motif_set.hpp \
     motif_tree.hpp \
     motif_tree_base.hpp \
     motif_tree_fast.hpp \
     motif_trie.hpp \
     ems1.hpp \
     ems2.hpp \
     ems2p.hpp

.PHONY: all clean dist-clean

all: ems

ems.o: ems.cpp $(SRCS)
	$(CXX) $(CXXFLAGS) -c $<

ems: ems.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

clean:
	$(RM) ems.o

dist-clean: clean
	$(RM) ems 

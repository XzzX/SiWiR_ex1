SHELL        = /bin/bash
CXX          = g++
 
FLAGS        = -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG
DEBUGFLAGS   = 
RELEASEFLAGS =
INCPATH      = 
LIBPATH      = 
LIBS         = 

# blas
INCPATH     += -I/usr/lib64/atlas/include/
LIBPATH     += -L/usr/lib64/atlas/
LIBS        += -lcblas -latlas

# likwid
FLAGS       += -DUSE_LIKWID -pthread
INCPATH     += -I/usr/local/likwid-3.1.2/include/
LIBPATH     += -L/usr/local/likwid-3.1.2/lib/
LIBS        += -llikwid
 
TARGET       = foomatic-widget
COMMON       = Timer.h
 
all: matmult

%.o: %.cpp $(COMMON)
	$(CXX) -c $(FLAGS) $(INCPATH) $<
 
matmult:
	$(CXX) $(FLAGS) $(INCPATH) $(LIBPATH) $(LIBS) -o matmult matmult.cpp
	
compare:
	$(CXX) $(FLAGS) $(INCPATH) -o compare compare.cpp

clean:
	rm -f *.o matmult

.PHONY : all clean
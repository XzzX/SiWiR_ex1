SHELL        = /bin/bash
CXX          = g++
 
FLAGS        = -march=native -mtune=native -std=c++0x -Wall -Wextra -Wshadow -Werror -O3 -DNDEBUG
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
 
matmult: matmult.cpp
	$(CXX) $(FLAGS) $(INCPATH) -o matmult matmult.cpp $(LIBPATH) $(LIBS)
	
compare: compare.cpp
	$(CXX) $(FLAGS) $(INCPATH) -o compare compare.cpp
	
test: matmult compare
	rm -f C1.out
	rm -f C2.out
	rm -f C3.out
	likwid-pin -c 2 ./matmult matrices/testMatrices/A.in matrices/testMatrices/B.in C1.out
	likwid-pin -c 2 ./matmult matrices/testMatrices/A2.in matrices/testMatrices/B2.in C2.out
	likwid-pin -c 2 ./matmult matrices/testMatrices/A3.in matrices/testMatrices/B3.in C3.out
	./compare C1.out matrices/testMatrices/C.out
	./compare C2.out matrices/testMatrices/C2.out
	./compare C3.out matrices/testMatrices/C3.out
	
perf: matmult
	rm -f perf.txt
	./matmult matrices/perfMatrices/32x32-1 matrices/perfMatrices/32x32-2 32x32-3 >> perf.txt
	./matmult matrices/perfMatrices/64x64-1 matrices/perfMatrices/64x64-2 64x64-3 >> perf.txt
	./matmult matrices/perfMatrices/128x128-1 matrices/perfMatrices/128x128-2 128x128-3 >> perf.txt
	./matmult matrices/perfMatrices/256x256-1 matrices/perfMatrices/256x256-2 256x256-3 >> perf.txt
	./matmult matrices/perfMatrices/512x512-1 matrices/perfMatrices/512x512-2 512x512-3 >> perf.txt
	./matmult matrices/perfMatrices/1024x1024-1 matrices/perfMatrices/1024x1024-2 1024x1024-3 >> perf.txt
	./matmult matrices/perfMatrices/2048x2048-1 matrices/perfMatrices/2048x2048-2 2048x2048-3 >> perf.txt
	ipython plot.py
	
perfTest: matmult
	rm -f 32x32-3
	rm -f 64x64-3
	rm -f 128x128-3
	rm -f 256x256-3
	rm -f 512x512-3
	rm -f 1024x1024-3
	rm -f 2048x2048-3
	likwid-pin -c 2 ./matmult matrices/perfMatrices/32x32-1 matrices/perfMatrices/32x32-2 32x32-3
	likwid-pin -c 2 ./matmult matrices/perfMatrices/64x64-1 matrices/perfMatrices/64x64-2 64x64-3
	likwid-pin -c 2 ./matmult matrices/perfMatrices/128x128-1 matrices/perfMatrices/128x128-2 128x128-3
	likwid-pin -c 2 ./matmult matrices/perfMatrices/256x256-1 matrices/perfMatrices/256x256-2 256x256-3
	likwid-pin -c 2 ./matmult matrices/perfMatrices/512x512-1 matrices/perfMatrices/512x512-2 512x512-3
	likwid-pin -c 2 ./matmult matrices/perfMatrices/1024x1024-1 matrices/perfMatrices/1024x1024-2 1024x1024-3
	likwid-pin -c 2 ./matmult matrices/perfMatrices/2048x2048-1 matrices/perfMatrices/2048x2048-2 2048x2048-3
	./compare 32x32-3 matrices/perfMatrices/32x32-3
	./compare 64x64-3 matrices/perfMatrices/64x64-3
	./compare 128x128-3 matrices/perfMatrices/128x128-3
	./compare 256x256-3 matrices/perfMatrices/256x256-3
	./compare 512x512-3 matrices/perfMatrices/512x512-3
	./compare 1024x1024-3 matrices/perfMatrices/1024x1024-3
	./compare 2048x2048-3 matrices/perfMatrices/2048x2048-3

clean:
	rm -f *.o matmult

.PHONY : all clean test perf

SHELL := /bin/bash

C = -O3 $(OPTIONS) -DOPENMP -fopenmp -Wall 
CC = g++-mp-6 

SOURCEDIR = ../src

I = -I/Users/hans/local/include
L = -L/Users/hans/local/lib -lfftw3

TARGETS := PoissonSolver
all: $(TARGETS)

# OBJECT FILES
POISSON_OBJS = poisson.o $(SOURCEDIR)/grid.o $(SOURCEDIR)/multigrid.o $(SOURCEDIR)/multigrid_solver.o

# HEADERS
HEADERS = poisson_solver.h $(SOURCEDIR)/grid.h $(SOURCEDIR)/multigrid.h $(SOURCEDIR)/multigrid_solver.h

PoissonSolver: $(POISSON_OBJS)
	$(CC) -std=c++11 -o $@ $^ $C $I $L

%.o: %.cpp $(HEADERS)
	$(CC) -std=c++11  -c -o $@ $< $C $I $L

clean:
	rm -rf $(TARGETS) *.o $(SOURCEDIR)/*.o


#ifndef _MULTIGRID_HEADER
#define _MULTIGRID_HEADER
#include <assert.h>
#include <iostream>
#include <vector>
#include <climits>
#include <complex>
#include "grid.h"

    //=========================================
    //                                       // 
    // A stack of _Nlevel grids with         //
    // N^NDIM / 2^Level cells in each level  //
    //                                       //
    // Bounds-check for array lookups:       //
    // #define _BOUNDSCHECK                  //
    //                                       //
    //=========================================

template<unsigned int NDIM, typename T>
class MultiGrid {
  private:

    unsigned int _N;                        // Number of cells per dim in domain-grid [0]
    unsigned int _Ntot;                     // Total number of cells in domain-grid [0]
    unsigned int _Nlevel;                   // Number of levels
    std::vector<unsigned int> _NinLevel;    // Number of cells per dim in each level
    std::vector<unsigned int> _NtotinLevel; // Total number of cells in each level
    std::vector<Grid<NDIM,T> > _y;          // The grid data

  public:

    // Constructors
    MultiGrid() {}
    MultiGrid(unsigned int N):  MultiGrid(N, int(log2(N)+1)) {}
    MultiGrid(unsigned int N, unsigned int Nlevel);
    MultiGrid(Grid<NDIM, T> &y, unsigned int Nlevel);
    MultiGrid(Grid<NDIM, T> &y);
    
    // Fetch a reference to the solution grid at a given level
    Grid<NDIM,T>& get_grid(unsigned int level = 0);

    // Fetch a pointer to the underlying array at each level
    T* operator[](unsigned int level);
    T* get_y(unsigned int level);

    // Fetch the value in the grid at a given level and index
    T get_y(unsigned int level, unsigned int i);
    
    // Fetch the value in the grid at a given level and coordinates (ix,iy...)
    T get_y(unsigned int level, std::vector<unsigned int>& coord_list);

    // Fetch info about the grid
    unsigned int get_N(unsigned int level = 0);
    unsigned int get_Ntot(unsigned int level = 0);
    unsigned int get_Ndim();
    unsigned int get_Nlevel();
    unsigned int get_Nmin();
  
    // Set the value of y at given level and index (save way to define value)
    void set_y(unsigned int level, unsigned int i, T value);
    
    // Gridindex from coordinate and vice versa
    unsigned int gridindex_from_coord(unsigned int level, std::vector<unsigned int>& coord_list);
    std::vector<unsigned int> coord_from_gridindex(unsigned int level, unsigned int i);

    // Free up all memory and reset all variables
    void clear();

    // Restrict down a grid 
    void restrict_down(unsigned int from_level); 
    void restrict_down(unsigned int from_level, Grid<NDIM, T> &to_grid); 
    void restrict_down_all();
};

#endif

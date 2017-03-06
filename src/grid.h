#ifndef _GRID_HEADER
#define _GRID_HEADER
#include <assert.h>
#include <cstring>  
#include <iostream> 
#include <fstream>
#include <complex>
#include <vector>
#include <climits>

    //=========================================
    //                                       // 
    // A simple multidimensional grid-class  //
    //                                       //
    // Bounds-check for array lookups        //
    // #define _BOUNDSCHECK                  //
    //                                       //
    //=========================================

template<unsigned int NDIM, typename T>
class Grid {

  private:

    unsigned int _N;      // Number of cells per dim in the grid
    unsigned int _Ntot;   // Total number of cells in the grid
    std::vector<T> _y;    // The grid data

  public:

    // Constructors
    Grid() : Grid(0, 0.0) {}
    Grid(unsigned int N) : Grid(N, 0.0) {}
    Grid(unsigned int N, T yini);

    // Get a pointer to the T-array
    T* get_y();

    // Allow syntax grid[i] to get/set the index = i'th element
    T& operator[](unsigned int i);

    // Fetch the index = i element in the grid
    T get_y(unsigned int i);

    // Assign value in the grid
    void set_y(std::vector<T> &y);
    void set_y(unsigned int i, T &value);

    // Grid-index -> coordinate list [ i = ix1 + N * ix2 + N^2 * ix3 + ... ]
    std::vector<unsigned int> index_list(unsigned int i);

    // Get some info about the grid
    unsigned int get_N();
    unsigned int get_Ntot();

    // Convert coordiates -> index in the grid
    unsigned int grid_index(std::vector<unsigned int> &index_list);
    unsigned int grid_index_3d(unsigned int ix, unsigned int iy, unsigned int iz);
    unsigned int grid_index_2d(unsigned int ix, unsigned int iy);

    // Dump a grid to file
    void dump_to_file(std::string filename);

    // Read a grid from file into the object
    void read_from_file(std::string filename);

    // Maximum (in norm)
    double max(){ 
      double maxval = std::norm(_y[0]);
#ifdef OPENMP
#pragma omp parallel for reduction(max: maxval)
#endif
      for(unsigned int i = 0; i < _Ntot; i++){
        double curval = std::norm(_y[i]);
        if(curval > maxval) maxval = curval;
      }
      return std::sqrt( maxval );
    }
    
    // Maximum (in norm)
    double min(){ 
      double minval = std::norm(_y[0]);
#ifdef OPENMP
#pragma omp parallel for reduction(min: minval)
#endif
      for(unsigned int i = 0; i < _Ntot; i++){
        double curval = std::norm(_y[i]);
        if(curval < minval) minval = curval;
      }
      return std::sqrt( minval );
    }

    // Free up all memory and reset all variables
    void clear();

    // Operator overloading: add two grids element by element
    template<unsigned int NNDIM, typename TT>
    Grid<NNDIM,TT>& operator+=(const Grid<NNDIM,TT>& rhs){
#ifdef _BOUNDSCHECK
      assert(this->_N == rhs._N);
#endif
#ifdef OPENMP
#pragma omp parallel for
#endif
      for(unsigned int i = 0; i < _Ntot; i++)
        this->_y[i] += rhs._y[i];
      return *this;      
    }
    
    // Operator overloading: subtract two grids element by element
    template<unsigned int NNDIM, typename TT>
    Grid<NNDIM,TT>& operator-=(const Grid<NNDIM,TT>& rhs){
#ifdef _BOUNDSCHECK
      assert(this->_N == rhs._N);
#endif
#ifdef OPENMP
#pragma omp parallel for
#endif
      for(unsigned int i = 0; i < _Ntot; i++)
        this->_y[i] -= rhs._y[i];
      return *this;      
    }

    // Operator overloading: multiply two grids element by element
    template<unsigned int NNDIM, typename TT>
    Grid<NNDIM,TT>& operator*=(const Grid<NNDIM,TT>& rhs){
#ifdef _BOUNDSCHECK
      assert(this->_N == rhs._N);
#endif
#ifdef OPENMP
#pragma omp parallel for
#endif
      for(unsigned int i = 0; i < _Ntot; i++)
        this->_y[i] *= rhs._y[i];
      return *this;      
    }
    
    // Operator overloading: multiply two grids element by element
    template<unsigned int NNDIM, typename TT>
    Grid<NNDIM,TT>& operator/=(const Grid<NNDIM,TT>& rhs){
#ifdef _BOUNDSCHECK
      assert(this->_N == rhs._N);
#endif
#ifdef OPENMP
#pragma omp parallel for
#endif
      for(unsigned int i = 0; i < _Ntot; i++)
        this->_y[i] /= rhs._y[i];
      return *this;      
    }
    
    // Operator overloading: multiply every element in grid by scalar
    Grid<NDIM,T>& operator *=(const T & rhs){ 
#ifdef OPENMP
#pragma omp parallel for
#endif
      for(unsigned int i = 0; i < _Ntot; i++) 
        this->_y[i] *= rhs; 
      return *this;
    } 
    
    // Operator overloading: divide every element in grid by scalar
    Grid<NDIM,T>& operator /=(const T & rhs){
#ifdef OPENMP
#pragma omp parallel for
#endif
      for(unsigned int i = 0; i < _Ntot; i++) 
        this->_y[i] /= rhs; 
      return *this;
    }

    // The rms-norm, sqrt[ Sum y[i]^2 / Ntot ], of the grid
    double rms_norm();

    // Check for NaN and exit if true
    void check_for_nan(bool exitifnan);
};
 
template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator+(Grid<NDIM,T> lhs, const Grid<NDIM,T>& rhs){
  lhs += rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator-(Grid<NDIM,T> lhs, const Grid<NDIM,T>& rhs){
  lhs -= rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator*(Grid<NDIM,T> lhs, const Grid<NDIM,T>& rhs){
  lhs *= rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator/(Grid<NDIM,T> lhs, const Grid<NDIM,T>& rhs){
  lhs /= rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator*(Grid<NDIM,T> lhs, const T& rhs){
  lhs *= rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator/(Grid<NDIM,T> lhs, const T& rhs){
  lhs /= rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator+(Grid<NDIM,T> lhs, const T& rhs){
  lhs += rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> operator-(Grid<NDIM,T> lhs, const T& rhs){
  lhs -= rhs;
  return lhs;
}

template<unsigned int NDIM, typename T>
Grid<NDIM,T> sqrt(Grid<NDIM,T> lhs){
  for(unsigned int i = 0; i < lhs.get_Ntot(); i++)
    lhs[i] = sqrt(fabs(lhs[i]));
  return lhs;
}

template<unsigned int NDIM, typename T>
double Grid<NDIM,T>::rms_norm(){
  double rms = 0.0;
#ifdef OPENMP
#pragma omp parallel for reduction(+:rms)
#endif
  for(unsigned int i = 0; i < _Ntot; i++){
    rms += std::norm(_y[i]);
  }
  rms = std::sqrt(rms / double(_Ntot));
  return rms;
}

#endif

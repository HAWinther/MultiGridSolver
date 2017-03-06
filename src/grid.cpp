#include "grid.h"

// Simple int-int a^b power-function
inline unsigned int power(unsigned int a, unsigned int b){
  unsigned int res = 1;
  for(unsigned int i = 0; i < b; i++) {
    res *= a;
  }
#ifdef _BOUNDSCHECK
  assert( (pow(1.0*a,1.0*b) - double(res)) < 0.5);
#endif
  return res;
}
    
template<unsigned int NDIM, typename T>
void Grid<NDIM,T>::check_for_nan(bool exitifnan){
  bool nanfound = false;
  for(unsigned int i = 0; i < _Ntot; i++){
    if(_y[i] != _y[i]){
      nanfound = true;
      break;
    }
  }

  if(nanfound){
    std::cout << "Warning: NaN found in grid" << (exitifnan ? " ...aborting!" : "" ) << std::endl;
    if(exitifnan) exit(1);
  }
}

// Constructor with intial value
template<unsigned int NDIM, typename T>
Grid<NDIM,T>::Grid(unsigned int N, T yini) : _N(N), _Ntot(power(_N, NDIM)), _y(std::vector<T>(_Ntot, yini)) {}

// Fetch pointer to grid
template<unsigned int NDIM, typename T>
T* Grid<NDIM,T>::get_y() { 
  return &_y[0]; 
}

// Allow to fetch value using f[i] syntax
template<unsigned int NDIM, typename T>
T& Grid<NDIM,T>::operator[](unsigned int i){ 
#ifdef _BOUNDSCHECK
  assert(i < _Ntot);
#endif
  return _y[i]; 
}

// Fetch value of grid-cell [i]
template<unsigned int NDIM, typename T>
T Grid<NDIM,T>::get_y(unsigned int i) { 
#ifdef _BOUNDSCHECK
  assert(i < _Ntot);
#endif
  return _y[i]; 
}

// Assign whole grid from vector
template<unsigned int NDIM, typename T>
void Grid<NDIM,T>::set_y(std::vector<T> &y){
#ifdef _BOUNDSCHECK
  assert(y.size() == _Ntot);
#endif
  _y = y;
}

// Assign the gridcell [i] with [value]
template<unsigned int NDIM, typename T>
void Grid<NDIM,T>::set_y(unsigned int i, T &value){
#ifdef _BOUNDSCHECK
  assert(i < _Ntot);
#endif
  _y[i] = value;
}

// Compute coordinates given a gridindex
template<unsigned int NDIM, typename T>
std::vector<unsigned int> Grid<NDIM,T>::index_list(unsigned int i){
  std::vector<unsigned int> ii(NDIM, 0);
  for(unsigned int j = 0, n = 1; j < NDIM; j++, n *= _N){
    ii[j] = i / n % _N;
  }
  return ii;
}

// Coordinates -> grid-index (index in the 1D _y vector)
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM,T>::grid_index(std::vector<unsigned int> &index_list){
  unsigned int index = 0;
  for(unsigned int j = 0, n = 1; j < NDIM; j++, n *= _N)
    index += index_list[j] * n;
#ifdef _BOUNDSCHECK
  assert(index < _Ntot);
#endif
  return index;
}
    
// Coordinate -> gridindex for 3D grid
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM,T>::grid_index_3d(unsigned int ix, unsigned int iy, unsigned int iz){
  return ix + _N*(iy + _N*iz);
}

// Coordinate -> gridindex for 2D grid
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM,T>::grid_index_2d(unsigned int ix, unsigned int iy){
  return ix + _N*iy;
}
// Returns number of cells per dim
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM,T>::get_N(){
  return _N;
}

// Return total number of cells
template<unsigned int NDIM, typename T>
unsigned int Grid<NDIM,T>::get_Ntot(){
  return _Ntot;
}

// Write a grid to file
template<unsigned int NDIM, typename T>
void Grid<NDIM,T>::dump_to_file(std::string filename){
  unsigned int ndim = NDIM;
  
  // Verbose
  std::cout << "==> Dumping grid to file [" << filename << "]" << std::endl;
  std::cout << "    Ndim: " << NDIM << " N: " << _N << " Ntot: " << _Ntot << std::endl;

  // Write header
  std::ofstream fout(filename.c_str(), std::ios::out | std::ios::binary);
  fout.write((char*)&_N, sizeof(unsigned int));
  fout.write((char*)&ndim, sizeof(unsigned int));

  // Write the grid-data
  fout.write((char*)&_y[0], _y.size() * sizeof(T));
}

// Read a grid from file (assumes specific format)
template<unsigned int NDIM, typename T>
void Grid<NDIM,T>::read_from_file(std::string filename){
  unsigned int ninfile = 0, ntot = 0, ndim = 0, size = 0;
  
  // Read header: N and NDIM
  std::ifstream input(filename.c_str(), std::ios::in | std::ifstream::binary);
  
  input.read((char *) &ninfile, sizeof(unsigned int));
  input.read((char *) &ndim, sizeof(unsigned int));

  ntot = power(ninfile, NDIM);
  size = sizeof(T) * ntot;
  
  // Checks
  assert(ndim == NDIM);
  assert(ninfile > 0 && ntot < INT_MAX);

  // Verbose
  std::cout << "==> Reading file into grid [" << filename << "]" << std::endl;
  std::cout << "    Ndim: " << ndim << " Nfile: " << ninfile << " Ntot: " << ntot << std::endl;

  // Read the data
  std::vector<char> tempvec(size, 0);
  input.read(&tempvec[0], tempvec.size());
 
  // Copy the grid-data and set parameters
  _N = ninfile;
  _Ntot = ntot;
  _y = std::vector<T>(ntot, 0.0);
  std::memcpy(&_y[0], &tempvec[0], size);
}

template<unsigned int NDIM, typename T>
void Grid<NDIM,T>::clear(){
  _N = _Ntot = 0;
  _y.clear();
}

// Explicit template specialisation
template class Grid<3,double>;
template class Grid<2,double>;
template class Grid<1,double>;
template class Grid<3,float>;
template class Grid<2,float>;
template class Grid<1,float>;
template class Grid<1,std::complex<double> >;


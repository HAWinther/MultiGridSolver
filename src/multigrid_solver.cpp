#include "multigrid_solver.h"

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

//================================================
// The L operator ( EOM is written as L(f) = 0 )
//================================================

template<unsigned int NDIM, typename T>
T MultiGridSolver<NDIM,T>::l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource){
  T h, l, source, kinetic;
  unsigned int i = index_list[0];

  // Gridspacing h
  h = 1.0/T( get_N(level) );

  // Compute the standard kinetic term [D^2 f]
  kinetic = -2.0 * NDIM * _f[level][ i ];
  for(unsigned int k = 1; k < 2*NDIM + 1; k++){
    kinetic += _f[level][ index_list[k] ];
  }

  // The right hand side of the PDE
  T *rho = _ext_field[0]->get_y(level);
  source = rho[i];

  // Add the source term arising from restricting the equation down to the lower level
  if( level > 0 && addsource ){
    source += _source[level][i];
  }

  // The full equation
  l = kinetic/(h*h) - source;

  return l;
}

//================================================
// The derivative of the L-operator dL/df
// or more accurately dL_{ijk..} / df_{ijk..}
//================================================

template<unsigned int NDIM, typename T>
T MultiGridSolver<NDIM,T>::dl_operator(unsigned int level, std::vector<unsigned int>& index_list){
  T dl, h;

  // Gridspacing
  h = 1.0/T( get_N(level) );

  // The derivtive dL/df
  dl = -2.0*NDIM/(h*h);

  // Sanity check
  if(fabs(dl) < 1e-10){
    std::cout << "Error: dl close to 0" << std::endl;
    exit(1);
  }

  return dl;
}

//================================================
// The driver routine for solving the PDE
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::solve(){
  // Init some variables
  _istep_vcycle = 0;
  _tot_sweeps_domain_grid = 0;
  _res_domain_array.clear();

  if(_verbose){
    std::cout << std::endl;
    std::cout << "===============================================================" << std::endl;
    std::cout << "==> Starting multigrid solver" << std::endl;
    std::cout << "===============================================================\n" << std::endl;
  }

  // Pre-solve on domaingrid
  solve_current_level(0);

  // Set the initial residual 
  _rms_res_i = _rms_res;
  _rms_res_old = 0.0;

  // Check if we already have convergence
  if(check_convergence()) return;

  // The V-cycle
  while(1) {
    ++_istep_vcycle;
    
    if(_verbose){
      std::cout << std::endl;
      std::cout << "===============================================================" << std::endl;
      std::cout << "==> Starting V-cycle istep = " << _istep_vcycle << " Res = " << _rms_res << std::endl;
      std::cout << "===============================================================\n" << std::endl;
    }

    // Go down to the bottom (from finest grid [0] to coarsest grid [_Nlevel-1])
    recursive_go_down(0);

    // Go up to the top
    recursive_go_up(_Nlevel-2);
    
    // Check for errors in the computation (NaN) and exit if true
    _f.get_grid(0).check_for_nan(true);
    
    // Check for convergence
    if(check_convergence()) break;
  }
}

template<unsigned int NDIM, typename T>
T* MultiGridSolver<NDIM,T>::get_y(unsigned int level){ 
  return _f.get_y(level); 
}

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::set_epsilon(double eps_converge){ 
  _eps_converge = eps_converge; 
}

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::add_external_grid(MultiGrid<NDIM,T> *field){ 
  assert(field->get_Nlevel() >= _Nlevel);
  _ext_field.push_back(field); 
}

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::set_maxsteps(unsigned int maxsteps){ 
  _maxsteps = maxsteps; 
}

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::set_convergence_criterion_residual(bool use_residual){
  _conv_criterion_residual = use_residual;
}

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::set_ngs_sweeps(unsigned int ngs_fine, unsigned int ngs_coarse){
  _ngs_fine = ngs_fine; 
  _ngs_coarse = ngs_coarse; 
}

template<unsigned int NDIM, typename T>
unsigned int  MultiGridSolver<NDIM,T>::get_N(unsigned int level){ 
  return _f.get_N(level); 
}

template<unsigned int NDIM, typename T>
unsigned int  MultiGridSolver<NDIM,T>::get_Ntot(unsigned int level){ 
  return _f.get_Ntot(level); 
}

template<unsigned int NDIM, typename T>
MultiGridSolver<NDIM,T>::MultiGridSolver(unsigned int N, unsigned int Nmin, bool verbose) :
  _N(N), _Ntot(power(N, NDIM)), _Nmin(Nmin), _Nlevel(int(log2(N) - _Nmin + 2)), _verbose(verbose), 
  _rms_res(0.0), _rms_res_i(0.0), _rms_res_old(0.0) {

    // Check that N is divisible by 2^{Nlevel - 1} which is required for the restriction to make sense
    assert( ( _N / power(2, _Nlevel - 1) ) * power(2, _Nlevel - 1) == _N); 

    // We should have atleast 1 level
    assert(_Nlevel >= 1);                 

    // Allocate memory
    _f      = MultiGrid<NDIM,T>(_N, _Nlevel);
    _source = MultiGrid<NDIM,T>(_N, _Nlevel);
    _res    = MultiGrid<NDIM,T>(_N, _Nlevel);
  }

//================================================
// The initial guess for the solver at the 
// domain level (level = 0)
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::set_initial_guess(T guess){
  T *f = _f.get_y(0);
  std::fill_n(f, _Ntot, guess);
}

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::set_initial_guess(T *guess){
  T *f = _f.get_y(0);
  std::copy( &guess[0], &guess[0] + _Ntot, &f[0] );
}

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::set_initial_guess(Grid<NDIM,T>& guessgrid){
  T *f = _f.get_y(0);
  T *guess = guessgrid.get_y();
  std::copy( &guess[0], &guess[0] + _Ntot, &f[0] );
}

//================================================
// Given a cell i = (ix,iy,iz, ...) it computes
// the grid-index of the 2NDIM neighboring cells
// 0: (ix  ,iy  , iz, ...)
// 1: (ix-1,iy  , iz, ...)
// 2: (ix+1,iy  , iz, ...)
// 3: (ix,  iy-1, iz, ...)
// 4: (ix,  iy+1, iz, ...)
// ...
//
// Assuming periodic boundary conditions
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::get_neighbor_gridindex(std::vector<unsigned int>& index_list, unsigned int i, unsigned int N){
  index_list = std::vector<unsigned int>(2*NDIM+1);
  index_list[0] = i;
  for(unsigned int j = 0, n = 1; j < NDIM; j++, n *= N){
    unsigned int ii = i/n % N;
    unsigned int iminus = ii >= 1   ? ii - 1 : N - 1;
    unsigned int iplus  = ii <= N-2 ? ii + 1 : 0;
    index_list[2*j+1] = i + (iminus - ii) * n;
    index_list[2*j+2] = i + (iplus - ii) * n;
  }
}

//================================================
// Calculates the residual in each cell at
// a given level and stores it in [res]. Returns 
// the rms-residual over the whole grid
//================================================

template<unsigned int NDIM, typename T>
double MultiGridSolver<NDIM,T>::calculate_residual(unsigned int level, Grid<NDIM,T> &res){
  unsigned int N    = get_N(level);
  unsigned int Ntot = get_Ntot(level);

  // Calculate and store (minus) the residual in each cell
#ifdef OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < Ntot; i++) {
    std::vector<unsigned int> index_list;
    get_neighbor_gridindex(index_list, i, N);
    res[i] = -l_operator(level, index_list, true);
  }

  // Calculate and return RMS residual
  return res.rms_norm();
}

//================================================
// Criterion for defining convergence.
// Standard ways are: based on residual or
// the ratio of the residual to the initial
// residual (err).
//================================================

template<unsigned int NDIM, typename T>
bool MultiGridSolver<NDIM,T>::check_convergence(){
  // Compute ratio of residual to initial residual
  double err = _rms_res_i != 0.0 ? _rms_res/_rms_res_i : 1.0;
  bool converged = false;

  // Print out some information
  if(_verbose){
    std::cout << "    Checking for convergence at step = " << _istep_vcycle << std::endl;
    std::cout << "        Residual = " << _rms_res << "  Residual_old = " <<  _rms_res_old << std::endl;
    std::cout << "        Residual_i = " << _rms_res_i << "  Err = " << err << std::endl;
  }

  // Convergence criterion
  if(_conv_criterion_residual) {

    // Convergence criterion based on the residual
    if(_rms_res < _eps_converge){
      if(_verbose || true){
        std::cout << std::endl;
        std::cout << "    The solution has converged res = " << _rms_res << " < " << _eps_converge << " istep = " << _istep_vcycle << "\n" << std::endl;
      }
      converged = true;
    } else {
      if(_verbose){
        std::cout << "    The solution has not yet converged res = " << _rms_res << " !< " << _eps_converge << std::endl; 
      }
    }

  } else {

    // Convergence criterion based on the ratio of the residual
    if(err < _eps_converge){
      if(_verbose || true){
        std::cout << std::endl;
        std::cout << "    The solution has converged err = " << err << " < " << _eps_converge << " ( res = " << _rms_res << " ) istep = " << _istep_vcycle << "\n" << std::endl;
      }
      converged = true;
    } else {
      if(_verbose){
       std::cout << "    The solution has not yet converged err = " << err << " !< " << _eps_converge << std::endl;
      }
    }
  }

  if(_verbose && (_rms_res > _rms_res_old && _istep_vcycle > 1) ){
    std::cout << "    Warning: Residual_old > Residual" << std::endl;
  }

  // Define converged if istep exceeds maxsteps to avoid infinite loop...
  if(_istep_vcycle >= _maxsteps){
    std::cout << "    WARNING: MultigridSolver failed to converge! Reached istep = maxsteps = " << _maxsteps << std::endl;
    std::cout << "    res = " << _rms_res << " res_old = " << _rms_res_old << " res_i = " << _rms_res_i << std::endl;
    converged  = true;
  }

  return converged;
}

//================================================
// Prolonge up solution phi from course grid
// to fine grid. Using trilinear prolongation
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::prolonge_up_array(unsigned int to_level, Grid<NDIM,T>& Bottom, Grid<NDIM,T>& Top){
  unsigned int twotondim = 1 << NDIM;
  unsigned int NTop      = get_N(to_level);
  unsigned int NtotTop   = get_Ntot(to_level);
  unsigned int NBottom   = NTop/2;
  
  // Compute NTop, Ntop^2, ... , Ntop^{Ndim-1} and similar for Nbottom
  std::vector<unsigned int> nBottomPow(NDIM, 1);
  std::vector<unsigned int> nTopPow(NDIM, 1);
  for(unsigned int j = 0, n = 1, m = 1; j < NDIM; j++, n *= NBottom, m *= NTop){
    nBottomPow[j] = n;
    nTopPow[j] = m;
  }

  // Trilinear prolongation
#ifdef OPENMP
#pragma omp parallel for  
#endif
  for (unsigned int i = 0; i < NtotTop; i++) {
    std::vector<double> fac(NDIM, 0.0);
    std::vector<int> iplus(NDIM, 0);

    double norm = 1.0;
    int iBottom = 0;

    // Compute the shift in index from iBottom to the cells corresponding to ix -> ix+1
    // The fac is the weight for the trilinear interpolation
    for(unsigned int j = 0; j < NDIM; j++){
      unsigned int ii = i/nTopPow[j] % NTop;
      unsigned int iiBottom = ii/2;
      iplus[j] = (iiBottom == NBottom-1 ? 1 - NBottom : 1) * nBottomPow[j];
      fac[j]   = ii % 2 == 0 ? 0.0 : 1.0;
      iBottom += iiBottom * nBottomPow[j];
      norm    *= (1.0 + fac[j]);
    }

    // Compute the sum Top[i] = Sum fac_i             * Top[iBottom + d_i] 
    //                        + Sum fac_i fac_j       * Top[iBottom + d_i + d_j] 
    //                        + Sum fac_i fac_j fac_k * Top[iBottom + d_i + d_j + d_k]
    //                        + ... +
    //                        + fac_1 ... fac_NDIM * Top[iBottom + d_1 + ... + d_NDIM]
    Top[i] = Bottom[iBottom];
    for(unsigned int k = 1; k < twotondim; k++){
      double termfac = 1.0;
      int iAdd = 0;
      std::bitset<NDIM> bits = std::bitset<NDIM>(k);
      for(unsigned int j = 0; j < NDIM; j++){
        iAdd = bits[j] * iplus[j];
        termfac *= 1.0 + bits[j] * (fac[j] - 1.0) ;
      }
      Top[i] += T(termfac) * Bottom[iBottom + iAdd];
    }
    Top[i] *= 1.0/T(norm);
  }
}

//================================================
// The Gauss-Seidel Sweeps with standard chess-
// board (first black then white) ordering of 
// gridnodes if _ngridcolours = 2 
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::GaussSeidelSweep(unsigned int level, unsigned int curcolor, T *f){
  unsigned int N    = get_N(level);
  unsigned int Ntot = get_Ntot(level);

#ifdef OPENMP
#pragma omp parallel for 
#endif
  for (unsigned int i = 0; i < Ntot; i++) {
    // Compute cell-color
    unsigned int color = 0;
    for(unsigned int j = 0, n = 1; j < NDIM; j++, n *= N){
      color += ( i / n % N );
    }
    color = color % _ngridcolours;

    // Only select cells with right color
    if( color == curcolor ){

      // Update the solution f = f - L / (dL/df)
      std::vector<unsigned int> index_list;
      get_neighbor_gridindex(index_list, i, N);
      T l  =  l_operator(level, index_list, true);
      T dl = dl_operator(level, index_list);
      f[i] -= l/dl;
    }
  }
}

//================================================
// Solve the equation on the current level
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::solve_current_level(unsigned int level){
  unsigned int ngs_sweeps;
 
  if(_verbose)
    std::cout << "    Performing Newton-Gauss-Seidel sweeps at level " << level << std::endl;

  // Number of sweeps we do
  if(level == 0)
    ngs_sweeps = _ngs_fine;
  else
    ngs_sweeps = _ngs_coarse;

  // Do N Gauss-Seidel Sweeps
  for (unsigned int i = 0; i < ngs_sweeps; i++) {
    if(level == 0)
      ++_tot_sweeps_domain_grid;

    // Sweep through grid according to sum of coord's mod _ngridcolours
    // Standard is _ngridcolours = 2 -> chess-board ordering
    for(unsigned int j = 0; j < _ngridcolours; j++)
      GaussSeidelSweep(level, j, _f[level]);

    // Calculate residual and output quite often.
    // For debug, but this is quite useful so keep it for now
    if(_verbose){
      if( (level > 0 && (i == 1 || i == ngs_sweeps-1) ) || (level == 0) ){
        std::cout << "        level = " << std::setw(5) << level << " NGS Sweep = " << std::setw(5) << i;
        std::cout << " Residual = " << std::setw(10) << calculate_residual(level, _res.get_grid(level)) << std::endl;
      }
    }

    // Compute and store the residual after every sweep on the domaingrid
    if(level == 0 && _store_all_residual){
      _res_domain_array.push_back( calculate_residual(level, _res.get_grid(level)) );
    }
  }
  if(_verbose) std::cout << std::endl;

  // Compute the residual
  double curres = calculate_residual(level, _res.get_grid(level));

  // Store domaingrid residual
  if (level == 0){
    _rms_res_old = _rms_res;
    _rms_res = curres;
  }

}

//================================================
// V-cycle go all the way up
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::recursive_go_up(unsigned int to_level){
  unsigned int from_level = to_level + 1;

  // Restrict down R[f] and store in _res (used as temp-array)
  _f.restrict_down(to_level, _res.get_grid(from_level));

  // Make prolongation array ready at from_level
  make_prolongation_array(_f.get_grid(from_level), _res.get_grid(from_level), _res.get_grid(from_level));

  // Prolonge up solution from-level to to-level and store in _res (used as temp array)
  if(_verbose)
    std::cout << "    Prolonge solution from level: " << to_level+1 << " -> " << to_level << std::endl;
  prolonge_up_array(to_level, _res.get_grid(from_level), _res.get_grid(to_level));

  // Correct solution at to_level (temp array _res contains the correction P[f-R[f]])
  _f.get_grid(to_level) += _res.get_grid(to_level);

  // Calculate new residual
  calculate_residual(to_level, _res.get_grid(to_level));

  // Solve on the level we just went up to
  solve_current_level(to_level);

  // Continue going up
  if(to_level > 0)
    recursive_go_up(to_level-1);
  else {
    return;
  }
}

//================================================
// Make the array we are going to prolonge up
// Assumes [Rf] contains the restiction of f
// from the upper level and returns [df]
// containing df = f - R[f]
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::make_prolongation_array(Grid<NDIM,T>& f, Grid<NDIM,T>& Rf, Grid<NDIM,T>& df){
  unsigned int Ntot = f.get_Ntot();
  df = f - Rf;
  return;
#ifdef OPENMP
#pragma omp parallel for
#endif
  for (unsigned int i = 0; i < Ntot; i++){
    df[i] = f[i] - Rf[i];
  }
}

//================================================
// Make new source
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::make_new_source(unsigned int level){
  unsigned int N    = get_N(level);
  unsigned int Ntot = get_Ntot(level);

  // Calculate the new source
#ifdef OPENMP
#pragma omp parallel for
#endif
  for(unsigned int i = 0; i < Ntot; i++){
    std::vector<unsigned int> index_list;
    get_neighbor_gridindex(index_list, i, N);
    T res = l_operator(level, index_list, false);
    _source[level][i] = _res[level][i] + res;
  }
}

//================================================
// V-cycle go all the way down
//================================================

template<unsigned int NDIM, typename T>
void MultiGridSolver<NDIM,T>::recursive_go_down(unsigned int from_level){
  unsigned int to_level = from_level + 1;

  // Check if we are at the bottom
  if(to_level >= _Nlevel) {
    if(_verbose) {
      std::cout << "    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << std::endl;
      std::cout << "    We have reached the bottom level = " << from_level << " Start going up." << std::endl;
      std::cout << "    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n" << std::endl;
    }
    return;
  }
  
  if(_verbose)
    std::cout << "    Going down from level " << from_level << " -> " << to_level << std::endl;

  // Restrict residual and solution
  _res.restrict_down(from_level, _res.get_grid(from_level + 1));
  _f.restrict_down(from_level, _f.get_grid(from_level + 1));

  // Make new source
  make_new_source(to_level);

  // Solve on current level
  solve_current_level(to_level);

  // Recursive call
  recursive_go_down(to_level);
}

template<unsigned int NDIM, typename T>
T MultiGridSolver<NDIM,T>::get_multigrid_source(unsigned int level, unsigned int i){ 
  return _source[level][i]; 
}

template<unsigned int NDIM, typename T>
T* MultiGridSolver<NDIM,T>::get_external_field(unsigned int level, unsigned int field){ 
  return _ext_field[field]->get_y(level);
}

template<unsigned int NDIM, typename T>
void  MultiGridSolver<NDIM,T>::clear() {
  _N = _Ntot = _Nlevel = 0;
  _f.clear();
  _res.clear();
  _source.clear();
  _ext_field.clear();
  _res_domain_array.clear();
}

// Explicit template specialisation
template class MultiGridSolver<3,double>;
template class MultiGridSolver<2,double>;
template class MultiGridSolver<1,double>;
template class MultiGridSolver<3,float>;
template class MultiGridSolver<2,float>;
template class MultiGridSolver<1,float>;
template class MultiGridSolver<1,std::complex<double> >;

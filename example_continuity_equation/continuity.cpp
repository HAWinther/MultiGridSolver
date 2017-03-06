#include "continuity_solver.h"
#include <fstream>
#include <fftw3.h>

//====================================
//
// This is an example of how to use the
// multigrid solver class for solving
// the continuity equation (assuming 
// ddelta/dt = f delta)
//                                       
//    fdelta + D[(1+delta) v] = 0  
// 
//====================================

//====================================
// The dim of the equation
//====================================
#define _NDIM 3

//====================================
// The type (double,float,complex,...)
//====================================
typedef double DType;
 
//====================================
// Parameters for the PDE
//====================================
DType rhomin_in_grad = 1.0;
DType rhomax_in_grad = 1.0;
DType growthfachubble = pow(0.26, 0.55);

//====================================
// Methods below
//====================================
void solve_with_fft(Grid<3,double>& drho, Grid<3,double> &sol);
void assign_to_grid(MultiGrid<_NDIM, DType>& f, std::string filename, std::string store_grid_file);
void compute_v(Grid<_NDIM, DType> phi, DType boxsize);
bool exists(std::string filename);

//==============================================================================
// Read particle-file and assign the particles to grid to compute rho(x,y,z,...)
//==============================================================================
void assign_to_grid(MultiGrid<_NDIM, DType>& f, std::string filename, std::string store_grid_filename){

  unsigned int npart = 0;
  unsigned int N    = f.get_N();
  unsigned int Ntot = f.get_Ntot();
  DType *y = f.get_y(0); 

  std::cout << "Grid N: " << N << " Ntot: " << Ntot << std::endl;

  // Open particle file
  std::ifstream fp;
  fp.open( filename.c_str() );

  fp >> npart;
  std::cout << "Npart: " << npart << std::endl;

  // Velocity grid
  Grid<_NDIM,DType> vel(f.get_grid(0));

  // Assigne particles to grid
  double vrms = 0.0;
  for(unsigned int i = 0; i < npart; i++){
    double xx, yy, zz;
    double vx, vy, vz, v2;

    fp >> xx;
    fp >> yy;
    fp >> zz;
    fp >> vx;
    fp >> vy;
    fp >> vz;

    v2 = vx*vx + vy*vy + vz*vz;
    vrms += v2;

    unsigned int ix = int(xx * N);
    unsigned int iy = int(yy * N);
    unsigned int iz = int(zz * N);
    if(ix == N) ix = 0;
    if(iy == N) iy = 0;
    if(iz == N) iz = 0;

    int index = ix + N*(iy + N*iz);

    vel[index] += sqrt(v2);
    y[index] += 1.0;
  }
  vrms = sqrt(vrms / double(npart));
  std::cout << "Particles read vrms = " << vrms << std::endl;

  // Normalize to delta = rho/rho_mean - 1
  double rhomean = double(npart) / double(Ntot);
  double sum = 0.0, max = -1e30, min = 1e30;
  for(unsigned int i = 0; i < Ntot; i++){
    y[i] = y[i] / rhomean - 1.0;
    if(y[i] > max) max = y[i];
    if(y[i] < min) min = y[i];
    sum += y[i];
  }
  std::cout << "Average: " << sum / double(Ntot) << " Max: " << max << " Min: " << min << std::endl;
  
  // Save grid to file...
  f.get_grid().dump_to_file(store_grid_filename);
}

//===============================================================================
// Compute the velocity from the solution to the PDE (sloopy way; use FFT instead)
//===============================================================================
void compute_v(Grid<_NDIM, DType> phi, DType boxsize){

  unsigned int Ntot = phi.get_Ntot();
  unsigned int N = phi.get_N();

  // This is (1/gridsize * H0 * B0) = (N * 100 * Boxsize_in_Mpc/h)  km/s
  DType norm = DType(N) * 100.0 * boxsize;

  double vrms = 0.0;
  for(unsigned int i = 0; i < Ntot; i++){
    std::vector<unsigned int> coord = phi.index_list(i);
    std::vector<unsigned int> iplus(_NDIM, 0);

    for(unsigned int j = 0; j < _NDIM; j++)
      iplus[j] = coord[j] + 1 < N ? coord[j] + 1 : 0;
   
    // Specialize to NDIM = 3
    unsigned int ixp = phi.grid_index_3d(iplus[0], coord[1], coord[2]);  
    unsigned int iyp = phi.grid_index_3d(coord[0], iplus[1], coord[2]);  
    unsigned int izp = phi.grid_index_3d(coord[0], coord[1], iplus[2]);  

    double v_x2 = std::norm((phi[ixp] - phi[i])*norm);
    double v_y2 = std::norm((phi[iyp] - phi[i])*norm);
    double v_z2 = std::norm((phi[izp] - phi[i])*norm);
    double v2 = v_x2 + v_y2 + v_z2;
    vrms += v2;
  }
  std::cout << "Vrms: " << std::sqrt(vrms/double(Ntot)) << " km/s" << std::endl;
}

//=====================================
// Solve linearized equation using FFT
//=====================================
void solve_with_fft(Grid<3,double>& drho, Grid<3,double> &sol){
  int n = drho.get_N(), nover2 = n/2, n2 = n*n, ntot = drho.get_Ntot();
  double fac = -1.0/(4.0*acos(-1)*acos(-1)) * 1.0/double(n * n * n);
  fftw_complex *out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n*n*n);
  fftw_plan fwd = fftw_plan_dft_3d(n, n, n, out, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan bwd = fftw_plan_dft_3d(n, n, n, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  // Set source
  double rhomean = 0.0;
  for(int i = 0; i < ntot; i++){
    out[i][0] = -growthfachubble * drho[i];
    out[i][1] = 0.0;
    rhomean += out[i][0];
  }
  std::cout << "Rhomean: " << rhomean / double(ntot) << std::endl;

  fftw_execute(fwd);

  // Divide by laplacian
  for(int i=0;i<n;i++){
    int ii = (i < nover2 ? i: i-n);
    for(int j=0;j<n;j++){
      int jj = (j < nover2 ? j: j-n);
      for(int k=0;k<n;k++){
        int kk = (k < nover2 ? k: k-n);
        int ind = i + j*n + k*n2;
        if(ind > 0){
          double facnow = fac/(ii*ii+jj*jj+kk*kk);
          out[ind][0] *=  facnow;
          out[ind][1] *=  facnow;
        }
      }
    }
  }
  out[0][0] = 0.0;
  out[0][1] = 0.0;
  
  fftw_execute(bwd);
  
  // Copy over solution
  for(int i = 0; i < ntot; i++){
    sol[i] = out[i][0];
  }

  fftw_free(out);
}

//==================================
// Check if file exists
//==================================
bool exists(std::string filename){
  std::ifstream fp(filename.c_str());
  return fp.good();
}

int main(){
  int NgridPerDim = 64;
  bool verbose = true;
  DType boxsize = 200.0;
  std::string prefix = "../test_data/rho_";
  std::string grid_filename = prefix + std::to_string(NgridPerDim) + ".dat";
  std::string part_filename = "../test_data/data.ascii";
  MultiGrid<_NDIM, DType> rho;

  // Read grid from file
  if(exists(grid_filename)){

    // Read grid from file 
    Grid<_NDIM, DType> tmp;
    tmp.read_from_file(grid_filename);
    NgridPerDim = tmp.get_N();

    // Make multigrid from grid read from file
    rho = MultiGrid<_NDIM, DType>(tmp);
    
    // Restict down density to all levels
    rho.restrict_down_all();

  } else {

    // Make a density-grid
    rho = MultiGrid<_NDIM, DType>(NgridPerDim);

    // Read particles from file and assign to grid
    assign_to_grid(rho, part_filename, grid_filename);

    // Restict down density to all levels
    rho.restrict_down_all();
  }

  // Solve with fft
  Grid<3,double> fft_sol = rho.get_grid(0);
  solve_with_fft(rho.get_grid(0), fft_sol);
  compute_v(fft_sol, boxsize);

  // Set up solver
  ContinuitySolver<_NDIM, DType> sol(NgridPerDim, verbose);

  // Set some options (rhomin = rhomax = 1 => linear eq. | f(a) H(a))/H0 = growthfachubble) 
  sol.set_rhomin(rhomin_in_grad);
  sol.set_rhomax(rhomax_in_grad);
  sol.set_growthfachubble(growthfachubble);

  sol.set_initial_guess(fft_sol);

  // Add pointer to the source grid to the solver
  sol.add_external_grid(&rho);

  // Set some parameters
  sol.set_ngs_sweeps(4, 10);
  sol.set_epsilon(1e-4);

  // Solve the equation
  sol.solve();

  Grid<3,double> &mgsol = sol.get_grid(0);
  Grid<3,double> err = (mgsol - fft_sol)/sqrt(fft_sol*fft_sol + mgsol*mgsol);
  std::cout << "Difference wrt linear solution: " << err.rms_norm() * 100.0 << " %" << std::endl;

  // Compute v
  compute_v(sol.get_grid(), boxsize);
}


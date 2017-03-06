#include <complex>
#include <fstream>
#include <iostream>
#include <fftw3.h>
#include "poisson_solver.h"

//===================================================
// 
// Example of how to use the multigrid solver 
// We solve the Poisson equation
//
//               D^2 u = S
// 
// for some source (set to give a desired solution)
//
// Select NDIM and DType above and define the two 
// functions f_analytical and laplacian_f_analytical 
//
//===================================================

//==========================
// Dimension of the equation
//==========================
#define _NDIM 1

//=======================================================
// The type of the equation (double,float,complex,...)
//=======================================================
typedef std::complex<double> DType;

//==========================
// Methods below
//==========================
void solve_with_fft(Grid<_NDIM, std::complex<double> >& drho, Grid<_NDIM, std::complex<double> > &sol);
void solve_with_fft(Grid<_NDIM, double>& drho, Grid<_NDIM, double> &sol);
void set_rho( MultiGrid<_NDIM, DType> &rho_mg, Grid<_NDIM, DType> &analytical_solution);
DType f_analytical(double x);
DType laplacian_f_analytical(double x);

//============================================================
// The analytical solution we try to recover
//============================================================

std::complex<double> f_analytical(double x){
  const double pi = std::acos(-1);
  const std::complex<double> I(0, 1);
  return std::exp(2.0 * pi * I * x);
}

std::complex<double> laplacian_f_analytical(double x){
  const double pi = std::acos(-1);
  const std::complex<double> I(0, 1);
  return - 4.0 * pi * pi * std::exp(2.0 * pi * I * x);
}

//============================================================
// Set source rho in D^2 phi = rho analytically
//============================================================

void set_rho( MultiGrid<_NDIM, DType> &rho_mg, Grid<_NDIM, DType> &analytical_solution){
  int N    = rho_mg.get_N();
  int Ntot = rho_mg.get_Ntot();
  DType *rho = rho_mg.get_y(0); 
  DType *yanal = analytical_solution.get_y();

  DType rhomean = 0.0;
  for(int i = 0; i < Ntot; i++){
    double xx = ((i % N) + 0.5) / double(N);

    // Set source rho such that solution to D^2 phi = rho is phi = e^{ 2 pi i x}
    rho[i] = laplacian_f_analytical(xx);
    yanal[i] = f_analytical(xx);
  
    // Box average of rho
    rhomean += rho[i];
  }
  rhomean = rhomean/double(Ntot);
  std::cout << "set_rho :: Box average of source |<rho>| = " << std::sqrt(std::norm(rhomean)) << std::endl;

  // Restrict down source to all levels
  rho_mg.restrict_down_all();
}

//============================================================
// Solve the equation D^2 phi = drho using FFTW
//============================================================

void solve_with_fft(Grid<_NDIM, double>& drho, Grid<_NDIM, double> &sol){
  int Ntot = drho.get_Ntot();
  Grid<_NDIM, std::complex<double> > drho_complex (Ntot);
  Grid<_NDIM, std::complex<double> > sol_complex (Ntot);
  for(int i = 0; i < Ntot; i++){
    drho_complex[i] = drho[i];
    sol_complex[i] = drho[i];
  }
  solve_with_fft(drho_complex, sol_complex);
  for(int i = 0; i < Ntot; i++){
    sol[i] = sol_complex[i].real();
  }
}

void solve_with_fft(Grid<_NDIM, std::complex<double> >& drho, Grid<_NDIM, std::complex<double> > &sol){
  int N = drho.get_N(), Nover2 = N/2, Ntot = drho.get_Ntot();

  fftw_complex *in  = reinterpret_cast<fftw_complex*>(drho.get_y());
  fftw_complex *out = reinterpret_cast<fftw_complex*>(sol.get_y());

#if _NDIM == 1
  fftw_plan fwd = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan bwd = fftw_plan_dft_1d(N, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
#elif _NDIM == 2
  fftw_plan fwd = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan bwd = fftw_plan_dft_2d(N, N, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
#elif _NDIM == 3
  fftw_plan fwd = fftw_plan_dft_3d(N, N, N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan bwd = fftw_plan_dft_3d(N, N, N, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
#else
  std::vector<int> inds(_NDIM, N);
  fftw_plan fwd = fftw_plan_dft(_NDIM, &inds[0], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan bwd = fftw_plan_dft(_NDIM, &inds[0], out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
  
  // Transform
  fftw_execute(fwd);

  // Divide by laplacian and normalize
  double fac = -1.0 / ( 4.0 * std::acos(-1) * std::acos(-1) * double(Ntot) );
  for(int i = 0; i < Ntot; i++){
    double k2 = 0.0;
    int ind = 0;
    for(int j = 0, Npow = 1; j < _NDIM; j++, Npow *= N){
      int coord = i / Npow % N;
      int kval = coord < Nover2 ? coord : coord - N;
      k2 += kval * kval;
      ind += coord * Npow;
    }
    if(i > 0){
      double facnow = fac/k2;
      out[ind][0] *= facnow;
      out[ind][1] *= facnow;
    } else {
      out[0][0] = 0.0;
      out[0][1] = 0.0;
    }
  }

  // Transform back
  fftw_execute(bwd);
}

int main(int argv, char ** argc){

  int NgridPerDim = 256;
  bool verbose = true;
  MultiGrid<_NDIM, DType> rho(NgridPerDim);
  Grid<_NDIM, DType> analytical_solution(NgridPerDim);

  // Make density grid
  set_rho(rho, analytical_solution);

  // Set up solver
  PoissonSolver<_NDIM, DType> sol(NgridPerDim, verbose);

  // Add pointer to the source grid to the solver
  sol.add_external_grid(&rho);

  // Set options
  sol.set_epsilon(1e-10);
  sol.set_ngs_sweeps(2,2);

  // Solve the equation
  sol.solve();

  // Copy over solution
  Grid<_NDIM, DType> solution = sol.get_grid();

  // Compute the fractional error wrt the analytical solution
  Grid<_NDIM, DType> err = (solution - analytical_solution) / (analytical_solution) * std::complex<double>(100.0);
  std::cout << "The rms error in the solution: " << err.rms_norm()  << " %" << std::endl;

  Grid<_NDIM,DType> fft_sol = rho.get_grid(0); 
  solve_with_fft(rho.get_grid(0), fft_sol);

  // Output solution(s)
  std::ofstream fp("test_poisson_solver_complex.txt");
  int N = sol.get_N(0);
  for(int i = 0; i < N; i++){
    fp << (i+0.5)/double(N) << " " << solution[i].real() << " " << fft_sol[i].real() << " " << analytical_solution[i].real() << std::endl;
  }

  // Clean up all memory
  sol.clear();
}


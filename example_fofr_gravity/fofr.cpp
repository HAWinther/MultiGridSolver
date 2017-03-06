#include "fofr_solver.h"
#include <fstream>

//=====================================
//
// This is an example of how to use the
// multigrid solver class for solving
//                                       
//     D[b(u) D f(u)] = g(u,rho)    
//                                       
//         where b(u) = e^u                                                 
// 
//=====================================

//=====================================
// The dim of the equation we solve
//=====================================
#define _NDIM 3

//=====================================
// Particle data container
//=====================================
class Particle{
  public:
    double x, y, z;
    Particle() {};
};

//=====================================
// Methods below
//=====================================
void read_particles(std::string filename, std::vector<Particle> &P);
void assign_to_grid(Grid<_NDIM, double>& f, std::vector<Particle> &P, std::string store_grid_filename);
bool exists(std::string filename);

//=====================================
// Read particles from ascii file
// 
// ------------------------------------
// Format:
// ------------------------------------
// npart
// x1 y1 z1 vx1 vy1 vz1
// x2 y2 z2 vx2 vy2 vz2
// ...
// xn yn zn vxn vyn vzn
// ------------------------------------
//
//=====================================
void read_particles(std::string filename, std::vector<Particle> &P){
  unsigned int npart = 0;
  
  // Check if file exists
  if( ! exists(filename) ){
    std::cout << "File [" << filename << "] does not exist" << std::endl;
    exit(1);
  } else {
    std::cout << "Reading particles from file [" << filename << "]" << std::endl;
  }

  // Open particle file
  std::ifstream fp;
  fp.open( filename.c_str() );

  fp >> npart;
  std::cout << "Npart: " << npart << std::endl;

  // Allocate memory
  P = std::vector<Particle>(npart);

  // Read particles into memory
  for(unsigned int i = 0; i < npart; i++){
    double vx, vy, vz;
    fp >> P[i].x;
    fp >> P[i].y;
    fp >> P[i].z;
    fp >> vx;
    fp >> vy;
    fp >> vz;
  }
  std::cout << "Particles read" << std::endl;
}

//=============================================
// Read particle-file and assign the particles 
// to grid to compute rho(x,y,z,...)
//=============================================
void assign_to_grid(Grid<_NDIM, double>& f, std::vector<Particle> &P, std::string store_grid_filename){
  unsigned int npart = P.size();
  unsigned int N    = f.get_N();
  unsigned int Ntot = f.get_Ntot();
  double *y = f.get_y(); 

  std::cout << "Assigning " << npart << " particles to grid with N: " << N << " Ntot: " << Ntot << std::endl;

  // Assign to grid
  for(unsigned int i = 0; i < npart; i++){
    unsigned int ix = int(P[i].x * N);
    unsigned int iy = int(P[i].y * N);
    unsigned int iz = int(P[i].z * N);
    
    if(ix == N)  ix = 0;
    if(iy == N)  iy = 0;
    if(iz == N)  iz = 0;

    unsigned int index = ix + N * (iy + N * iz);
    
    y[index] += 1.0;
  }

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
  f.dump_to_file(store_grid_filename);
}

//=============================================
// Check if file exists
//=============================================
bool exists(std::string filename){
  std::ifstream fp(filename.c_str());
  return fp.good();
}

int main(int argc, char **argv){
  int NgridPerDim = 64;
  bool verbose = true;
  std::string prefix = "../test_data/rho_";
  std::string grid_filename = prefix + std::to_string(NgridPerDim) + ".dat";
  std::string part_filename = "../test_data/data.ascii";
  MultiGrid<_NDIM, double> rho;

  // Read grid from file or read particles from file and then assign to grid
  if(exists(grid_filename)){

    // Read domain grid from file 
    Grid<_NDIM, double> tmp;
    tmp.read_from_file(grid_filename);

    // Make multigrid from grid read from file
    rho = MultiGrid<_NDIM, double>(tmp);

    // Read data from other levels if availiable
    for(unsigned int i = 1; i < rho.get_Nlevel(); i++){
      int Ntot = rho.get_N(i);
      Grid<_NDIM, double> &r = rho.get_grid(i);

      std::string cur_filename = prefix + std::to_string(Ntot) + ".dat";
      if(exists(cur_filename)) {
        r.read_from_file(prefix + std::to_string(Ntot) + ".dat");
      } else {
        // File do not exist: restict down from upper level
        rho.restrict_down(i-1);
      }
    }

    // Alternative to the above: restict down density to all levels
    // rho.restrict_down_all();

  } else {
    
    // Read particles
    std::vector<Particle> P;
    read_particles(part_filename, P);

    // Allocate the density-grid
    rho = MultiGrid<_NDIM, double>(NgridPerDim);

    // Assign particles to grid
    for(unsigned int i = 0; i < rho.get_Nlevel(); i++) {
      int Ntot = rho.get_N(i);
      assign_to_grid(rho.get_grid(i), P, prefix + std::to_string(Ntot) + ".dat");
    }
  }

  // Set up solver
  FofrSolver<_NDIM, double> sol(NgridPerDim, verbose);

  // Parameters
  double boxsize = 200.0;
  double omegam  = 0.3;
  double aexp    = 1.0;
  double nfofr   = 1.0;
  double fofr0   = 1e-5;
  
  // Set some options in the solver 
  sol.set_parameters(boxsize, omegam, aexp, nfofr, fofr0);
  sol.set_ngs_sweeps(30, 10);
  sol.set_epsilon(1e-8);

  // Make initial guess
  sol.set_initial_guess( log(fofr0) );

  // Add pointer to the source grid to the solver
  sol.add_external_grid(&rho);

  // Solve the equation
  sol.solve();

  // Fetch a reference to the solution
  Grid<_NDIM,double> solution = sol.get_grid(0);

  // Output max and min values for |f_R / f_R0|
  std::cout << "Solution min/max: " << exp(-solution.min()) / fofr0 << "  " << exp(-solution.max()) / fofr0 << std::endl;
}


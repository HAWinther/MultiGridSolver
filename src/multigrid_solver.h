#ifndef _MULTIGRIDSOLVER_HEADER
#define _MULTIGRIDSOLVER_HEADER
#include <bitset>
#include <assert.h>
#include <iostream>
#include <vector>
#include <climits>
#include <iomanip>
#include <complex>
#include "grid.h"
#include "multigrid.h"

    //=============================================
    //                                           //
    // A general multigrid solver to solve       //
    // PDEs in a domain with periodic BC         //
    // in any dimension.                         //
    //                                           //
    // Implement: l_operator, dl_operator        //
    // and (optional) check_convergence          //
    //                                           //
    // The standard equation that is implemented // 
    // is Poissons eq.: D^2phi  = rho, see       //
    // poisson_solver.h for this and how one can //
    // use this class to solve general equations //
    //                                           //
    // Any external fields needed to define the  //
    // PDE can be added to a grid-list           //
    // _ext_fields by running add_external_field // 
    //                                           //
    // _MAXSTEPS defines how many V-cycles we    //
    // do before giving up if convergence is not //
    // reached. Change by running [set_maxsteps] //
    //                                           //
    // _EPS_CONVERGE is a parameters defining    //
    // convergence. In standard implementation   //
    // we define convergence as _rms_res < eps   //
    // Change by running [set_epsilon]           //
    //                                           //
    // _NGRIDCOLOURS defines in which order we   //
    // sweep through the grid: sum of int-coord  //
    // mod _NGRIDCOLOURS. For 2 we have standard //
    // chess-board ordering                      //
    //                                           //
    //===========================================//

template<unsigned int NDIM, typename T>
class MultiGridSolver {
  private:

    unsigned int _N;      // The number of cells per dim in the main grid
    unsigned int _Ntot;   // Total number of cells in the main grid
    unsigned int _Nmin;   // The number of cells per dim in the smallest grid
    unsigned int _Nlevel; // Number of levels
    bool _verbose;        // Turn on verbose while solving

    // All the grids needed
    MultiGrid<NDIM, T> _f;      // The solution
    MultiGrid<NDIM, T> _res;    // The residual
    MultiGrid<NDIM, T> _source; // The multigrid source (restriction of residual)

    // If the source of the equation depends on fields external to the solver they can
    // be added by running add_ext_field and then used in l_operator etc.
    std::vector<MultiGrid<NDIM,T> * > _ext_field;

    // Convergence criterion (if the convergence check is not overwritten)
    bool _conv_criterion_residual = true;  // [True]: residual < eps [False]: residual/residual_i < eps
    double _eps_converge          = 1e-5;  // Convergence criterion
    
    // Solver parameters
    unsigned int _ngs_coarse      = 2;     // Number of NGS sweeps on coarse grid
    unsigned int _ngs_fine        = 2;     // Number of NGS sweeps on the main grid
    unsigned int _maxsteps        = 1000;  // Maximum number of V-cycles
    unsigned int _ngridcolours    = 2;     // The order we go through the grid: 
                                           // [Sum_i coord[i] % ngridcolour == j for j = 0,1,..,ngridcolour-1]
    
    // Book-keeping variables
    unsigned int _istep_vcycle = 0;           // The number of V-cycles we are currenlty at
    unsigned int _tot_sweeps_domain_grid = 0; // Total number of sweeps on the domaingrid (level = 0)

    // Residual information
    double _rms_res;                       // Residual
    double _rms_res_i;                     // The initial residual
    double _rms_res_old;                   // The residual at the old step

    // Internal methods:
    double calculate_residual(unsigned int level, Grid<NDIM,T>& res);
    void   get_neighbor_gridindex(std::vector<unsigned int>& index_list, unsigned int i, unsigned int ngrid);
    void   prolonge_up_array(unsigned int to_level, Grid<NDIM,T>& BottomGrid, Grid<NDIM,T>& TopGrid);
    void   make_prolongation_array(Grid<NDIM,T>& f, Grid<NDIM,T>& Rf, Grid<NDIM,T>& df);
    void   GaussSeidelSweep(unsigned int level, unsigned int curcolor, T *f);
    void   solve_current_level(unsigned int level);
    void   recursive_go_up(unsigned int to_level);
    void   recursive_go_down(unsigned int from_level);
    void   make_new_source(unsigned int level);

  public:

    bool _store_all_residual = false;         // Store the residual after every sweep
    std::vector<double> _res_domain_array;    // Array with the residuals after each step
 
    // Constructors
    MultiGridSolver() {}
    MultiGridSolver(unsigned int N) : MultiGridSolver(N, true) {}
    MultiGridSolver(unsigned int N, bool verbose) : MultiGridSolver(N, 2, verbose) {}
    MultiGridSolver(unsigned int N, unsigned int Nmin, bool verbose);

    // Get a pointer to the solution array / grid
    T *get_y(unsigned int level = 0);
    Grid<NDIM, T> &get_grid(unsigned int level = 0){ return _f.get_grid(level); };

    // Fetch values in externally added fields
    T* get_external_field(unsigned int level, unsigned int field);
    
    // Get values of the multigrid-source used to store the restricted residual during the solve step
    T get_multigrid_source(unsigned int level, unsigned int i);

    // Set precision parameters
    void set_epsilon(double eps_converge);
    void set_maxsteps(unsigned int maxsteps);
    void set_ngs_sweeps(unsigned int ngs_fine, unsigned int ngs_coarse);
    void set_convergence_criterion_residual(bool use_residual);
    
    // Fetch info about the grids
    unsigned int get_N(unsigned int level = 0);
    unsigned int get_Ntot(unsigned int level = 0);

    // Add a pointer to an external grid if this grid is needed to define the PDE
    void add_external_grid(MultiGrid<NDIM,T> *field);

    // Set the initial guess (uniform or from a grid)
    void set_initial_guess(T  guess);
    void set_initial_guess(T *guess);
    void set_initial_guess(Grid<NDIM,T>& guess);

    // Solve the PDE
    void solve();

    // Free up all memory
    void clear();

    // Methods that must be implemented by user
    virtual T l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource);
    virtual T dl_operator(unsigned int level, std::vector<unsigned int>& index_list);
    virtual bool   check_convergence();
};

#endif

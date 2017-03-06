#ifndef _POISSONSOLVER_HEADER
#define _POISSONSOLVER_HEADER
#include "../src/multigrid_solver.h"
    
    //=============================================
    //                                           //
    //  An example of how to use the multigrid   //
    //  solver. We implement the Poisson like eq //
    //                                           //
    //         D^2 phi = m^2 phi + C * drho      //
    //                                           //
    // T can be any float-type: float, double,   //
    //  complex<float>, complex<double>, ...     //
    //                                           //
    // NB: we need <drho> = 0.0 for consistency  //
    //=============================================

template<unsigned int NDIM, typename T>
class PoissonSolver : public MultiGridSolver<NDIM,T> {
  private:

    T _mass = 0.0;    // The mass-term 'm'
    T _rhofac = 1.0;  // The prefactor 'C' to drho

  public:
  
    // Constructors
    PoissonSolver() {}
    PoissonSolver(unsigned int N, bool verbose = true) : MultiGridSolver<NDIM,T>(N, verbose) {} ;
    PoissonSolver(unsigned int N, unsigned int Nmin, bool verbose = true) : MultiGridSolver<NDIM,T>(N, Nmin, verbose) {} ;
    PoissonSolver(MultiGrid<NDIM,T> &source, bool verbose = true) : MultiGridSolver<NDIM,T>(source.get_N(), source.get_Nmin(), verbose) {
      MultiGridSolver<NDIM,T>::add_external_grid(&source);
    }

    // Change the mass-term
    void set_mass(double mass){
      _mass = mass;
    }
    
    // Change the rho-factor
    void set_rhofac(double rhofac){
      _rhofac = rhofac;
    }
    
    // The discretized equation of motion L_{ijk...}(phi) (= 0 when phi is the correct solution)
    T l_operator(unsigned int level, std::vector<int>& index_list, bool addsource){ 
      unsigned int i = index_list[0];
      
      // Gridspacing
      T h = 1.0/T( MultiGridSolver<NDIM,T>::get_N(level) ); 

      // Solution and pde-source grid at current level
      T *phi  = MultiGridSolver<NDIM,T>::get_y(level);
      T *drho = MultiGridSolver<NDIM,T>::get_external_field(level, 0);

      // Compute the standard kinetic term [D^2 phi] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
      T kinetic = -2.0 * NDIM * phi[ i ];
      for(unsigned int k = 1; k < 2*NDIM + 1; k++){
        kinetic += phi[ index_list[k] ];
      }

      // The right hand side of the PDE with the source term arising 
      // from restricting the equation down to the lower level
      T source = _mass * _mass * phi[ i ] + _rhofac * drho[ i ];
      if( level > 0 && addsource ){
        source += MultiGridSolver<NDIM,T>::get_multigrid_source(level, i);
      }

      return kinetic/(h*h) - source;
    }

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(unsigned int level, std::vector<int>& index_list){ 
      T h = 1.0/T( MultiGridSolver<NDIM,T>::get_N(level) );
     
      // Derivative of kinetic term
      T dkinetic = -2.0*NDIM;

      // Derivative of source
      T dsource = _mass * _mass;

      T dl = dkinetic/(h*h) - dsource;
      return dl;
    }
};

#endif

#ifndef _FOFRSOLVER_HEADER
#define _FOFRSOLVER_HEADER
#include "../src/multigrid_solver.h"
    
    //=============================================
    //                                           //
    //  An example of how to use the multigrid   //
    //  solver. We implement the Hu-Sawicky      //
    //  f(R) gravity equation:                   //
    //                                           //
    //          D[b(u) D f(u)] = g(u,rho)        //
    //                                           //
    //              where b(u) = e^u             //                                         
    //=============================================

template<unsigned int NDIM, typename T>
class FofrSolver : public MultiGridSolver<NDIM,T> {
  private:

    // Parameters
    double _boxsize = 100.0;    // Boxsize in units of Mpc/h
    double _omegam  = 0.3;      // Density parameter Omega_matter
    double _aexp    = 1.0;      // Scale factor
    double _nfofr   = 1.0;      // Hu-Sawicky f(R) model n-index
    double _fofr0   = 1.0e-5;   // Hu-Sawicky f(R) model f(R0) value

    // Derived parameters
    double _fac1    = 10.33333333;
    double _fac2    = 0.032676900;
    double _prefac  = 0.000333778;

  public:

    // Constructors
    FofrSolver() {}
    FofrSolver(unsigned int N, bool verbose = true) : MultiGridSolver<NDIM,T>(N, verbose) {} ;
    FofrSolver(unsigned int N, int Nmin, bool verbose = true) : MultiGridSolver<NDIM,T>(N, Nmin, verbose) {} ;
    FofrSolver(MultiGrid<NDIM,T> &source, bool verbose = true) : MultiGridSolver<NDIM,T>(source.get_N(), source.get_Nmin(), verbose) {
      MultiGridSolver<NDIM,T>::add_source_grid(&source);
    }

    void set_parameters(double boxsize, double omegam, double aexp, double nfofr, double fofr0){
      std::cout << "==> FofrSolver::set_parameters | Boxsize: " << boxsize << " Omegam: " << omegam;
      std::cout << " aexp: " << aexp    << " nfofr: " << nfofr << " fofr0: " << fofr0 << std::endl;

      // Set parameters
      _boxsize = boxsize;
      _omegam  = omegam;
      _aexp    = aexp;
      _nfofr   = nfofr;
      _fofr0   = fofr0;
      
      // Derived parameters
      _fac1 = 1.0 + 4.0 * _aexp * _aexp * _aexp * (1.0 - _omegam) / _omegam;
      _fac2 = (1.0 + 4.0 * (1.0 - _omegam) / _omegam ) * _aexp * _aexp * _aexp * std::pow( _fofr0 * _aexp, 1.0 / (_nfofr+1.0) );
      _prefac  = pow( _boxsize / 2998.0 , 2) * _omegam * _aexp;
    }

    inline double bfunc(double x){ return std::exp(x); }

    // The dicretized equation L(phi)
    T  l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource){ 
      unsigned int i = index_list[0];
      
      // Gridspacing
      T h = 1.0/T( MultiGridSolver<NDIM,T>::get_N(level) ); 

      // Solution and pde-source grid at current level
      T *phi = MultiGridSolver<NDIM,T>::get_y(level);
      T *drho = MultiGridSolver<NDIM,T>::get_external_field(level, 0);
 
      // The kinetic term is D[ b[f] Df ] with b = exp(x)
      T kinetic = 0.0;
      T f  = phi[i];
      T boff = bfunc(f);
      for(unsigned int j = 0; j < NDIM; j++){
        unsigned int iminus = index_list[2*j+1];
        unsigned int iplus  = index_list[2*j+2];
        T fminus = phi[iminus];
        T fplus  = phi[iplus];
        T bminus = 0.5 * ( bfunc(fminus) + boff );
        T bplus  = 0.5 * ( bfunc(fplus)  + boff );
        kinetic += bplus * (fplus - f) - bminus * (f - fminus);
      }

      // The right hand side of the PDE 
      T source = _prefac * ( drho[ i ] + _fac1 - _fac2 * std::exp(-phi[ i ]/(_nfofr+1.0)) );
      if( level > 0 && addsource ){
        source += MultiGridSolver<NDIM,T>::get_multigrid_source(level, i);
      }

      // The discretized equation of motion L_{ijk...}(phi) = 0
      return kinetic/(h*h) - source;
    }

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(unsigned int level, std::vector<unsigned int>& index_list){ 
      unsigned int i = index_list[0];
      
      // Gridspacing
      T h = 1.0/T( MultiGridSolver<NDIM,T>::get_N(level) );
      
      // Solution and pde-source grid at current level
      T *phi = MultiGridSolver<NDIM,T>::get_y(level);

      // The derivative of the kinetic term D[ b[f] Df ]
      T dkinetic = 0.0;
      T f        = phi[i];
      T boff     = bfunc(f);
      T dboff    = boff;
      for(unsigned int j = 0; j < NDIM; j++){
        unsigned int iminus = index_list[2*j+1];
        unsigned int iplus  = index_list[2*j+2];
        T fminus = phi[iminus];
        T fplus  = phi[iplus];
        T bminus = 0.5 * ( bfunc(fminus) + boff );
        T bplus  = 0.5 * ( bfunc(fplus)  + boff );
        
        // Add up kinetic term
        dkinetic += 0.5 * dboff * (fplus + fminus - 2.0 * f) - (bplus + bminus);
      }
     
      // The derivative of the right hand side
      T dsource = _prefac * ( _fac2 * std::exp(-phi[ i ] / (_nfofr+1.0) ) ) / (1.0+_nfofr);
        
      return dkinetic/(h*h) - dsource;
    }
};

#endif

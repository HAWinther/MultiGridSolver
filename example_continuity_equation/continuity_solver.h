#ifndef _CONTINUITYSOLVER_HEADER
#define _CONTINUITYSOLVER_HEADER
#include "../src/multigrid_solver.h"
    
    //=============================================
    //  Solving the continuity equation:         //
    //                                           //
    //       fH delta + D[(1+delta) v] = 0       //
    //                                           //
    // where delta' = fHdelta (linear approx)    //
    // and v = H0 B0 D_codephi) is the velocity  // 
    // field.                                    //
    //                                           //
    // Define f(a) H(a)/H0 by setting            //
    // _growthfacHubble                          //
    //                                           //
    // To avoid dividing by 0 we restrict        //
    // (1+delta) >= _rhomin                      //
    //                                           //
    // Taking _rhomax = _rhomin = 1.0 we get     //
    // the completely linear equation            //
    //=============================================

template<unsigned int NDIM, typename T>
class ContinuitySolver : public MultiGridSolver<NDIM,T> {
  private:

    double _rhomin = 0.1;
    double _rhomax = 1e30;
    double _growthfachubble = 0.5;

  public:
  
    ContinuitySolver() {}
    ContinuitySolver(unsigned int N, bool verbose = true) : MultiGridSolver<NDIM,T>(N, verbose) {} ;
    ContinuitySolver(unsigned int N, unsigned int Nmin, bool verbose = true) : MultiGridSolver<NDIM,T>(N, Nmin, verbose) {} ;
    ContinuitySolver(MultiGrid<NDIM,T> &source, bool verbose = true) : MultiGridSolver<NDIM,T>(source.get_N(), source.get_Nmin(), verbose) {
      MultiGridSolver<NDIM,T>::add_external_grid(&source);
    }

    void set_rhomax(double rhomax){
      _rhomax = rhomax;
    }
   
    void set_rhomin(double rhomin){
      _rhomin = rhomin;
    }
    
    void set_growthfachubble(double growthfachubble){ 
      _growthfachubble = growthfachubble; 
    }
    
    T bfunc(T delta_rho){
      return std::min( std::max(1.0 + delta_rho, _rhomin) , _rhomax);
    }

    // Define the PDE. l_operator is the discretized EOM [ L(phi) = 0 ]
    T l_operator(unsigned int level, std::vector<unsigned int>& index_list, bool addsource){ 
      unsigned int i = index_list[0];
    
      // Gridspacing
      T h = 1.0/T( MultiGridSolver<NDIM,T>::get_N(level) ); 

      // Solution and pde-source grid at current level
      T *phi = MultiGridSolver<NDIM,T>::get_y(level);
      T *delta = MultiGridSolver<NDIM,T>::get_external_field(level, 0);

      // Compute the kinetic term D[b Df] (in 1D this is phi''_i =  phi_{i+1} + phi_{i-1} - 2 phi_{i} )
      T kinetic = 0.0;
      T f    = phi[i];
      T boff = bfunc(delta[i]);
      for(unsigned int j = 0; j < NDIM; j++){
        unsigned int iminus = index_list[2*j+1];
        unsigned int iplus  = index_list[2*j+2];
        T fminus = phi[iminus];
        T fplus  = phi[iplus];
        T bminus = 0.5 * ( bfunc(delta[iminus]) + boff );
        T bplus  = 0.5 * ( bfunc(delta[iplus])  + boff );
        kinetic += bplus * (fplus - f) - bminus * (f - fminus);
      }

      // The right hand side of the PDE with the source term arising 
      // from restricting the equation down to the lower level
      T source = -_growthfachubble * delta[i];
      if( level > 0 && addsource ){
        source += MultiGridSolver<NDIM,T>::get_multigrid_source(level, i);
      }

      // The discretized equation of motion L_{ijk...}(phi)
      return kinetic/(h*h) - source;
    }

    // Differential of the L operator: dL_{ijk...}/dphi_{ijk...}
    T dl_operator(unsigned int level, std::vector<unsigned int>& index_list){ 
      unsigned int i = index_list[0];
      
      // Gridspacing
      T h = 1.0/T( MultiGridSolver<NDIM,T>::get_N(level) );
      
      // Solution and pde-source grid at current level
      T *delta = MultiGridSolver<NDIM,T>::get_external_field(level, 0);
  
      // Derivative of the kinetic term
      T dkinetic = 0.0;
      T boff = bfunc(delta[i]);
      for(unsigned int j = 0; j < NDIM; j++){
        unsigned int iminus = index_list[2*j+1];
        unsigned int iplus  = index_list[2*j+2];
        T bminus = 0.5 * ( bfunc(delta[iminus]) + boff );
        T bplus  = 0.5 * ( bfunc(delta[iplus])  + boff );
        dkinetic += -(bplus + bminus);
      }
   
      // The derivative of the source is zero as there is no field dependence
      T dsource = 0.0;

      // Check for error
      T dl = dkinetic/(h*h) - dsource;
      if(fabs(dl) < 1e-5){
        std::cout << "Error: dl = " << dl << " delta_rho = " << delta[i] << " b = " << boff << std::endl;  
        exit(1);
      }

      return dl;
    }
};

#endif

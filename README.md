# MultiGridSolver

A multigrid solver for general PDEs in [0,1]^n with periodic boundary conditions.

Author: Hans A. Winther, University of Oxford (2016)

 - Simple multigrid-solver for general (linear or non-linear) PDEs in any dimension.

 - Templated on the type: float, double, complex<double>, etc. 

 - The boundary conditions implemented are periodic. Can be extended with not too much work, but have not added this yet (e.g. add a mask field)

 - Made to make it easy to implement new equations by defining a new class that extends the MG-solver class.

 - Make the class [MysSolver : public MultiGridSolver] and implement the functions [l_operator] and [dl_operator] 

 - One can also implement the convergence criterion in [check_convergence] (rms-residual < epsilon is standard). 

 - Any external grids needed to define the PDE can be added through [add_external_grid]

Examples included:
 -  Poisson solver     D^2 phi = S
 -  f(R) solver        D[b(phi) D phi] = c(phi)
 -  Continuity solver  f delta + D[(1+delta]v] = 0 for velocity field v


#ifndef GLOBALS_H
#define GLOBALS_H
extern int  verbose;     
extern Real tolerance_rel;
extern Real tolerance_abs;
extern int  maxiter; 
extern int  plot_rhs;
extern int  plot_beta;
extern int  plot_soln;
extern int  plot_asol;
extern int  plot_err;
extern int  comp_norm;

extern Real dx[BL_SPACEDIM];

enum solver_t {BoxLib_C, BoxLib_F, Hypre, HPGMG, All};
enum bc_t {Periodic = 0,
	   Dirichlet = LO_DIRICHLET, 
	   Neumann = LO_NEUMANN};
extern solver_t solver_type;
extern bc_t     bc_type;
extern int domain_boundary_condition;

extern int Ncomp;
#endif

#include <globals.H>

void solve(MultiFab& soln, const MultiFab& anaSoln, 
	   Real a, Real b, MultiFab& alpha, MultiFab beta[], 
	   MultiFab& rhs, const BoxArray& bs, const Geometry& geom,
	   solver_t solver);

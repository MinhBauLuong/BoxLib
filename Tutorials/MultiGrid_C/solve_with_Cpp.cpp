#include <MultiFab.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LO_BCTYPES.H>
#include <BndryData.H>
#include <ABecLaplacian.H>
#include <MultiGrid.H>

#include <globals.H>
#include <set_boundary.H>

void solve_with_Cpp(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BL_PROFILE("solve_with_Cpp()");
  BndryData bd(bs, 1, geom);
  set_boundary(bd, rhs, 0);

  ABecLaplacian abec_operator(bd, dx);
  abec_operator.setScalars(a, b);
  abec_operator.setCoefficients(alpha, beta);

  MultiGrid mg(abec_operator);
  mg.setVerbose(verbose);
  mg.solve(soln, rhs, tolerance_rel, tolerance_abs);
} 

#ifdef USEHYPRE
void solve_with_hypre(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		      MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BL_PROFILE("solve_with_hypre()");
  BndryData bd(bs, 1, geom);
  set_boundary(bd, rhs);

  HypreABecLap hypreSolver(bs, geom);
  hypreSolver.setScalars(a, b);
  hypreSolver.setACoeffs(alpha);
  hypreSolver.setBCoeffs(beta);
  hypreSolver.setVerbose(verbose);
  hypreSolver.solve(soln, rhs, tolerance_rel, tolerance_abs, maxiter, bd);
}
#endif

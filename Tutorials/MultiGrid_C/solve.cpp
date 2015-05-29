#include <MultiFab.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LO_BCTYPES.H>

#include <globals.H>
#include <solve_with_Cpp.H>
#include <solve_with_F90.H>
#include <writePlotFile.H>

void solve(MultiFab& soln, const MultiFab& anaSoln, 
	   Real a, Real b, MultiFab& alpha, MultiFab beta[], 
	   MultiFab& rhs, const BoxArray& bs, const Geometry& geom,
	   solver_t solver)
{
  BL_PROFILE("solve()");
  std::string ss;

  soln.setVal(0.0);

  const Real run_strt = ParallelDescriptor::second();

  if (solver == BoxLib_C) {
    ss = "CPP";
    solve_with_Cpp(soln, a, b, alpha, beta, rhs, bs, geom);
  }
#ifdef USE_F90_SOLVERS
  else if (solver == BoxLib_F) {
    ss = "F90";
    solve_with_F90(soln, a, b, alpha, beta, rhs, bs, geom);
  }
#endif
#ifdef USEHYPRE
  else if (solver == Hypre) {
    ss = "Hyp";
    solve_with_hypre(soln, a, b, alpha, beta, rhs, bs, geom);
  }
#endif
  else {
    BoxLib::Error("Invalid solver");
  }

  Real run_time = ParallelDescriptor::second() - run_strt;

  ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "   Run time      : " << run_time << std::endl;
  }

  if (plot_soln) {
    writePlotFile("SOLN-"+ss, soln, geom);
  }

  if (plot_err || comp_norm) {
    soln.minus(anaSoln, 0, Ncomp, 0); // soln contains errors now
    MultiFab& err = soln;

    if (plot_err) {
      writePlotFile("ERR-"+ss, soln, geom);
    }
    
    if (comp_norm) {
      Real twoNorm = err.norm2();
      Real maxNorm = err.norm0();
      
      err.setVal(1.0);
      Real vol = err.norm2();
      twoNorm /= vol;
      
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "   2 norm error  : " << twoNorm << std::endl;
	std::cout << "   max norm error: " << maxNorm << std::endl;
      }
    }
  }
}

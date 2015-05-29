#include <MultiFab.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <LO_BCTYPES.H>

#include <writePlotFile.H>
#include <RHS_F.H>
#include <globals.H>

void setup_rhs(MultiFab& rhs, const Geometry& geom)
{
  BL_PROFILE("setup_rhs()");
  ParmParse pp;

  Real a, b, sigma, w;
  pp.get("a",  a);
  pp.get("b",  b);
  pp.get("sigma", sigma);
  pp.get("w", w);

  int ibnd = static_cast<int>(bc_type);

  // We test the sum of the RHS to check solvability
  Real sum_rhs = 0.;

  for ( MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
    const int* rlo = rhs[mfi].loVect();
    const int* rhi = rhs[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_SET_RHS(rhs[mfi].dataPtr(),ARLIM(rlo),ARLIM(rhi),
                 bx.loVect(),bx.hiVect(),dx, a, b, sigma, w, ibnd);
  }

  for ( MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
    sum_rhs += rhs[mfi].sum(0,1);
  }

  for (int n=0; n < BL_SPACEDIM; n++) {
    sum_rhs *= dx[n];
  }

  ParallelDescriptor::ReduceRealSum(sum_rhs, ParallelDescriptor::IOProcessorNumber());

  if (ParallelDescriptor::IOProcessor()) {
     std::cout << "Sum of RHS      : " << sum_rhs << std::endl;
  }

  if (plot_rhs == 1) {
    writePlotFile("RHS", rhs, geom);
  }
}

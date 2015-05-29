#include <MultiFab.H>
#include <LO_BCTYPES.H>
#include <ArrayLim.H>

#include <globals.H>
#include <COEF_F.H>

void compute_analyticSolution(MultiFab& anaSoln)
{
  BL_PROFILE("compute_analyticSolution()");
  int ibnd = static_cast<int>(bc_type); 

  for (MFIter mfi(anaSoln); mfi.isValid(); ++mfi) {
    const int* alo = anaSoln[mfi].loVect();
    const int* ahi = anaSoln[mfi].hiVect();
    const Box& bx = mfi.validbox();

    FORT_COMP_ASOL(anaSoln[mfi].dataPtr(), ARLIM(alo), ARLIM(ahi),
		   bx.loVect(),bx.hiVect(),dx, ibnd);
  }
}


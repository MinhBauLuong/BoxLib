#include <BoxArray.H>
#include <MultiFab.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <LO_BCTYPES.H>
#include <ArrayLim.H>

#include <globals.H>
#include <COEF_F.H>
#include <writePlotFile.H>

void setup_coeffs(BoxArray& bs, MultiFab& alpha, MultiFab beta[], const Geometry& geom, MultiFab& cc_coef)
{
  BL_PROFILE("setup_coeffs()");
  ParmParse pp;

  Real sigma, w;
  pp.get("sigma", sigma);
  pp.get("w", w);

  for ( MFIter mfi(alpha); mfi.isValid(); ++mfi ) {
    const Box& bx = mfi.validbox();

    const int* alo = alpha[mfi].loVect();
    const int* ahi = alpha[mfi].hiVect();
    FORT_SET_ALPHA(alpha[mfi].dataPtr(),ARLIM(alo),ARLIM(ahi),
    		   bx.loVect(),bx.hiVect(),dx);

    const int* clo = cc_coef[mfi].loVect();
    const int* chi = cc_coef[mfi].hiVect();
    FORT_SET_CC_COEF(cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),
		     bx.loVect(),bx.hiVect(),dx, sigma, w);
  }

  // convert cell-centered beta to edges
  for ( int n=0; n<BL_SPACEDIM; ++n ) {
    for ( MFIter mfi(beta[n]); mfi.isValid(); ++mfi ) {
      int i = mfi.index();
      Box bx(bs[i]);
      const int* clo = cc_coef[mfi].loVect();
      const int* chi = cc_coef[mfi].hiVect();
      const int* edgelo = beta[n][mfi].loVect();
      const int* edgehi = beta[n][mfi].hiVect();

      FORT_COEF_TO_EDGES(&n,beta[n][mfi].dataPtr(),ARLIM(edgelo),ARLIM(edgehi),
			 cc_coef[mfi].dataPtr(),ARLIM(clo),ARLIM(chi),
			 bx.loVect(),bx.hiVect());
    }
  }

  if (plot_beta == 1) {
    writePlotFile("BETA", cc_coef, geom);
  }
}


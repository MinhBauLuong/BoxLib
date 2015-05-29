#ifdef USE_F90_SOLVERS
#include <MultiFab.H>
#include <Geometry.H>
#include <LO_BCTYPES.H>
#ifdef USE_F90_SOLVERS
#include <MGT_Solver.H>
#include <stencil_types.H>
#endif
#include <MacBndry.H>

#include <globals.H>

void solve_with_F90(MultiFab& soln, Real a, Real b, MultiFab& alpha, MultiFab beta[], 
		    MultiFab& rhs, const BoxArray& bs, const Geometry& geom)
{
  BL_PROFILE("solve_with_F90()");
  // Translate into F90 solver
  std::vector<BoxArray> bav(1);
  bav[0] = bs;
  std::vector<DistributionMapping> dmv(1);
  dmv[0] = rhs.DistributionMap();
  std::vector<Geometry> fgeom(1);
  fgeom[0] = geom;

  int mg_bc[2*BL_SPACEDIM];

  Array< Array<Real> > xa(1);
  Array< Array<Real> > xb(1);

  xa[0].resize(BL_SPACEDIM);
  xb[0].resize(BL_SPACEDIM);

  // Set the boundary conditions to live exactly on the faces of the domain --
  //   this is only relevant for Dirichlet and Neumann bc's
  for ( int n = 0; n < BL_SPACEDIM; ++n ) {
    xa[0][n] = 0.;
    xb[0][n] = 0.;
  }

  if (bc_type == Periodic) {
    // Define the type of boundary conditions to be periodic
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_PER;
      mg_bc[2*n + 1] = MGT_BC_PER;
    }
  }
  else if (bc_type == Neumann) {
    // Define the type of boundary conditions to be Neumann
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_NEU;
      mg_bc[2*n + 1] = MGT_BC_NEU;
    }
  }
  else if (bc_type == Dirichlet) {
    // Define the type of boundary conditions to be Dirichlet
    for ( int n = 0; n < BL_SPACEDIM; ++n ) {
      mg_bc[2*n + 0] = MGT_BC_DIR;
      mg_bc[2*n + 1] = MGT_BC_DIR;
    }
  }

  MGT_Solver mgt_solver(fgeom, mg_bc, bav, dmv, false, CC_CROSS_STENCIL);

  MultiFab* soln_p[1]; soln_p[0] = &soln;
  MultiFab*  rhs_p[1];  rhs_p[0] = &rhs;

  PArray<MultiFab> acoeffs(1, PArrayManage);
  acoeffs.set(0, new MultiFab(bs,Ncomp,0,Fab_allocate));
  acoeffs[0].copy(alpha);
  acoeffs[0].mult(a); 

  Array< PArray<MultiFab> > bcoeffs(BL_SPACEDIM);
  for (int n = 0; n < BL_SPACEDIM ; n++) {
    bcoeffs[n].resize(1,PArrayNoManage);
    bcoeffs[n].set(0, &beta[n]);
  }

  // The coefficients are set such that we will solve
  //  (a alpha - b del dot beta grad) soln = rhs
  //  written in the form 
  //  (acoeffs - b del dot bcoeffs grad) soln = rhs
  mgt_solver.set_visc_coefficients(acoeffs,bcoeffs,b,xa,xb);

  BCRec phys_bc;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    phys_bc.setLo(n,0);
    phys_bc.setHi(n,0);
  }

  MacBndry bndry(bs,1,geom);
  bndry.setBndryValues(soln,0,0,1,phys_bc);
  
  Real final_resnorm;
  mgt_solver.solve(soln_p, rhs_p, tolerance_rel, tolerance_abs, bndry, final_resnorm);
}
#endif

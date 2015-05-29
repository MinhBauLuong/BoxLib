// We solve (a alpha - b del dot beta grad) soln = rhs
// where a and b are scalars, alpha and beta are arrays

#include <fstream>
#include <iomanip>

#include <Utility.H>
#include <ParmParse.H>
#include <LO_BCTYPES.H>
#include <BndryData.H>
#include <MultiGrid.H>
#include <CGSolver.H>
#include <Laplacian.H>
#include <ABecLaplacian.H>
#include <ParallelDescriptor.H>
#include <MacBndry.H>
#ifdef USE_F90_SOLVERS
#include <MGT_Solver.H>
#include <stencil_types.H>
#endif

#ifdef USEHYPRE
#include <HypreABecLap.H>
#endif

#include <COEF_F.H>
#include <RHS_F.H>
#include <writePlotFile.H>
#include <compute_analyticSolution.H>
#include <setup_coeffs.H>
#include <setup_rhs.H>
#include <set_boundary.H>
#include <solve.H>
#include <solve_with_Cpp.H>
#ifdef USE_F90_SOLVERS
#include <solve_with_F90.H>
#endif
#ifdef USEHYPRE
#include <solve_with_hypre.H>
#endif
#ifdef USEHPGMG
#include <solve_with_HPGMG.H>
#endif
#include <globals.H>

int verbose       = 2;
Real tolerance_rel = 1.e-10; // converged if ||b-Ax|| / ||b|| < tolerance_rel
Real tolerance_abs = 0.0; // converged if ||D^{-1}(b-Ax)|| < tolerance_abs
int maxiter       = 100;
int plot_rhs      = 0;
int plot_beta     = 0;
int plot_soln     = 0;
int plot_asol     = 0;
int plot_err      = 0;
int comp_norm     = 1;
int Ncomp = 1;

Real dx[BL_SPACEDIM];

solver_t solver_type;
bc_t     bc_type;
int      domain_boundary_condition;

int main(int argc, char* argv[])
{
  BoxLib::Initialize(argc,argv);

  BL_PROFILE_VAR("main()", pmain);

  std::cout << std::setprecision(15);

  ParmParse ppmg("mg");
  ppmg.query("v", verbose);

  ParmParse pp;

  {
    std::string solver_type_s;
    pp.get("solver_type",solver_type_s);
    if (solver_type_s == "BoxLib_C") {
      solver_type = BoxLib_C;
    }
    else if (solver_type_s == "BoxLib_F") {
#ifdef USE_F90_SOLVERS
      solver_type = BoxLib_F;
#else
      BoxLib::Error("Set USE_FORTRAN=TRUE in GNUmakefile");
#endif
    }
    else if (solver_type_s == "Hypre") {
#ifdef USEHYPRE
      solver_type = Hypre;
#else
      BoxLib::Error("Set USE_HYPRE=TRUE in GNUmakefile");
#endif
    }
    else if (solver_type_s == "All") {
      solver_type = All;
    }
    else {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Don't know this solver type: " << solver_type << std::endl;
      }
      BoxLib::Error("");
    }
  }

  {
    std::string bc_type_s;
    pp.get("bc_type",bc_type_s);
    if (bc_type_s == "Dirichlet") {
      bc_type = Dirichlet;
      domain_boundary_condition = BC_DIRICHLET;
    }
    else if (bc_type_s == "Neumann") {
      bc_type = Neumann;
#ifdef USEHPGMG
      BoxLib::Error("HPGMG does not support Neumann boundary conditions");
#endif
    }
    else if (bc_type_s == "Periodic") {
      bc_type = Periodic;
      domain_boundary_condition = BC_PERIODIC;
    }
    else {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Don't know this boundary type: " << bc_type << std::endl;
      }
      BoxLib::Error("");
    }
  }

  pp.query("tol_rel", tolerance_rel);
  pp.query("tol_abs", tolerance_abs);
  pp.query("maxiter", maxiter);
  pp.query("plot_rhs" , plot_rhs);
  pp.query("plot_beta", plot_beta);
  pp.query("plot_soln", plot_soln);
  pp.query("plot_asol", plot_asol);
  pp.query("plot_err", plot_err);
  pp.query("comp_norm", comp_norm);

  Real a, b;
  pp.get("a",  a);
  pp.get("b",  b);

  int n_cell;
  int max_grid_size;

  pp.get("n_cell",n_cell);
  pp.get("max_grid_size",max_grid_size);

  // Define a single box covering the domain
  IntVect dom_lo(D_DECL(0,0,0));
  IntVect dom_hi(D_DECL(n_cell-1,n_cell-1,n_cell-1));
  Box domain(dom_lo,dom_hi);

  // Initialize the boxarray "bs" from the single box "bx"
  BoxArray bs(domain);

  // Break up boxarray "bs" into chunks no larger than "max_grid_size" along a
  // direction
  bs.maxSize(max_grid_size);

  // This defines the physical size of the box.  Right now the box is [0,1] in
  // each direction.
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++) {
    real_box.setLo(n, 0.0);
    real_box.setHi(n, 1.0);
  }

  // This says we are using Cartesian coordinates
  int coord = 0;

  // This sets the boundary conditions to be periodic or not
  int is_per[BL_SPACEDIM];

  if (bc_type == Dirichlet || bc_type == Neumann) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Using Dirichlet or Neumann boundary conditions." << std::endl;
    }
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 0;
  }
  else {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "Using periodic boundary conditions." << std::endl;
    }
    for (int n = 0; n < BL_SPACEDIM; n++) is_per[n] = 1;
  }

  // This defines a Geometry object which is useful for writing the plotfiles
  Geometry geom(domain,&real_box,coord,is_per);

  for ( int n=0; n<BL_SPACEDIM; n++ ) {
    dx[n] = (geom.ProbHi(n)-geom.ProbLo(n))/domain.length(n);
  }

  if (ParallelDescriptor::IOProcessor()) {
     std::cout << "Domain size     : " << n_cell << std::endl;
     std::cout << "Max_grid_size   : " << max_grid_size << std::endl;
     std::cout << "Number of grids : " << bs.size() << std::endl;
  }

  // Allocate and define the right hand side.
  MultiFab rhs(bs, Ncomp, 0, Fab_allocate);
  setup_rhs(rhs, geom);

  int my_rank = 0, num_ranks = 1;

  // Set up the Helmholtz operator coefficients.
  MultiFab alpha(bs, Ncomp, 0, Fab_allocate);
  MultiFab beta[BL_SPACEDIM];
  for (int n=0; n<BL_SPACEDIM; ++n) {
    BoxArray bx(bs);
    beta[n].define(bx.surroundingNodes(n), Ncomp, 1, Fab_allocate);
  }
  MultiFab beta_cc(bs,Ncomp,1); // cell-centered beta
  setup_coeffs(bs, alpha, beta, geom, beta_cc);

  BndryData bd(bs, 1, geom);
  set_boundary(bd, rhs, 0);

#ifdef USEHPGMG
  int minCoarseDim;
  if (domain_boundary_condition == BC_PERIODIC)
  {
    minCoarseDim = 2; // avoid problems with black box calculation of D^{-1} for poisson with periodic BC's on a 1^3 grid
  }
  else
  {
    minCoarseDim = 1; // assumes you can drop order on the boundaries
  }

  level_type level_h;
  mg_type MG_h;
  int numVectors = 12;
#ifdef BL_USE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &num_ranks);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
#endif /* BL_USE_MPI */

  const double h0 = dx[0];
  // Create the geometric structure of the HPGMG grid using the RHS MultiFab as
  // a template. This doesn't copy any actual data.
  MultiFab::CreateHPGMGLevel(&level_h, rhs, n_cell, max_grid_size, my_rank, num_ranks, domain_boundary_condition, numVectors, dx[0]);

  // Set up the coefficients for the linear operator L.
  MultiFab::SetupHPGMGCoefficients(a, b, alpha, beta_cc, &level_h);

  // Now that the HPGMG grid is built, populate it with RHS data.
  MultiFab::ConvertToHPGMGLevel(rhs, n_cell, max_grid_size, &level_h, VECTOR_F);

#ifdef USE_HELMHOLTZ
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Creating Helmholtz (a=" << a << ", b=" << b << ") test problem" << std::endl;;
  }
#else
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Creating Poisson (a=" << a << ", b=" << b << ") test problem" << std::endl;;
  }
#endif /* USE_HELMHOLTZ */
  const int boxes_in_i = n_cell / max_grid_size;

  if (level_h.boundary_condition.type == BC_PERIODIC)
  {
    double average_value_of_f = mean (&level_h, VECTOR_F);
    if (average_value_of_f != 0.0)
    {
      if (ParallelDescriptor::IOProcessor())
      {
        std::cerr << "WARNING: Periodic boundary conditions, but f does not sum to zero... mean(f)=" << average_value_of_f << std::endl;
      }
      //shift_vector(&level_h,VECTOR_F,VECTOR_F,-average_value_of_f);
    }
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rebuild_operator(&level_h,NULL,a,b);    // i.e. calculate Dinv and lambda_max
  MGBuild(&MG_h,&level_h,a,b,minCoarseDim); // build the Multigrid Hierarchy
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (ParallelDescriptor::IOProcessor())
      std::cout << std::endl << std::endl << "===== STARTING SOLVE =====" << std::endl << std::flush;

  MGResetTimers (&MG_h);
  zero_vector (MG_h.levels[0], VECTOR_U);
#ifdef USE_FCYCLES
  FMGSolve (&MG_h, 0, VECTOR_U, VECTOR_F, a, b, tolerance_abs, tolerance_rel);
#else
  MGSolve (&MG_h, 0, VECTOR_U, VECTOR_F, a, b, tolerance_abs, tolerance_rel);
#endif /* USE_FCYCLES */

  MGPrintTiming (&MG_h, 0);   // don't include the error check in the timing results
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (ParallelDescriptor::IOProcessor())
    std::cout << std::endl << std::endl << "===== Performing Richardson error analysis ==========================" << std::endl;
  // solve A^h u^h = f^h
  // solve A^2h u^2h = f^2h
  // solve A^4h u^4h = f^4h
  // error analysis...
  MGResetTimers(&MG_h);
  const double dtol = tolerance_abs;
  const double rtol = tolerance_rel;
  int l;for(l=0;l<3;l++){
    if(l>0)restriction(MG_h.levels[l],VECTOR_F,MG_h.levels[l-1],VECTOR_F,RESTRICT_CELL);
           zero_vector(MG_h.levels[l],VECTOR_U);
    #ifdef USE_FCYCLES
    FMGSolve(&MG_h,l,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #else
     MGSolve(&MG_h,l,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
    #endif
  }
  richardson_error(&MG_h,0,VECTOR_U);

  // Now convert solution from HPGMG back to rhs MultiFab.
  for (MFIter mfi(rhs); mfi.isValid(); ++mfi)
  {
      const Box &bx = mfi.validbox();
      double *fab_data = rhs[mfi].dataPtr();

      // First find the HPGMG box corresponding to this BoxLib box.
      const int *loVect = bx.loVect();
      int box;
      for (box = 0; box < level_h.num_my_boxes; ++box)
      {
        if ((level_h.my_boxes[box].low.i == loVect[0]) &&
            (level_h.my_boxes[box].low.j == loVect[1]) &&
            (level_h.my_boxes[box].low.k == loVect[2]))
          break;
      }

      const Box &fabbox = mfi.fabbox();

      // Found the matching boxes, now fill the data.
      const int dim_i = level_h.my_boxes[box].dim;
      const int dim_j = level_h.my_boxes[box].dim;
      const int dim_k = level_h.my_boxes[box].dim;
      const int ghosts = level_h.my_boxes[box].ghosts;
      const int jStride = level_h.my_boxes[box].jStride;
      const int kStride = level_h.my_boxes[box].kStride;
      const int BoxLib_ghosts = rhs.nGrow();

      for(int k=0;k<dim_k;k++){
      for(int j=0;j<dim_j;j++){
      for(int i=0;i<dim_i;i++){

        const int ijk_HPGMG = (i+ghosts) + (j+ghosts)*jStride + (k+ghosts)*kStride;

	// WARNING: this indexing stride works for FABs *ONLY* if we have ONE
	// component in the FAB. If we have more than one we have to stride
	// over the components in the outermost loop (outside of k).
        const int BL_jStride = fabbox.length(0);
        const int BL_kStride = fabbox.length(0) * fabbox.length(1);
        const int ijk_BoxLib = (i+BoxLib_ghosts) + (j+BoxLib_ghosts)*BL_jStride + (k+BoxLib_ghosts)*BL_kStride;

        fab_data[ijk_BoxLib] = level_h.my_boxes[box].vectors[VECTOR_U][ijk_HPGMG];
      }}}
  }
  rhs.FillBoundary();
  geom.FillPeriodicBoundary(rhs);

  const double norm_from_HPGMG = norm(&level_h, VECTOR_U);
  const double mean_from_HPGMG = mean(&level_h, VECTOR_U);
  const Real norm0 = rhs.norm0();
  const Real norm2 = rhs.norm2();
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "mean from HPGMG: " << mean_from_HPGMG << std::endl;
    std::cout << "norm from HPGMG: " << norm_from_HPGMG << std::endl;
    std::cout << "norm0 of RHS copied to MF: " << norm0 << std::endl;
    std::cout << "norm2 of RHS copied to MF: " << norm2 << std::endl;
  }

  // Write the MF to disk for comparison with the in-house solver
  writePlotFile("MF_HPGMG", rhs, geom);

  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  MGDestroy(&MG_h);
  destroy_level(&level_h);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#else // use in-house multigrid solver

  MultiFab anaSoln;
  if (comp_norm || plot_err || plot_asol) {
    anaSoln.define(bs, Ncomp, 0, Fab_allocate);
    compute_analyticSolution(anaSoln);

    if (plot_asol) {
      writePlotFile("ASOL", anaSoln, geom);
    }
  }

  // Allocate the solution array
  // Set the number of ghost cells in the solution array.
  MultiFab soln(bs, Ncomp, 1, Fab_allocate);

#ifdef USEHYPRE
  if (solver_type == Hypre || solver_type == All) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Solving with Hypre " << std::endl;
    }

    solve(soln, anaSoln, a, b, alpha, beta, rhs, bs, geom, Hypre);
  }
#endif /* USEHYPRE */

  if (solver_type == BoxLib_C || solver_type == All) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Solving with BoxLib C++ solver " << std::endl;
    }

    solve(soln, anaSoln, a, b, alpha, beta, rhs, bs, geom, BoxLib_C);
  }

#ifdef USE_F90_SOLVERS
  if (solver_type == BoxLib_F || solver_type == All) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Solving with BoxLib F90 solver " << std::endl;
    }

    solve(soln, anaSoln, a, b, alpha, beta, rhs, bs, geom, BoxLib_F);
  }
#endif /* USE_F90_SOLVERS */

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "----------------------------------------" << std::endl;
  }

  const Real norm0 = soln.norm0();
  const Real norm1 = soln.norm1();
  const Real norm2 = soln.norm2();
  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "norm0 of RHS with in-house multigrid solver: " << norm0 << std::endl;
    std::cout << "norm1 of RHS with in-house multigrid solver: " << norm1 << std::endl;
    std::cout << "norm2 of RHS with in-house multigrid solver: " << norm2 << std::endl;
  }

  // Write the MF to disk for comparison with the HPGMG solver
  writePlotFile("MF_BoxLib", soln, geom);

#endif /* USEHPGMG */
  BL_PROFILE_VAR_STOP(pmain);
  BoxLib::Finalize();
}

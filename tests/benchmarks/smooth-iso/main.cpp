#include "hermes2d.h"
#include "solver_umfpack.h"

//  This test makes sure that the benchmark "smooth" works correctly.

int P_INIT = 1;                   // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;     // This is a quantitative parameter of the adapt(...) function and
                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;           // Adaptive strategy:
                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                  //   error is processed. If more elements have similar errors, refine
                                  //   all to keep the mesh symmetric.
                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                  //   than THRESHOLD times maximum element error.
                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                  //   than THRESHOLD.
                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const RefinementSelectors::AllowedCandidates ADAPT_TYPE = RefinementSelectors::H2DRS_CAND_HP;         // Type of automatic adaptivity.
const bool ISO_ONLY = false;      // Isotropic refinement flag (concerns quadrilateral elements only).
                                  // ISO_ONLY = false ... anisotropic refinement of quad elements
                                  // is allowed (default),
                                  // ISO_ONLY = true ... only isotropic refinements of quad elements
                                  // are allowed.
const int MESH_REGULARITY = -1;   // Maximum allowed level of hanging nodes:
                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                  // Note that regular meshes are not supported, this is due to
                                  // their notoriously bad performance.
const double ERR_STOP = 1e-4;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 400;       // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// exact solution
static double fn(double x, double y)
{
  return sin(x)*sin(y);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  dx = cos(x)*sin(y);
  dy = sin(x)*cos(y);
  return fn(x, y);
}

// boundary condition types
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// function values for Dirichlet boundary conditions
scalar bc_values(int marker, double x, double y)
{
  return fn(x, y);
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  return 2 * sin(x) * sin(y);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);
  if(P_INIT == 1) P_INIT++;  // this is because there are no degrees of freedom
                             // on the coarse mesh lshape.mesh if P_INIT == 1

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form), SYM);
  wf.add_liform(0, callback(linear_form));

  // visualize solution and mesh
  //ScalarView sview("Coarse solution", 0, 0, 500, 400);
  //OrderView  oview("Polynomial orders", 505, 0, 500, 400);

  // matrix solver
  UmfpackSolver solver;

  // prepare selector
  RefinementSelectors::H1NonUniformHP selector(ISO_ONLY, ADAPT_TYPE, 1.0, H2DRS_DEFAULT_ORDER, &shapeset);

  // DOF and CPU convergence graphs
  SimpleGraph graph_dof_est, graph_dof_exact, graph_cpu_est, graph_cpu_exact;

  // adaptivity loop
  int it = 1, ndofs;
  bool done = false;
  double cpu = 0.0;
  Solution sln_coarse, sln_fine;
  do
  {
    info("\n---- Adaptivity step %d ---------------------------------------------\n", it++);

    // time measurement
    begin_time();

    // solve the coarse mesh problem
    LinSystem ls(&wf, &solver);
    ls.set_spaces(1, &space);
    ls.set_pss(1, &pss);
    ls.assemble();
    ls.solve(1, &sln_coarse);

    // time measurement
    cpu += end_time();

    // calculate error wrt. exact solution
    ExactSolution exact(&mesh, fndd);
    double error = h1_error(&sln_coarse, &exact) * 100;
    info("\nExact solution error: %g%%", error);

    // view the solution
    //sview.show(&sln_coarse);
    //oview.show(&space);

    // time measurement
    begin_time();

    // solve the fine mesh problem
    RefSystem rs(&ls);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate error estimate wrt. fine mesh solution
    H1AdaptHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Estimate of error: %g%%", err_est);

    // add entries to DOF convergence graphs
    graph_dof_exact.add_values(space.get_num_dofs(), error);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof_est.dat");

    // add entries to CPU convergence graphs
    graph_cpu_exact.add_values(cpu, error);
    graph_cpu_exact.save("conv_cpu_exact.dat");
    graph_cpu_est.add_values(cpu, err_est);
    graph_cpu_est.save("conv_cpu_est.dat");

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      done = hp.adapt(THRESHOLD, STRATEGY, &selector, MESH_REGULARITY);
      ndofs = space.assign_dofs();
      if (ndofs >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu += end_time();

    // wait for keyboard or mouse input
    //sview.wait_for_keypress("Click into the mesh window and press any key to proceed.");
  }
  while (done == false);
  verbose("Total running time: %g sec", cpu);

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 49;
  printf("n_dof_actual = %d\n", ndofs);
  printf("n_dof_allowed = %d\n", n_dof_allowed); // ndofs was 49 at the time this test was created
  if (ndofs <= n_dof_allowed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}



#include "hermes2d.h"
#include "solver_umfpack.h"

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear time-dependent PDE discretized implicitly in time
//  using implicit Euler or Crank-Nicolson.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt - div[lambda(u)grad u] = f
//
//  Domain: square (-10,10)^2
//
//  BC:  Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int P_INIT = 2;             // Initial polynomial degree of all mesh elements.
const int PROJ_TYPE = 1;          // 1 for H1 projections, 0 for L2 projections
const int TIME_DISCR = 2;         // 1 for implicit Euler, 2 for Crank-Nicolson
const double TAU = 0.5;           // Time step
const double T_FINAL = 5.0;       // Time interval length
const int INIT_GLOB_REF_NUM = 2;  // Number of initial uniform mesh refinements

// Adaptivity
const int UNREF_FREQ = 1;         // Every UNREF_FREQth time step the mesh is unrefined.
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
const int ADAPT_TYPE = 0;         // Type of automatic adaptivity:
                                  // ADAPT_TYPE = 0 ... adaptive hp-FEM (default),
                                  // ADAPT_TYPE = 1 ... adaptive h-FEM,
                                  // ADAPT_TYPE = 2 ... adaptive p-FEM.
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
const double ERR_STOP = 1.0;      // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations
const bool NEWTON_ON_COARSE_MESH = false;  // true... Newton is done on coarse mesh in every adaptivity step
                                           // false...Newton is done on coarse mesh only once, then projection
                                           // of the fine mesh solution to coarse mesh is used

// Thermal conductivity (temperature-dependent)
// Note: for any u, this function has to be positive
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition
scalar initial_condition(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary condition type
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Dirichlet boundary condition values
scalar bc_values(int marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y')
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // load the mesh and perform initial refinements
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create finite element space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // enumerate basis functions
  space.assign_dofs();

  // Solutions for the time stepping and the Newton's method
  Solution u_prev_time, u_prev_newton;

  // initialize the weak formulation
  WeakForm wf(1);
  if(TIME_DISCR == 1) {
    wf.add_biform(0, 0, callback(J_euler), UNSYM, ANY, 1, &u_prev_newton);
    wf.add_liform(0, callback(F_euler), ANY, 2, &u_prev_newton, &u_prev_time);
  }
  else {
    wf.add_biform(0, 0, callback(J_cranic), UNSYM, ANY, 1, &u_prev_newton);
    wf.add_liform(0, callback(F_cranic), ANY, 2, &u_prev_newton, &u_prev_time);
  }

  // initialize the nonlinear system and solver
  UmfpackSolver solver;
  NonlinSystem nls(&wf, &solver);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // visualize solution and mesh
  ScalarView view("Initial condition", 0, 0, 700, 600);
  view.fix_scale_width(80);
  OrderView ordview("Initial mesh", 700, 0, 700, 600);

  // error estimate and discrete problem size as a function of physical time
  SimpleGraph graph_time_err, graph_time_dof;

  // project the function initial_condition() on the coarse mesh
  nls.set_ic(initial_condition, &mesh, &u_prev_time, PROJ_TYPE);
  u_prev_newton.copy(&u_prev_time);

  // view initial guess for Newton's method
  // satisfies BC conditions
  view.show(&u_prev_newton);
  ordview.show(&space);

  // Newton's loop on the coarse mesh
  info("---- Time step 1, Newton solve on coarse mesh:\n");
  if (!nls.solve_newton_1(&u_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
    error("Newton's method did not converge.");

  // store the result in sln_coarse
  Solution sln_coarse, sln_fine;
  sln_coarse.copy(&u_prev_newton);

  // time stepping loop
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int n = 1; n <= nstep; n++)
  {
    // periodic global derefinements
    if (n > 1 && n % UNREF_FREQ == 0) {
      info("---- Time step %d, global derefinement.\n", n);
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
      space.assign_dofs();

      // project the fine mesh solution on the globally derefined mesh
      info("---- Time step %d, projecting fine mesh solution on globally derefined mesh:\n", n);
      nls.set_ic(&sln_fine, &u_prev_newton, PROJ_TYPE);

      if (NEWTON_ON_COARSE_MESH) {
        // Newton's loop on the globally derefined mesh
        info("---- Time step %d, Newton solve on globally derefined mesh:\n", n);
        if (!nls.solve_newton_1(&u_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
          error("Newton's method did not converge.");
      }

      // store the result in sln_coarse
      sln_coarse.copy(&u_prev_newton);
    }

    // adaptivity loop
    bool done = false;
    double err_est;
    int a_step = 1;
    do
    {
      info("---- Time step %d, adaptivity step %d, Newton solve on fine mesh:\n", n, a_step);

      // reference nonlinear system
      RefNonlinSystem rnls(&nls);
      rnls.prepare();

      // set initial condition for the Newton's method on the fine mesh
      if (a_step == 1) rnls.set_ic(&sln_coarse, &u_prev_newton);
      else rnls.set_ic(&sln_fine, &u_prev_newton);

      // Newton's method on fine mesh
      if (!rnls.solve_newton_1(&u_prev_newton, NEWTON_TOL_FINE, NEWTON_MAX_ITER))
        error("Newton's method did not converge.");

      // store the result in sln_fine
      sln_fine.copy(&u_prev_newton);

      // visualization of intermediate solution
      // and mesh during adaptivity
      char title[100];
      sprintf(title, "Solution, time level %d, adapt step %d", n, a_step);
      view.set_title(title);
      view.show(&sln_fine);
      sprintf(title, "Fine mesh, time level %d, adapt step %d", n, a_step);
      ordview.set_title(title);
      ordview.show(rnls.get_space(0));

      // calculate error estimate wrt. fine mesh solution
      H1OrthoHP hp(1, &space);
      err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
      info("Error estimate: %g%", err_est);

      // if err_est too large, adapt the mesh
      if (err_est < ERR_STOP) done = true;
      else {
        hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
        int ndof = space.assign_dofs();
        if (ndof >= NDOF_STOP) done = true;

        // project the fine mesh solution on the new coarse mesh
        info("---- Time step %d, adaptivity step %d, projecting fine mesh solution on new coarse mesh:\n",
          n, a_step);
        nls.set_ic(&sln_fine, &u_prev_newton, PROJ_TYPE);

        // moving to new adaptivity step
        a_step++;

        if (NEWTON_ON_COARSE_MESH) {
          // Newton's loop on the coarse mesh
          info("---- Time step %d, adaptivity step %d, Newton solve on new coarse mesh:\n", n, a_step);
          if (!nls.solve_newton_1(&u_prev_newton, NEWTON_TOL_COARSE, NEWTON_MAX_ITER))
            error("Newton's method did not converge.");
        }

        // store the result in sln_coarse
        sln_coarse.copy(&u_prev_newton);
      }
    }
    while (!done);

    // add entries to convergence graphs
    graph_time_err.add_values(n*TAU, err_est);
    graph_time_err.save("time_error.dat");
    graph_time_dof.add_values(n*TAU, space.get_num_dofs());
    graph_time_dof.save("time_dof.dat");

    // copying the new time level solution into u_prev_time
    u_prev_time.copy(&sln_fine);
  }

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}

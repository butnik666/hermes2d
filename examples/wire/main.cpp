#include "hermes2d.h"
#include "solver_umfpack.h"  // defines the class UmfpackSolver

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -Laplace A + ii*omega*gamma*mu*A = mu *J_ext
//
//  Domain: Rectangle of height 0.003 and width 0.004. Different
//  materials for the wire, air, and iron (see mesh file domain2.mesh).
//
//  BC: Zero Dirichlet on the top and right edges, zero Neumann
//  elsewhere.
//
//  The following parameters can be changed:

const int P_INIT = 1;             // Initial polynomial degree of all mesh elements.
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
const double ERR_STOP = 0.01;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 50000;      // Adaptivity process stops when the number of degrees of freedom grows
                                  // over this limit. This is to prevent h-adaptivity to go on forever.


// Problem parameters
double mu_0 = 4.0*3.141592654E-7;
double J_wire = 5000000.0;
double freq = 5E3;
double omega = 2*3.141592654*freq;
double gamma_iron = 6E6;
double mu_iron = 1000*mu_0;

int bc_types(int marker)
{
  if (marker==1) {return BC_NATURAL;}
  if (marker==2) {return BC_ESSENTIAL;}
  if (marker==3) {return BC_ESSENTIAL;}
  if (marker==4) {return BC_ESSENTIAL;}
}

cplx dir_bc_values(int marker, double x, double y)
{
  return cplx(0.0,0.0);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_iron(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  scalar ii = cplx(0.0, 1.0);
  return 1./mu_iron * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + ii*omega*gamma_iron*int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_wire(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 1./mu_0 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_air(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return 1./mu_0 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v); // conductivity gamma is zero
}

template<typename Real, typename Scalar>
Scalar linear_form_wire(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return J_wire * int_v<Real, Scalar>(n, wt, v);
}

int main(int argc, char* argv[])
{
  // load the mesh
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain2.mesh", &mesh);

  // initial uniform subdivision
  //mesh.refine_all_elements();

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(dir_bc_values);
  space.set_uniform_order(P_INIT);

  // enumerate basis functions
  space.assign_dofs();

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form_iron), SYM, 3);
  wf.add_biform(0, 0, callback(bilinear_form_wire), SYM, 2);
  wf.add_biform(0, 0, callback(bilinear_form_air), SYM, 1);
  wf.add_liform(0, callback(linear_form_wire), 2);

  // visualize solution and mesh
  ScalarView view("Vector potential A", 0, 0, 1000, 600);
  OrderView  oview("Polynomial orders", 1100, 0, 900, 600);

  // matrix solver
  UmfpackSolver solver;

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
    LinSystem sys(&wf, &solver);
    sys.set_spaces(1, &space);
    sys.set_pss(1, &pss);

    // time measurement
    cpu += end_time();

    // assemble the stiffness matrix and solve the system
    sys.assemble();
    sys.solve(1, &sln_coarse);

    // visualize the solution
    view.show(&sln_coarse, EPS_HIGH);
    oview.show(&space);

    // time measurement
    begin_time();

    // solve the fine mesh problem
    RefSystem rs(&sys);
    rs.assemble();
    rs.solve(1, &sln_fine);

    // calculate element errors and total error estimate
    H1OrthoHP hp(1, &space);
    double err_est = hp.calc_error(&sln_coarse, &sln_fine) * 100;
    info("Error estimate: %g%%", err_est);

    // add entries to DOF convergence graph
    graph_dof_est.add_values(space.get_num_dofs(), err_est);
    graph_dof_est.save("conv_dof.dat");

    // add entries to CPU convergence graph
    graph_cpu_est.add_values(cpu, err_est);
    graph_cpu_est.save("conv_cpu.dat");

    // if err_est too large, adapt the mesh
    if (err_est < ERR_STOP) done = true;
    else {
      hp.adapt(THRESHOLD, STRATEGY, ADAPT_TYPE, ISO_ONLY, MESH_REGULARITY);
      ndofs = space.assign_dofs();
      if (ndofs >= NDOF_STOP) done = true;
    }

    // time measurement
    cpu += end_time();
  }
  while (done == false);
  verbose("Total running time: %g sec", cpu);

  // show the fine solution - this is the final result
  view.set_title("Final solution");
  view.show(&sln_fine);

  // wait for keyboard or mouse input
  View::wait("Waiting for all views to be closed.");
  return 0;
}

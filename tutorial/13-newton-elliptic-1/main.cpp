#include "hermes2d.h"
#include "solver_umfpack.h"
#include "function.h"

//  This example shows an introductory application of the Newton's
//  method to a nonlinear elliptic problem. We use zero Dirichlet boundary
//  conditions and a constant initial guess for the Newton's method.
//  The treatment of nonzero Dirichlet BC and a more general initial guess
//  will be shown in the next example newton-elliptic-2.
//
//  PDE: stationary heat transfer equation with nonlinear thermal
//  conductivity, - div[lambda(u)grad u] = 0
//
//  Domain: unit square (-10,10)^2
//
//  BC: Zero Dirichlet
//
//  The following parameters can be changed:

const int P_INIT = 2;             // Initial polynomial degree
const int PROJ_TYPE = 1;          // For the projection of the initial condition
                                  // on the initial mesh: 1 = H1 projection,
                                  // 0 = L2 projection
const double NEWTON_TOL = 1e-6;   // Stopping criterion for the Newton's method
const int NEWTON_MAX_ITER = 100;  // Maximum allowed number of Newton iterations
const int INIT_GLOB_REF_NUM = 3;  // Number of initial uniform mesh refinements
const int INIT_BDY_REF_NUM = 5;   // Number of initial refinements towards boundary

// Thermal conductivity (temperature-dependent)
// Note: for any u, this function has to be positive
template<typename Real>
Real lam(Real u) { return 1 + pow(u, 4); }

// Derivative of the thermal conductivity with respect to 'u'
template<typename Real>
Real dlam_du(Real u) { return 4*pow(u, 3); }

// Boundary condition type (essential = Dirichlet)
int bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Heat sources (can be a general function of 'x' and 'y')
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Jacobian matrix
template<typename Real, typename Scalar>
Scalar jac(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (dlam_du(u_prev->val[i]) * u->val[i] * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
                       + lam(u_prev->val[i]) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));

  return result;
}

// Fesidual vector
template<typename Real, typename Scalar>
Scalar res(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Func<Scalar>* u_prev = ext->fn[0];
  for (int i = 0; i < n; i++)
    result += wt[i] * (lam(u_prev->val[i]) * (u_prev->dx[i] * v->dx[i] + u_prev->dy[i] * v->dy[i])
		       - heat_src(e->x[i], e->y[i]) * v->val[i]);
  return result;
}

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // initial mesh refinements
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1,INIT_BDY_REF_NUM);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_uniform_order(P_INIT);
  space.assign_dofs();

  // previous solution for the Newton's iteration
  Solution u_prev;

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(jac), UNSYM, ANY, 1, &u_prev);
  wf.add_liform(0, callback(res), ANY, 1, &u_prev);

  // initialize the nonlinear system and solver
  UmfpackSolver umfpack;
  NonlinSystem nls(&wf, &umfpack);
  nls.set_spaces(1, &space);
  nls.set_pss(1, &pss);

  // use a constant function as the initial guess
  u_prev.set_const(&mesh, 3.0);
  nls.set_ic(&u_prev, &u_prev, PROJ_TYPE);

  // Newton's loop
  if (!nls.solve_newton_1(&u_prev, NEWTON_TOL, NEWTON_MAX_ITER)) error("Newton's method did not converge.");

  // visualise the solution and mesh
  ScalarView sview("Solution", 0, 0, 800, 600);
  OrderView oview("Mesh", 820, 0, 800, 600);
  char title[100];
  sprintf(title, "Solution");
  sview.set_title(title);
  sview.show(&u_prev);
  sprintf(title, "Mesh");
  oview.set_title(title);
  oview.show(&space);

  // wait for keyboard or mouse input
  View::wait();
  return 0;
}


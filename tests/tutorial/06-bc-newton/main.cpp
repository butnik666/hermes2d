#include "hermes2d.h"
#include "solver_umfpack.h"

// This test makes sure that example 06-bc-newton works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

double T1 = 30.0;            // prescribed temperature on Gamma_3
double T0 = 20.0;            // outer temperature on Gamma_1
double H  = 0.05;            // heat flux on Gamma_1
int UNIFORM_REF_LEVEL = 1;   // number of initial uniform mesh refinements
int CORNER_REF_LEVEL = 3;   // number of mesh refinements towards the re-entrant corner

// boundary condition types
int bc_types(int marker)
  { return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL; }

// function values for Dirichlet boundary markers
scalar bc_values(int marker, double x, double y)
  { return T1; }

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return H * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return T0 * H * int_v<Real, Scalar>(n, wt, v);
}


int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);
  for(int i=0; i<UNIFORM_REF_LEVEL; i++) mesh.refine_all_elements();
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // initialize the shapeset and the cache
  H1Shapeset shapeset;
  PrecalcShapeset pss(&shapeset);

  // create an H1 space
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_bc_values(bc_values);

  // initialize the weak formulation
  WeakForm wf(1);
  wf.add_biform(0, 0, callback(bilinear_form));
  wf.add_biform_surf(0, 0, callback(bilinear_form_surf), 1);
  wf.add_liform_surf(0, callback(linear_form_surf), 1);

  // initialize the linear system and solver
  UmfpackSolver umfpack;
  LinSystem sys(&wf, &umfpack);
  sys.set_spaces(1, &space);
  sys.set_pss(1, &pss);

  // testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  int success = 1;
  for (int p_init = 1; p_init <= 10; p_init++) {
    printf("********* p_init = %d *********\n", p_init);
    space.set_uniform_order(p_init);
    space.assign_dofs();

    // assemble the stiffness matrix and solve the system
    Solution sln;
    sys.assemble();
    sys.solve(1, &sln);

    scalar *sol_vector;
    int n_dof;
    sys.get_solution_vector(sol_vector, n_dof);
    printf("n_dof = %d\n", n_dof);
    double sum = 0;
    for (int i=0; i < n_dof; i++) sum += sol_vector[i];
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the
    // current shapeset. If you change the shapeset,
    // you need to correct these numbers.
    if (p_init == 1 && fabs(sum - 1146.15) > 1e-1) success = 0;
    if (p_init == 2 && fabs(sum - 1145.97) > 1e-1) success = 0;
    if (p_init == 3 && fabs(sum - 1145.97) > 1e-1) success = 0;
    if (p_init == 4 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 5 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 6 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 7 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 8 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 9 && fabs(sum - 1145.96) > 1e-1) success = 0;
    if (p_init == 10 && fabs(sum - 1145.96) > 1e-1) success = 0;
  }

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

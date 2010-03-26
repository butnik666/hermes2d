#include "hermes2d.h"
#include "solver_umfpack.h"
#include <iostream>
using namespace std;


//example
// The following parameters can be changed:

const int P_INIT = 1;

// projected function
double F(double x, double y)
{
  return x*x*x + y*y*y;
}

// bilinear and linear form defining the projection
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, u, v);
}

// return the value \int v dx
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * ((pow(e->x[i], 3) + pow(e->y[i], 3)) * v->val[i]);
  return result;
}

// boundary conditions
int bc_types(int marker)
{
   return BC_NONE;
}

int main(int argc, char* argv[])
{
	// load the mesh file
	Mesh mesh;
	H2DReader mloader;
	mloader.load("domain.mesh", &mesh);


	// perform some sample initial refinements//do nothing,  means that there is only one dividing
	mesh.refine_all_elements(2); // Refines all elements.
  mesh.refine_towards_vertex(2, 1);    // Refines mesh towards vertex #3 (4x).
  mesh.refine_towards_boundary(2, 1);  // Refines all elements along boundary 2 (4x).
  mesh.refine_all_elements();          // Refines element #86 isotropically.
  mesh.refine_element(18, 0);         // Refines element #18 isotropically.
//  mesh.refine_all_elements();
  // initialize the shapeset and the cache
   L2Shapeset shapeset;
   PrecalcShapeset pss(&shapeset);

   // create the L2 space
   L2Space space(&mesh, &shapeset);
   space.set_bc_types(bc_types);

   // set uniform polynomial degrees
   space.set_uniform_order(P_INIT);

   // enumerate basis functions
   space.assign_dofs();

   // display the mesh
    MeshView mview("Hello world!", 100, 100, 500, 500);  // (100, 100) is the upper left corner position
    mview.show(&mesh);                                   // 500 x 500 is the window size


   Solution sln;

   // matrix solver
   UmfpackSolver umfpack;

   // initialize the weak formulation
   WeakForm wf(1);
   wf.add_biform(0, 0, callback(bilinear_form));
   wf.add_liform(0, callback(linear_form));

   // assemble and solve the finite element problem
   LinSystem sys(&wf, &umfpack);
   sys.set_spaces(1, &space);
   sys.set_pss(1, &pss);
   sys.assemble();
   sys.solve(1, &sln);



  //my code
   Element* e;
   Neighbor* neighb;
   for_all_active_elements(e, &mesh){
     neighb = new Neighbor(e, &sln);
     for (int i = 0;  i < e->nvert; i++)
     {
    	 if(e->en[i]->bnd == 0)
    	 {
    		 neighb->set_active_edge(i);
    	 }
     }
     vector<int>* sousedi;
     sousedi = neighb->get_neighbors();
     cout <<"\n"<< sousedi->size();
  	 delete neighb;
    }

   // visualize the solution
//	  ScalarView view1("Solution 1");
//	  view1.show(&sln);
//	  view1.wait_for_keypress();

	  // wait for keyboard or mouse input
	  View::wait("Waiting for keyboard or mouse input.");

	return 0;
}


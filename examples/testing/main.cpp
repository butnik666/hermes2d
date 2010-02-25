#include "hermes2d.h"
#include "solver_umfpack.h"
#include <iostream>
using namespace std;

//way up
void finding_act_elem(Mesh* mesh, Element* elem, int edge_num)
{
	Node* edge = NULL;
	int p1, p2; //id of parents of the edge

	p1 = elem->vn[edge_num]->id;
	p2 = elem->vn[(edge_num + 1) % elem->nvert]->id;
	cout << p1 << " " << p2 << "\n";
	edge = mesh->peek_edge_node(p1, p2);

	if (edge == NULL)
	{
		finding_act_elem(mesh, elem->parent, edge_num);
	} else
		for (int i = 0; i < 2; i++)
		{
			if ((edge->elem[i] != NULL) && (edge->elem[i]->active == 1))
				//do something
				cout << "way up neighb:" << edge->elem[i]->id << "\n";
		}
}

//way down
void finding_act_elem(Mesh* mesh, Node* vertex)
{
	int son;
	int parents[2];

	Node* edge = NULL;
	Node* n = NULL;
	son = vertex->id;

	parents[0] = vertex->p1;
	parents[1] = vertex->p2;

	for (int i = 0; i < 2; i++)
	{
		edge = mesh->peek_edge_node(son, parents[i]);
		//test if edge is active, means on one of sides can!! be active element
		if (edge == NULL)
		{
			n = mesh->peek_vertex_node(son, parents[i]);
			if(n == NULL)
				error("wasn't able to find middle vertex");
			else
				finding_act_elem(mesh, n);
		} else
			//test if on one of sides is active element
			for (int j = 0; j < 2; j++)
			{
				if (edge->elem[j] != NULL)
					if (edge->elem[j]->active == 1)
						//do something
						cout << "id: " << edge->elem[j]->id << "\n";
			}
	}
}

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

	// perform some sample initial refinements
	mesh.refine_all_elements(2); // Refines all elements.
  mesh.refine_towards_vertex(2, 1);    // Refines mesh towards vertex #3 (4x).
  mesh.refine_towards_boundary(2, 1);  // Refines all elements along boundary 2 (4x).
  mesh.refine_all_elements();          // Refines element #86 isotropically.
  mesh.refine_element(30, 0);         // Refines element #112 isotropically.
  mesh.refine_all_elements();
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

   //testing of transformations

/*   Element* center = mesh.get_element(50);
   sln.set_active_element(center);
   sln.push_transform(3);
   Quad2D* quad = sln.get_quad_2d();
   int eo = quad->get_edge_points(1);
   sln.set_quad_order(eo);
   double *fn = sln.get_fn_values();

   // print the values
   printf("%d: ", 1);
   int np = quad->get_num_points(eo);
   for (int i = 0; i < np; i++)
   printf(" % lf", fn[i]);
   printf("\n");
*/

  //my code
	Element* e = NULL;
	Element* neighb;
	for_all_active_elements(e, &mesh)
			{
		printf("vertices: %d %d %d %d\n", e->vn[0]->id,	e->vn[1]->id, e->vn[2]->id, e->vn[3]->id);
		for (int i = 0; i < e->nvert; i++)
				{
					if (e->en[i]->bnd == 0)
					{
//						printf("element id:%d , edge: %d, vertices: %d %d \n", e->id,	e->en[i]->id, e->en[i]->p1, e->en[i]->p2);

						neighb == NULL;
						neighb = e->get_neighbor(i);
						if (neighb != NULL)
						{
							for (int j = 0; j < neighb->nvert; j++)
							{
								if (e->en[i] == neighb->en[j])
									cout << "aktivni elem" << neighb->id << "\n";
							}
						} else
						{
							Node* vertex = NULL;
							Element* parent = NULL;
							vertex = mesh.peek_vertex_node(e->en[i]->p1, e->en[i]->p2);

							if (vertex == NULL)
							{
								//way up
								int edge_num = i;
								parent = e->parent;
								finding_act_elem(&mesh, parent, edge_num);
							} else
								//way down
								finding_act_elem(&mesh, vertex);
						}
					}
				}
			}

	 // visualize the solution
	  ScalarView view1("Solution 1");
	//  view1.show(&sln);
	//  view1.wait_for_keypress();

	  // wait for keyboard or mouse input
	  View::wait("Waiting for keyboard or mouse input.");
	return 0;
}


#include "hermes2d.h"
#include "solver_umfpack.h"
#include <iostream>
using namespace std;
/*
//way up
void finding_act_elem(Solution* sln, Mesh* mesh, Element* elem, int edge_num, int* orig_vertex_id, Node** road_vertices, int n_road_vertices)
{

	Node* edge = NULL;
	Node* vertex = NULL;
	int p1, p2; //id of parents of the edge

	p1 = elem->vn[edge_num]->id;
	p2 = elem->vn[(edge_num + 1) % elem->nvert]->id;
	
	edge = mesh->peek_edge_node(p1, p2);
	vertex = mesh->peek_vertex_node(p1, p2);
	road_vertices[n_road_vertices] = vertex;
	n_road_vertices = n_road_vertices++;

	if (edge == NULL)
	{
		finding_act_elem(sln, mesh, elem->parent, edge_num, orig_vertex_id, road_vertices, n_road_vertices);

	} else
		for (int i = 0; i < 2; i++)
		{
			if ((edge->elem[i] != NULL) && (edge->elem[i]->active == 1)){

				//getting to correct edge
	//			cout << "way up neighb:" << edge->elem[i]->id << "\n";
				Element* neighb = edge->elem[i];
				int neighb_edge = -1;
				for(int j = 0; j < neighb->nvert; j++)
					if(neighb->en[j] == edge)
						neighb_edge = j;
				if(neighb_edge == -1) error("edge wasn't found");
			
				int v1, v2;
				v1 = p1;
				v2 = p2;

				Node* n = NULL;
				//set active the neighbour
				sln->set_active_element(neighb);

				// go threw between elements and set correct transformation
				for(int j = n_road_vertices; j > 0; j-- ){
					if(road_vertices[j] == NULL)
						if(j > 0){
						//do nothing,  means that there is only one dividing
							continue;
						}
						else
							error("shouldn't be here neighb doesn't have dividing");
					else{
						n = mesh->peek_vertex_node(road_vertices[j]->id, v1);
						if(n == NULL){
							n = mesh->peek_vertex_node(road_vertices[j]->id, v2);
							sln->push_transform((neighb_edge + 1) % neighb->nvert);
							v1 = n->id;
						}
						else{
								if(n->id == road_vertices[j-1]->id){
									sln->push_transform((neighb_edge + 1) % neighb->nvert);
									v2 = n->id;
								}
								else{
									n = mesh->peek_vertex_node(road_vertices[j]->id, v2);
									sln->push_transform(neighb_edge);
									v1 = n->id;												
								}
						}
					}
				}

				cout << orig_vertex_id[0]<<" "<< orig_vertex_id[1]<<" "<< road_vertices[0]->id<<"\n";				
				// final transformation on active element
				int test = 0;
				if (orig_vertex_id[0] == road_vertices[0]->id)
					test = 1;

				if(test == 1){
					sln->push_transform(neighb_edge);
				}
				else{
					sln->push_transform((neighb_edge + 1) % neighb->nvert);
				}

				Quad2D *quad = sln->get_quad_2d();
				int eo = quad->get_edge_points(neighb_edge);
				sln->set_quad_order(eo);
				double *fn = sln->get_fn_values();
				int np = quad->get_num_points(eo);
/*				for (int i = 0; i < np; i++)
				printf(" % lf", fn[i]);
				printf("\n");
*/
				//reset transformations
/*				sln->reset_transform();
				
			}
		}
}
*/
//way down
/*
void finding_act_elem(Solution* sln, Mesh* mesh, Node* vertex, int* par_vertex_id, Neighbor* neighbs, int* road, int n_road, int use_edge, int n_vert, int active_order, int* n_neighbs)
{
	int son;
	int parents[2];

	Node* edge = NULL;
	Node* n = NULL;
	Element* neighb = NULL;
	son = vertex->id;

	parents[0] = par_vertex_id[0];
	parents[1] = par_vertex_id[1];

	for (int i = 0; i < 2; i++)
	{
		road[n_road] = (use_edge + i) % n_vert;

		edge = mesh->peek_edge_node(son, parents[i]);
		//test if edge is active, means on one of sides can!! be active element
		if (edge == NULL)
		{
			n = mesh->peek_vertex_node(son, parents[i]);
			if(n == NULL)
				error("wasn't able to find middle vertex");
			else{
				if(i == 0) par_vertex_id[1] = son;
				else par_vertex_id[0] = son;

				n_road = n_road++;
				finding_act_elem(sln, mesh, n, par_vertex_id, neighbs, road, n_road, use_edge, n_vert, active_order, n_neighbs);
			}
		} else
			//test if on one of sides is active element
			for (int j = 0; j < 2; j++)
			{
				if (edge->elem[j] != NULL)
					if (edge->elem[j]->active == 1){
						//do something
							cout << "way down neighb: " << edge->elem[j]->id << "\n";
							neighb = mesh->get_element(edge->elem[j]->id);
							sln->set_active_element(neighb);
							Quad2D *quad = sln->get_quad_2d();
							int neighb_edge = -1;
							for(int k = 0; k < neighb->nvert; k++)
								if(neighb->en[k] == edge)
									neighb_edge = k;
							if(neighb_edge == -1) error("edge wasn't found");
							int eo = quad->get_edge_points(neighb_edge);
							int np = quad->get_num_points(eo);

							sln->set_quad_order(eo);
							double *fn = sln->get_fn_values();
							for (int k = 0; k < np; k++)
							printf(" % lf", fn[k]);
							printf("\n");

		
							//filling neighbors

							for(int k = 0; k <= n_road; k++) neighbs[*n_neighbs].transformations[k] = road[k];
							int neighb_order = sln->get_fn_order();

							neighbs[*n_neighbs].max_order = std::max(active_order, neighb_order);
							//in future need to set correct np and fna_values acc. max_order
							//add correct direction
							neighbs[*n_neighbs].n_fn_values = np;
							for(int k = 0; k <= np; k++) neighbs[*n_neighbs].fn_values[k] = fn[k];
							*n_neighbs = (*n_neighbs)++;
					}
			}
	}
}
*/
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

	Quad2D* quad = sln.get_quad_2d();
  int max_order = quad->get_max_order();
  double max_integ_n_points = (double)quad->get_num_points(max_order);
  max_integ_n_points = sqrt(max_integ_n_points);

	Element* e = NULL;


	e = mesh.get_element(12);
	Neighbor neighb(e, &sln);
	neighb.set_active_edge(1);

	neighb.set_active_edge(2);
	 // visualize the solution
//	  ScalarView view1("Solution 1");
//	  view1.show(&sln);
//	  view1.wait_for_keypress();

	  // wait for keyboard or mouse input
	  View::wait("Waiting for keyboard or mouse input.");
	return 0;
}


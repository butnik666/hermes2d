#include "neighbor.h"

Neighbor::Neighbor(Element* e, Solution* sln){
	central_el = e;
	sol = sln;
	quad = sln->get_quad_2d();
	mesh = sln->get_mesh();
	sol->set_active_element(central_el);
	central_order = sol->get_fn_order();
	for(int i = 0;  i < max_n_trans; i++)
	{
		fn_values[i] = NULL;
		fn_values_neighbor[i] = NULL;
	}
	n_neighbors = 0;
	neighbors_id.reserve(20 * e->nvert);
	solution_flag = 1;
};

Neighbor::Neighbor(Element* e, Mesh* given_mesh)
{
	central_el = e;
	sol = NULL;
	quad = NULL;
	mesh = given_mesh;
	central_order = -1;
	for(int i = 0;  i < max_n_trans; i++)
	{
		fn_values[i] = NULL;
		fn_values_neighbor[i] = NULL;
	}
	n_neighbors = 0;
	neighbors_id.reserve(20 * e->nvert);
	solution_flag = 0;
};


Neighbor::~Neighbor()
{
	for(int i = 0; i < max_n_trans; i++){
		if(fn_values[i] != NULL)
			delete[] fn_values[i];
		if(fn_values_neighbor[i] != NULL)
			delete[] fn_values_neighbor[i];
	}
};




//How works the finding neighbors
/*! \brief We will call the active  element as central element.
 *  If we have irregular mesh, there can be three options in relation of central element and neighbor element at common edge.
 *  First, the neighbor is same "size" as central, so the edge has active elements on both sides. This option is tested by function get_neighbor().
 *  Second, the neighbor is "bigger" then central. Then we to go "way up".
 *  Third, the neighbor is "smaller", we have more neighbors against the edge. This solves "way down".
The choice between way up or way down is made by testing if we can find vertex in the middle of the edge. If we can
then we go way down.

Also at every way we fill function values of central and neighbor elements threw function set_fn_values(). Last step is
possible change of order of neighbor's function values to correspond function values of central element at same points.

We also need transform solution either on neighbor or central element to get points at correct part of the edge. We use method "push_transform"
and use only range [0-3]. These types of transformation are common for triangles and quads and choosing right transformation can
be derived from local numbers of edges.

For numbering and ordering of edges, vertices and sons of an element look into mesh.cpp
*/
void Neighbor::set_active_edge(int edge)
{
	//erase all data from previous edge or element
	clean_all();

	active_edge = edge;

	debug_log("central element: %d \n\n\n", central_el->id);
	if (central_el->en[active_edge]->bnd == 0)
	{
		neighb_el = central_el->get_neighbor(active_edge);
		// test if on the other side of the edge is active element
		if (neighb_el != NULL)
		{
			debug_log("active neighbor el: %d \n", neighb_el->id);
			for (int j = 0; j < neighb_el->nvert; j++)
			{
				if (central_el->en[active_edge] == neighb_el->en[j])
				{
					neighbor_edge = j;

					if(solution_flag == 1){
						// fill function values of central a neighbor element
						set_fn_values(H2D_NO_TRANSF);

						// set same direction as central element
						int p1 = central_el->vn[active_edge]->id;
						int p2 = central_el->vn[(active_edge + 1) % central_el->nvert]->id;
						set_correct_direction(p1, p2, 0);
					}

					// raise the number of neighbors
					n_neighbors = n_neighbors++;

					// add neighbor id to neighbors_id
					neighbors_id.push_back(neighb_el->id);
				}
			}
		} else
		{
			Node* vertex = NULL;
			vertex = mesh->peek_vertex_node(central_el->en[active_edge]->p1,	central_el->en[active_edge]->p2);
			int orig_vertex_id[2];
			orig_vertex_id[0] = central_el->vn[active_edge]->id;
			orig_vertex_id[1]	= central_el->vn[(active_edge + 1) % central_el->nvert]->id;
			if (vertex == NULL)
			{
				// way up

				Element* parent = NULL;
				parent = central_el->parent;

				Node** road_vertices;
				road_vertices = new Node*[max_n_trans];

				int n_road_vertices = 0; // number of used vertices

				for (int j = 0; j < max_n_trans; j++)
					road_vertices[j] = NULL;

				finding_act_elem(parent, active_edge, orig_vertex_id, road_vertices, n_road_vertices);

				delete[] road_vertices;
			} else
			{
				// way down

				int road[max_n_trans]; //array for temporal transformation
				int n_road = 0; //number of used transformations

				finding_act_elem( vertex, orig_vertex_id, road, n_road,	active_edge, central_el->nvert);

				debug_log("number of neighbors: %d \n", n_neighbors);
			}
		}
	}
	else
		error("The given edge isn't inner");
};



//way up
/*! \brief Function for finding "bigger" neighbor
*If the neighbor is "bigger" then this means central element is descendant of some inactive elements. We go threw this parents and
*stop when against an edge, which has same local number as the original edge, we have active element.
*Important is that all sons have same orientation as parent, so local number of the edge is same.

*/
void Neighbor::finding_act_elem( Element* elem, int edge_num, int* orig_vertex_id, Node** road_vertices, int n_road_vertices)
{
	Node* edge = NULL;
	Node* vertex = NULL;
	int p1, p2; //id of parents of the edge

	// order parents in direction of parent element (needed for transformation of solution)
	p1 = elem->vn[edge_num]->id;
	p2 = elem->vn[(edge_num + 1) % elem->nvert]->id;
	
	int id_of_par_orient_1 = p1;
	int id_of_par_orient_2 = p2;

	// find if between parents p1 and p2 is active edge (is used by neighbor element)
	edge = mesh->peek_edge_node(p1, p2);

	// When we are on parent, we take middle vertex on the edge and add it to road_vertices. This is for consequent transformation of solution
	// on neighbor element.
	vertex = mesh->peek_vertex_node(p1, p2);
	if(vertex != NULL){
		if (n_road_vertices == 0){
			road_vertices[n_road_vertices] = vertex;
			n_road_vertices = n_road_vertices++;
		}
		else
			if(road_vertices[n_road_vertices - 1]->id != vertex->id){
				road_vertices[n_road_vertices] = vertex;
				n_road_vertices = n_road_vertices++;
			}
	}
	
	if ((edge == NULL) || (central_el->en[edge_num]->id == edge->id)){
		finding_act_elem(elem->parent, edge_num, orig_vertex_id, road_vertices, n_road_vertices);
	}
	else
		for (int i = 0; i < 2; i++)
		{
			// this condition test if on one of sides is some element and if the element is active, because it may happen that
			// something is found even thought it's not an active element
			if ((edge->elem[i] != NULL) && (edge->elem[i]->active == 1)){

				//getting to correct edge
				debug_log("way up neighbor: %d \n", edge->elem[i]->id);
				neighb_el = edge->elem[i];
				neighbor_edge = -1;
				for(int j = 0; j < neighb_el->nvert; j++)
					if(neighb_el->en[j] == edge)
						neighbor_edge = j;
				if(neighbor_edge == -1) error("edge wasn't found");

				Node* n = NULL;

				//not sure about it !!!!!!!!!!!!!!!
				n_trans[n_neighbors] = n_road_vertices;
				for(int k = 0 ; k < n_road_vertices; k++)
					debug_log("vertices on the way: %d ", road_vertices[k]->id);
				debug_log("\n");

				// go threw between elements and set correct transformation
				for(int j = n_road_vertices; j > 0; j-- ){
					if(road_vertices[j] == NULL){
//						if(j > 0)
							continue;
					}
					else{
						n = mesh->peek_vertex_node(road_vertices[j]->id, p1);
						if(n == NULL){
							n = mesh->peek_vertex_node(road_vertices[j]->id, p2);
							transformations[n_neighbors][n_road_vertices - j - 1] = neighbor_edge;
//							sol->push_transform(neighbor_edge);
							p1 = road_vertices[j]->id;
						}
						else{
								if(n->id == road_vertices[j-1]->id){
									transformations[n_neighbors][n_road_vertices - j - 1] = (neighbor_edge + 1) % neighb_el->nvert;
//									sol->push_transform((neighbor_edge + 1) % neighb_el->nvert);
									p2 = road_vertices[j]->id;
								}
								else{
									n = mesh->peek_vertex_node(road_vertices[j]->id, p2);
									transformations[n_neighbors][n_road_vertices - j - 1] = neighbor_edge;
//									sol->push_transform(neighbor_edge);
									p1 = road_vertices[j]->id;
								}
						}
					}
				}

				// final transformation on active element
				int test = 0;
				if (orig_vertex_id[0] == road_vertices[0]->id)
					test = 1;

				if(test == 1){
					transformations[n_neighbors][n_road_vertices - 1] = neighbor_edge;
//					sln->push_transform(neighbor_edge);
				}
				else{
					transformations[n_neighbors][n_road_vertices - 1] = (neighbor_edge + 1) % neighb_el->nvert;
//					sln->push_transform((neighbor_edge + 1) % neighb_el->nvert);
				}

				if(solution_flag == 1){
					// fill function values of central a neighbor element
					set_fn_values(H2D_WAY_UP);

					// set same direction as central element
					set_correct_direction(id_of_par_orient_1, id_of_par_orient_2, i);
				}
				// raise the number of neighbors
				n_neighbors = n_neighbors++;

				// add neighbor id to neighbors_id
				neighbors_id.push_back(neighb_el->id);
			}
		}
};

/*! \brief On active edge we have more neighbors. Gives us information from all neighbors.
 *
 *	We use recurrence in this way. In every step we take middle vertex of the edge (starting with active edge). This vertex split the edge
 *	on two parts. On every part (an edge) we test if the new edge is active. If not, the middle vertex is found and the method is called
 *	again with this new vertex on this part.
 */


//way down
void Neighbor::finding_act_elem( Node* vertex, int* par_vertex_id, int* road, int n_road, int use_edge, int n_vert)
{
	int son;
	int parents[2];

	Node* edge = NULL;
	Node* n = NULL;
	neighb_el = NULL;
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

				int n_road_next = n_road + 1;
				finding_act_elem( n, par_vertex_id, road, n_road_next, use_edge, n_vert);
			}
		} else
			//test if on one of sides is active element
			for (int j = 0; j < 2; j++)
			{
				if (edge->elem[j] != NULL)
					if (edge->elem[j]->active == 1){

						  debug_log("way down neighbor: %d \n", edge->elem[j]->id);
							neighb_el = mesh->get_element(edge->elem[j]->id);

							// getting to correct edge
							neighbor_edge = -1;
							for(int k = 0; k < neighb_el->nvert; k++)
								if(neighb_el->en[k] == edge)
									neighbor_edge = k;
							if(neighbor_edge == -1) error("edge wasn't found");

							// filling transformation
							for(int k = 0; k <= n_road; k++) transformations[n_neighbors][k] = road[k];

							// + 1 is because how to n_road is computed it's one less then number of transformations
							n_trans[n_neighbors] = n_road + 1;

							if(solution_flag == 1){
								// fill function values of central and neighbors elements
								set_fn_values(H2D_WAY_DOWN);

								// set same direction as central element
								set_correct_direction(parents[0], parents[1], i);
							}
							// raise number of neighbors
							n_neighbors = n_neighbors++;

							// add neighbor id to neighbors_id
							neighbors_id.push_back(neighb_el->id);
							}
			}
	}
};

/*! \brief Fill function values of central a neighbors elements
 *
 *	The flag distinguish ways and according the way it is chosen on what element are applied transformations
 */

void Neighbor::set_fn_values(Trans_flag flag){

	int number_integ_points = 0;

	switch(flag){
		case H2D_NO_TRANSF:
			{
				sol->set_active_element(neighb_el);
				neighbor_order = sol->get_fn_order();
				int max_order = std::max(central_order, neighbor_order);

				int eo = quad->get_edge_points(neighbor_edge, max_order);
				number_integ_points = quad->get_num_points(eo);

				scalar* local_fn_values_c = new scalar[number_integ_points];
				scalar* local_fn_values_n = new scalar[number_integ_points];

				// fill function values of neighbor
				sol->set_quad_order(eo);
				for(int i = 0; i < number_integ_points ; i++) local_fn_values_n[i] = sol->get_fn_values()[i];
				fn_values_neighbor[n_neighbors] = local_fn_values_n;


				//fill the central
				sol->set_active_element(central_el);
				eo = quad->get_edge_points(active_edge, max_order);
				sol->set_quad_order(eo);
				for(int i = 0; i < number_integ_points ; i++) local_fn_values_c[i] = sol->get_fn_values()[i];
				fn_values[n_neighbors] = local_fn_values_c;

				break;
			}
		case H2D_WAY_DOWN:
			{
			sol->set_active_element(neighb_el);
			neighbor_order = sol->get_fn_order();

			int max_order = std::max(central_order, neighbor_order);

			int eo = quad->get_edge_points(neighbor_edge, max_order);
			number_integ_points = quad->get_num_points(eo);

			scalar* local_fn_values_c = new scalar[number_integ_points];
			scalar* local_fn_values_n = new scalar[number_integ_points];

			// fill function values of neighbor
			sol->set_quad_order(eo);
			for(int i = 0; i < number_integ_points; i++) local_fn_values_n[i] = sol->get_fn_values()[i];

			fn_values_neighbor[n_neighbors] = local_fn_values_n;

			sol->set_active_element(central_el);

			//transform central on appropriate part
			for(int i = 0; i < n_trans[n_neighbors]; i++){
				sol->push_transform(transformations[n_neighbors][i]);
			}
			//fill the central
			eo = quad->get_edge_points(active_edge, max_order);
			sol->set_quad_order(eo);

			for(int i = 0; i < number_integ_points; i++) local_fn_values_c[i] = sol->get_fn_values()[i];

			fn_values[n_neighbors] = local_fn_values_c;

			break;
			}
		case H2D_WAY_UP:
			{
			sol->set_active_element(neighb_el);

			//transform neighbor on appropriate part
			for(int i = 0; i < n_trans[n_neighbors]; i++){
				sol->push_transform(transformations[n_neighbors][i]);
			}

			neighbor_order = sol->get_fn_order();

			int max_order = std::max(central_order, neighbor_order);

			int eo = quad->get_edge_points(neighbor_edge, max_order);
			number_integ_points = quad->get_num_points(eo);

			scalar* local_fn_values_c = new scalar[number_integ_points];
			scalar* local_fn_values_n = new scalar[number_integ_points];

			// fill function values of neighbor
			sol->set_quad_order(eo);

			for(int i = 0; i < number_integ_points; i++) local_fn_values_n[i] = sol->get_fn_values()[i];
			fn_values_neighbor[n_neighbors] = local_fn_values_n;

	
			//fill the central
			sol->set_active_element(central_el);
			eo = quad->get_edge_points(active_edge, max_order);
			sol->set_quad_order(eo);

			for(int i = 0; i < number_integ_points; i++) local_fn_values_c[i] = sol->get_fn_values()[i];

			fn_values[n_neighbors] = local_fn_values_c;

			break;
			}
		default:
			error("wasn't find scheme for getting correctly transformed function values");
	}



	// test if number of function values was assigned
	if(number_integ_points == 0)
		error("number of integration points is 0");

	np[n_neighbors] = number_integ_points;

	//reset transformations
	sol->reset_transform();

};


// correct direction, means the orientation of function values on neighbor has to be same
// as on central element

void Neighbor::set_correct_direction(int parent1, int parent2, int part_of_edge)
{

	int test = 0;
	int neighb_first_vertex = neighb_el->vn[neighbor_edge]->id;
	if (part_of_edge == 0)
	{
		if (neighb_first_vertex != parent1)
			test = 1; // means the orientation is conversely
	} else
	{
		if (neighb_first_vertex == parent2)
			test = 1; // means the orientation is conversely
	}

	if (test == 1)
	{

		scalar local_fn_value = 0;

		for (int i = 0; i < np[n_neighbors] / 2; i++){
			local_fn_value = fn_values_neighbor[n_neighbors][i];
			fn_values_neighbor[n_neighbors][i] = fn_values_neighbor[n_neighbors][np[n_neighbors] - i - 1];
			fn_values_neighbor[n_neighbors][np[n_neighbors] - i - 1] = local_fn_value;
		}
	}
};

void Neighbor::clean_all()
{
	active_edge = -1;

	n_neighbors = 0;
	neighb_el = NULL;
	neighbor_edge =-1;
	neighbor_order = -1;

	for(int i = 0; i < max_n_trans; i++)
	{
		n_trans[i] = 0;
		np[i] = 0;

		if(fn_values[i] != NULL)
			delete[] fn_values[i];
		if(fn_values_neighbor[i] != NULL)
			delete[] fn_values_neighbor[i];

		fn_values[i] = NULL;
		fn_values_neighbor[i]	= NULL;

		for(int j = 0; j < max_n_trans; j++)
			transformations[i][j] = -1;
	}
};


int Neighbor::number_of_neighbs()
{
	if(n_neighbors == 0) error("called before setting common edge");
	else return n_neighbors;
};

int* Neighbor::get_transformations(int part_edge)
{
	return transformations[part_edge];
};

scalar* Neighbor::get_fn_values_central(int part_edge)
{
	return fn_values[part_edge];
};
scalar* Neighbor::get_fn_values_neighbor(int part_edge)
{
	return fn_values_neighbor[part_edge];
};


int Neighbor::get_n_integ_points(int part_edge)
{
	return np[part_edge];
};


std::vector<int>* Neighbor::get_neighbors()
{
	return &neighbors_id;
};



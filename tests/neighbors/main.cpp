#include "hermes2d.h"
#include <iostream>

// This test tests if class neighbor finds for every active element all its neighbors


#define ERROR_SUCCESS   0
#define ERROR_FAILURE  -1

int main(int argc, char* argv[])
{
	// load mesh file
	Mesh mesh;
	H2DReader mloader;
	mloader.load("domain.mesh", &mesh);

	// perform refinement to get complex mesh
	mesh.refine_element(0, 1);
	mesh.refine_element(2);
	mesh.refine_element(3, 1);
	mesh.refine_element(7, 2);
	mesh.refine_element(13, 2);
	mesh.refine_element(8);
	mesh.refine_element(9);
	mesh.refine_element(15);
	mesh.refine_element(17);
	mesh.refine_element(24);
	mesh.refine_element(28);

  // display the mesh
//	 MeshView mview("neighbors_test", 100, 100, 500, 500);
//	 mview.show(&mesh);


	Element* e = NULL;
	Neighbor* neighb = NULL;
	std::vector<int>* neighbors_id;
	std::map<int, std::vector<int> > all_neighbors;

	// for every element we save all its neighbors into "all_neighhbors".
	// For every neighbor, his id is compared with already inserted elements into all_neighbors. If is found then
	// in his own vector of his neighbors id the id of active element is searched. Failure of the test is if
	// is no found.

	e = mesh.get_element(35);
	neighb = new Neighbor(e, &mesh);
	std::cout << e->vn[0]->id<<" "<< e->vn[1]->id<<" "<< e->vn[2]->id<<"\n";
	std::cout << e->en[1]->p1 <<" "<< e->en[1]->p2 << "\n";

	neighb->set_active_edge(1);
getchar();
	e = NULL;
	for_all_active_elements(e, &mesh)
	{
		neighb = new Neighbor(e, &mesh);
		for(int i = 0; i < e->nvert; i++){
			if(e->en[i]->bnd == 0)
				neighb->set_active_edge(i);
		}
		neighbors_id = neighb->get_neighbors();
		all_neighbors[e->id] = *neighbors_id;
		delete neighb;
	}

	int size = all_neighbors.size();
	printf("%d \n", size);

  // wait for keyboard or mouse input
  View::wait("Waiting for keyboard or mouse input.");




  // if you return this, the test will succeed:
  return ERROR_SUCCESS;
  // if you return this, the test will fail:
  //return ERROR_FAILURE;
}

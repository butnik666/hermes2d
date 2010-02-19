#include "hermes2d.h"
#include <iostream>
using namespace std;
// This example shows how to load a mesh, perform various types
// of "manual"  element refinements, and use keyboard and mouse
// controls.

void finding_act_elem(Mesh* mesh,Element* e, int edge)
{
	int son1, son2;
	int sons[2];
	int count;
	if(e->active){
		count = 0;
		printf("act neighb element %d \n", e->id);
	}
	else
	{
		count = 2;
		mesh->get_edge_sons(e, edge, son1, son2);
		sons[0] = son1;
		sons[1] = son2;
	}
	int i = 0;
	while(i < count)
	{
		finding_act_elem(mesh, e->sons[i], edge);
		i++;
	}
}



int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // perform some sample initial refinements
  mesh.refine_all_elements(2);          // Refines all elements.
//  mesh.refine_towards_vertex(2, 1);    // Refines mesh towards vertex #3 (4x).
//  mesh.refine_towards_boundary(2, 1);  // Refines all elements along boundary 2 (4x).
//  mesh.refine_all_elements();          // Refines element #86 isotropically.
//  mesh.refine_element(112, 0);         // Refines element #112 isotropically.
    mesh.refine_element(2, 1);
    mesh.refine_element(4, 0);
    // Refines element #84 anisotropically.
//	  mesh.refine_element(114, 1);         // Refines element #114 anisotropically.

  int edge;
  Element* e = NULL;
  Element* neighb;
 // for_all_active_elements(e, &mesh)
 // {
//	  for(int i = 0; i < 4; i++){
//		  if(e->en[i]->bnd == 0){
//			  printf("element id:%d , edge: %d, elem:%d \n", e->id, e->en[i]->id, e->en[i]->elem[0]->id);
//			  std::cout << e->en[i]->elem[0]<<" " << e->en[i]->elem[1]<<"\n";
/*			  neighb == NULL;
			  neighb = e->get_neighbor(i);
			  if(neighb == NULL) error("neighbor wasn't found");
			  for(int j = 0; j < 4; j++){
				  if (e->en[i] == neighb->en[j])
					  edge = j;
			  }
			  finding_act_elem(&mesh, neighb, edge);
		  }
	  }*/
//	  printf("elem: %d \n", e->id);
  //}

  e = mesh.get_element(3);
/*  Node* n;
  n = mesh.peek_vertex_node(e->en[1]->p1, e->en[1]->p2);
  cout << e->vn[0]->p1 << " "<<e->vn[0]->p2<<"\n";
  cout << mesh.peek_edge_node(e->vn[0]->p1, e->vn[0]->p2)<<"\n";
  */
//  cout << n->id << "\n";
  // display the mesh
  MeshView mview("Hello world!", 100, 100, 500, 500);  // (100, 100) is the upper left corner position
  mview.show(&mesh);                                   // 500 x 500 is the window size
  // wait for keyboard or mouse input
  View::wait("Waiting for keyboard or mouse input.");
  return 0;
}




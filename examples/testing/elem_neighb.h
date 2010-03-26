#ifndef ELEM_NEIGHB_H_
#define ELEM_NEIGHB_H_

struct Neighbor
{
	int transformations[20]; // transformations of solution on active (way down) or neighbor(way up) element
													 // the size (20) is according to max allowed transformations on one element
	int max_order; 					 // maximum of orders of solutions from both elements
	double fn_values[20];		 // function values of solution of neighbor in same orientation as active element
	int n_fn_values; 				 // number of fn_values;
};





#endif

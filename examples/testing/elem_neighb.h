struct Neighbor{
	int transformations[20];   // road for active of neighbor element to get solution on correct part
	int max_order;						 // maximum of orders from both elements 
	double fn_values[20];			 // function values from neighbor(way down) or active(way up) element
	int n_fn_values;
	};

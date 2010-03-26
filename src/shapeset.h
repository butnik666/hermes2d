// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __HERMES2D_SHAPESET_H
#define __HERMES2D_SHAPESET_H

#include "common.h"


#define check_mode      assert(mode == MODE_TRIANGLE || mode == MODE_QUAD)
#define check_vertex    assert(vertex >= 0 && vertex < nvert)
#define check_edge      assert(edge >= 0 && edge < nvert)
#define check_order(o)  assert((o) >= 0 && (o) <= max_order)
#define check_part      assert(part >= 0)
#define check_index     assert(index >= 0 && index <= max_index[mode])
#define check_component assert(component >= 0 && component < num_components)


/// \brief Defines a set of shape functions.
///
/// This class stores mainly the definitions of the polynomials for all shape functions,
/// but also their polynomial degrees, their types (vertex, edge, bubble), and contains
/// mechanisms for the calculation and storage of constrained shape functions.
///
/// The class returns shape function values for both triangles and quads, depending on
/// what mode it is in. Use the function set_mode() to switch between MODE_TRIANGLE and
/// MODE_QUAD.
///
/// Each shape function is assigned a unique number - 'index'. For standard shape functions,
/// index is positive and not greater than the value returned by get_max_index(). Negative
/// index values are reserved for constrained shape functions. These are special edge
/// functions designed to fit standard edge functions along a portion of an edge, used in
/// the construction of FE spaces on meshes with hanging nodes.
///
/// The usage of Shapeset is simple: you first obtain the index of the shape function you
/// are interested in, be it a vertex function for a given vertex, or an edge function of
/// a given order, etc. Then you call one of the functions get_fn_value(), get_dx_value(), etc.
/// All shape functions are defined on the reference domain. For triangles, this is the
/// standard triangle (-1,-1), (1,-1), (-1,1), and for quads this is the square (-1,1)^2.
///
/// The polynomial degree (or 'order') is an integer typically in the range [1-10] for H1
/// shapesets and [0-10] for H(curl) shapesets. Quadrilaterals are allowed to have different
/// orders in the x and y directions (of the reference domain). The 'order' for quads thus
/// has to be formed with the macro make_quad_order(), see common.h.
///
/// Vertex shape functions in H1 shapesets are also regarded as edge functions of orders 0
/// and 1. This simplifies constraint calculations and BC projections.
///
/// Shape functions are always real-valued.
///
class HERMES2D_API Shapeset
{
public:

  ~Shapeset() { free_constrained_edge_combinations(); }

  /// Selects MODE_TRIANGLE or MODE_QUAD.
  void set_mode(int mode)
  {
    check_mode;
    this->mode = mode;
    nvert = (mode == MODE_TRIANGLE) ? 3 : 4;
  }

  /// Returns the current mode.
  int get_mode() const { return mode; }

  /// Returns the maximum poly degree for all shape functions.
  int get_max_order() const { return max_order; }

  /// Returns the highest shape function index.
  int get_max_index() const { return max_index[mode]; }

  /// Returns 2 if this is a vector shapeset, 1 otherwise.
  int get_num_components() const { return num_components; }


  /// Returns the index of a vertex shape function associated with the specified vertex.
  int get_vertex_index(int vertex) const
  {
    check_vertex;
    return vertex_indices[mode][vertex];
  }

  /// Returns the index of an edge function associated with the specified edge and of the
  /// requested order. 'ori' can be 0 or 1 and determines edge orientation (this is for
  /// shapesets with non-symmetric edge functions).
  int get_edge_index(int edge, int ori, int order) const
  {
    check_edge; check_order(order); assert(ori == 0 || ori == 1);
    return edge_indices[mode][edge][2*order + ori];
  }

  /// Returns a complete set of indices of bubble functions for an element of the given order.
  int* get_bubble_indices(int order) const
  {
    check_order(get_h_order(order));
    check_order(get_v_order(order));
    return bubble_indices[mode][order];
  }

  /// Returns the number of bubble functions for an element of the given order.
  int get_num_bubbles(int order) const
  {
    check_order(get_h_order(order));
    check_order(get_v_order(order));
    return bubble_count[mode][order];
  }

  /// Returns the index of a constrained edge function. 'part' is 0 or 1 for edge
  /// halves, 2, 3, 4, 5 for edge quarters, etc. See shapeset.cpp.
  int get_constrained_edge_index(int edge, int order, int ori, int part) const
  {
    check_edge; check_order(order); check_part;
    assert(order <= order_mask);
    return -1 - ((part << 7) + (order << 3) + (edge << 1) + ori);
  }

  /// Returns the polynomial degree of the specified shape function.
  /// If on quads, it returns encoded orders. The orders has to be decoded through macros
  /// get_h_order and get_v_order.
  int get_order(int index) const
  {
    if (index >= 0) {
      check_index;
      return index_to_order[mode][index];
    }
    else return ((-1 - index) >> 3) & 15;
  }


  /// Obtains the value of the given shape function. (x,y) is a coordinate in the reference
  /// domain, component is 0 for scalar shapesets and 0 or 1 for vector shapesets.
  inline double get_value(int n, int index, double x, double y, int component)
  {
    if (index >= 0)
    {
      check_index; check_component;
      return shape_table[n][mode][component][index](x, y);
    }
    else
      return get_constrained_value(n, index, x, y, component);
  }

  inline double get_fn_value (int index, double x, double y, int component) { return get_value(0, index, x, y, component); }
  inline double get_dx_value (int index, double x, double y, int component) { return get_value(1, index, x, y, component); }
  inline double get_dy_value (int index, double x, double y, int component) { return get_value(2, index, x, y, component); }
  inline double get_dxx_value(int index, double x, double y, int component) { return get_value(3, index, x, y, component); }
  inline double get_dyy_value(int index, double x, double y, int component) { return get_value(4, index, x, y, component); }
  inline double get_dxy_value(int index, double x, double y, int component) { return get_value(5, index, x, y, component); }


  /// Returns the coordinates of the reference domain vertices.
  double2* get_ref_vertex(int vertex)
  {
    return &ref_vert[mode][vertex];
  }

  /// Shape-function function type. Internal.
  typedef double (*shape_fn_t)(double, double);

  /// Returns shapeset identifier. Internal.
  virtual int get_id() const = 0;


protected:

  int mode;
  int nvert;

  shape_fn_t*** shape_table[6];

  int**  vertex_indices;
  int*** edge_indices;
  int*** bubble_indices;
  int**  bubble_count;
  int**  index_to_order;

  double2 ref_vert[2][4];
  int max_order;
  int max_index[2];
  int num_components;

  int ebias; ///< 2 for H1 shapesets, 0 for H(curl) shapesets. It is the order of the
             ///< first edge function.

  double** comb_table;
  int table_size;

  double* calculate_constrained_edge_combination(int order, int part, int ori);
  double* get_constrained_edge_combination(int order, int part, int ori, int& nitems);

  void    free_constrained_edge_combinations();

  double get_constrained_value(int n, int index, double x, double y, int component);

};


/// Creates a new shapeset. Internal. ("Shapeset factory")
//Shapeset* hermes2d_create_shapeset(int id);


// TODO : promyslet moznost ulozeni shapesetu jako tabulky monomialnich koeficientu
// - mozna efektivnejsi nez stavajici kod
// - monolitictejsi, elegantnejsi
// - pri ulozeni koeficientu jako long double mozna i presnejsi
// - moznost ulozeni jen zakladnich polynomu, derivace lze dopocitat automaticky


#undef check_mode
#undef check_vertex
#undef check_edge
#undef check_order
#undef check_part
#undef check_index
#undef check_component

#endif

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

#include "common.h"
#include "linear.h"
#include "refmap.h"


#include "linear_data.cpp"

// vertices
static int*      ord_np[2]          = { num_vert_tri, num_vert_quad };
static double3*  ord_tables_tri[2]  = { vert_tri0, vert_tri1 };
static double3*  ord_tables_quad[2] = { vert_quad0, vert_quad1 };
static double3** ord_tables[2]      = { ord_tables_tri, ord_tables_quad };

// triangles
static int*      num_elem[2]        = { num_elem_tri, num_elem_quad};
static int3*     ord_elem_tri[2]    = { elem_tri0, elem_tri1 };
static int3*     ord_elem_quad[2]   = { elem_quad0,  elem_quad1 };
static int3**    ord_elem[2]        = { ord_elem_tri, ord_elem_quad };

// edges
static int*      num_edge[2]        = { num_edge_tri, num_edge_quad};
static int3*     ord_edge_tri[2]    = { edge_tri0, edge_tri1 };
static int3*     ord_edge_quad[2]   = { edge_quad0,  edge_quad1 };
static int3**    ord_edge[2]        = { ord_edge_tri, ord_edge_quad };


static class Quad2DOrd : public Quad2D
{
public:

  Quad2DOrd()
  {
    mode = MODE_TRIANGLE;
    max_order[0] = max_order[1] = 1;
    num_tables[0] = num_tables[1] = 2;
    tables = ord_tables;
    np = ord_np;
  };

  virtual void dummy_fn() {}

} quad_ord;


//// Orderizer /////////////////////////////////////////////////////////////////////////////////////

Orderizer::Orderizer()
         : Linearizer()
{
  ltext = NULL;
  lvert = NULL;
  lbox = NULL;

  nl = cl1 = cl2 = cl3 = 0;

  for (int i = 0, p = 0; i <= 10; i++)
  {
    for (int j = 0; j <= 10; j++)
    {
      assert(p < sizeof(buffer)-5);
      sprintf(buffer+p, (i == j) ? "%d" : "%d|%d", i, j);
      labels[i][j] = buffer+p;
      p += strlen(buffer+p) + 1;
    }
  }
}


void Orderizer::process_solution(Space* space)
{
  if (!space->is_up_to_date())
    error("The space is not up to date.");

  int type = 1;

  nv = nt = ne = nl = 0;
  del_slot = -1;

  // estimate the required number of vertices and triangles
  Mesh* mesh = space->get_mesh();
  int nn = mesh->get_num_active_elements();
  int ev = 77 * nn, et = 64 * nn, ee = 16 * nn, el = nn + 10;

  // reuse or allocate vertex, triangle and edge arrays
  lin_init_array(verts, double3, cv, ev);
  lin_init_array(tris, int3, ct, et);
  lin_init_array(edges, int3, ce, ee);
  lin_init_array(lvert, int, cl1, el);
  lin_init_array(ltext, char*, cl2, el);
  lin_init_array(lbox, double2, cl3, el);
  info = NULL;

  int oo, o[6];

  RefMap refmap;
  refmap.set_quad_2d(&quad_ord);

  // make a mesh illustrating the distribution of polynomial orders over the space
  Element* e;
  for_all_active_elements(e, mesh)
  {
    oo = o[4] = o[5] = space->get_element_order(e->id);
    for (unsigned int k = 0; k < e->nvert; k++)
      o[k] = space->get_edge_order(e, k);

    refmap.set_active_element(e);
    double* x = refmap.get_phys_x(type);
    double* y = refmap.get_phys_y(type);

    double3* pt = quad_ord.get_points(type);
    int np = quad_ord.get_num_points(type);
    int id[80];
    assert(np <= 80);

    #define make_vert(index, x, y, val) \
      { (index) = add_vertex(); \
      verts[index][0] = (x); \
      verts[index][1] = (y); \
      verts[index][2] = (val); }

    int mode = e->get_mode();
    if (e->is_quad())
    {
      o[4] = get_h_order(oo);
      o[5] = get_v_order(oo);
    }
    make_vert(lvert[nl], x[0], y[0], o[4]);

    for (int i = 1; i < np; i++)
      make_vert(id[i-1], x[i], y[i], o[(int) pt[i][2]]);

    for (int i = 0; i < num_elem[mode][type]; i++)
      add_triangle(id[ord_elem[mode][type][i][0]], id[ord_elem[mode][type][i][1]], id[ord_elem[mode][type][i][2]]);

    for (int i = 0; i < num_edge[mode][type]; i++)
    {
      if (e->en[ord_edge[mode][type][i][2]]->bnd || (y[ord_edge[mode][type][i][0] + 1] < y[ord_edge[mode][type][i][1] + 1]) ||
          ((y[ord_edge[mode][type][i][0] + 1] == y[ord_edge[mode][type][i][1] + 1]) &&
           (x[ord_edge[mode][type][i][0] + 1] <  x[ord_edge[mode][type][i][1] + 1])))
      {
        add_edge(id[ord_edge[mode][type][i][0]], id[ord_edge[mode][type][i][1]], 0);
      }
    }

    double xmin = 1e100, ymin = 1e100, xmax = -1e100, ymax = -1e100;
    for (unsigned int k = 0; k < e->nvert; k++)
    {
      if (e->vn[k]->x < xmin) xmin = e->vn[k]->x;
      if (e->vn[k]->x > xmax) xmax = e->vn[k]->x;
      if (e->vn[k]->y < ymin) ymin = e->vn[k]->y;
      if (e->vn[k]->y > ymax) ymax = e->vn[k]->y;
    }
    lbox[nl][0] = xmax - xmin;
    lbox[nl][1] = ymax - ymin;
    ltext[nl++] = labels[o[4]][o[5]];
  }

  refmap.set_quad_2d(&g_quad_2d_std);
}


Orderizer::~Orderizer()
{
  lin_free_array(lvert, nl, cl1);
  lin_free_array(ltext, nl, cl2);
  lin_free_array(lbox, nl, cl3);
}


//// save & load ///////////////////////////////////////////////////////////////////////////////////

void Orderizer::save_data(const char* filename)
{
  FILE* f = fopen(filename, "wb");
  if (f == NULL) error("Could not open %s for writing.", filename);
  lock_data();

  AUTOLA_OR(int, orders, nl);
  int vo, ho;
  for (int i = 0; i < nl; i++)
  {
    if (strchr(ltext[i], '|'))
      sscanf(ltext[i], "%d|%d", &ho, &vo);
    else
      { sscanf(ltext[i], "%d", &ho); vo = ho; }
    orders[i] = make_quad_order(ho, vo);
  }

  if (fwrite("H2DO\001\000\000\000", 1, 8, f) != 8 ||
      fwrite(&nv, sizeof(int), 1, f) != 1 ||
      fwrite(verts, sizeof(double3), nv, f) != nv ||
      fwrite(&nt, sizeof(int), 1, f) != 1 ||
      fwrite(tris, sizeof(int3), nt, f) != nt ||
      fwrite(&ne, sizeof(int), 1, f) != 1 ||
      fwrite(edges, sizeof(int3), ne, f) != ne ||
      fwrite(&nl, sizeof(int), 1, f) != 1 ||
      fwrite(lvert, sizeof(int), nl, f) != nl ||
      fwrite(lbox, sizeof(double2), nl, f) != nl ||
      fwrite(orders, sizeof(int), nl, f) != nl)
  {
    error("Error writing data to %s", filename);
  }

  unlock_data();
  fclose(f);
}


void Orderizer::load_data(const char* filename)
{
  FILE* f = fopen(filename, "rb");
  if (f == NULL) error("Could not open %s for reading.", filename);
  lock_data();

  struct { char magic[4]; int ver; } hdr;
  if (fread(&hdr, sizeof(hdr), 1, f) != 1)
    error("Error reading %s", filename);

  if (hdr.magic[0] != 'H' || hdr.magic[1] != '2' || hdr.magic[2] != 'D' || hdr.magic[3] != 'O')
    error("%s is not a Hermes2D Orderizer file.", filename);
  if (hdr.ver > 1)
    error("%s -- unsupported file version.", filename);

  #define read_array(array, type, n, c, what) \
    if (fread(&n, sizeof(int), 1, f) != 1) \
      error("Error reading the number of " what " from %s", filename); \
    lin_init_array(array, type, c, n); \
    if (fread(array, sizeof(type), n, f) != n) \
      error("Error reading " what " from %s", filename);

  read_array(verts, double3, nv, cv,  "vertices");
  read_array(tris,  int3,    nt, ct,  "triangles");
  read_array(edges, int3,    ne, ce,  "edges");
  read_array(lvert, int,     nl, cl1, "label vertices");

  lin_init_array(lbox, double2, cl3, nl);
  if (fread(lbox, sizeof(double2), nl, f) != nl)
    error("Error reading label bounding boxes from %s", filename);

  AUTOLA_OR(int, orders, nl);
  if (fread(orders, sizeof(int), nl, f) != nl)
    error("Error reading element orders from %s", filename);

  lin_init_array(ltext, char*, cl2, nl);
  for (int i = 0; i < nl; i++)
    ltext[i] = labels[get_h_order(orders[i])][get_v_order(orders[i])];

  find_min_max();
  unlock_data();
  fclose(f);
}

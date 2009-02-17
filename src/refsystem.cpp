// This file is part of Hermes2D.
//
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@utep.edu>
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
#include "space.h"
#include "weakform.h"
#include "refsystem.h"


  
RefSystem::RefSystem(LinSystem* base, int order_increase, int refinement)
         : LinSystem(base->wf, base->solver)
{
  this->base = base;
  order_inc = order_increase;
  this->refinement = refinement;
  ref_meshes = NULL;
  ref_spaces = NULL;
}

RefSystem::~RefSystem()
{
  free_ref_data();
}


void RefSystem::set_spaces(int n, ...)
{
  error("set_spaces must not be called in RefSystem.");
}

void RefSystem::set_pss(int n, ...)
{
  error("set_pss must not be called in RefSystem.");
}


void RefSystem::assemble(bool rhsonly)
{
  int i, j;
  
  // get rid of any previous data
  free_ref_data();
  
  ref_meshes = new Mesh*[wf->neq];
  ref_spaces = new Space*[wf->neq];
  
  // copy meshes from the coarse problem, refine them
  for (i = 0; i < wf->neq; i++)
  {
    Mesh* mesh = base->spaces[i]->get_mesh();
    
    // check if we already have the same mesh
    for (j = 0; j < i; j++)
      if (mesh->get_seq() == base->spaces[j]->get_mesh()->get_seq())
        break;
      
    if (j < i) // yes
    {
      ref_meshes[i] = ref_meshes[j];
    }
    else // no, copy and refine the coarse one
    {
      Mesh* rmesh = new Mesh;
      rmesh->copy(mesh);
      if (refinement ==  1) rmesh->refine_all_elements();
      if (refinement == -1) rmesh->unrefine_all_elements();
      ref_meshes[i] = rmesh;
    }
  }
  
  // duplicate spaces from the coarse problem, assign reference orders and dofs
  int dofs = 0;
  for (i = 0; i < wf->neq; i++)
  {
    ref_spaces[i] = base->spaces[i]->dup(ref_meshes[i]);
    if (refinement == -1)
    {
      Element* re;
      for_all_active_elements(re, ref_meshes[i])
      {
        Mesh* mesh = base->spaces[i]->get_mesh();
        Element* e = mesh->get_element(re->id);
        int max_o = 0;
        if (e->active) 
          max_o = get_h_order(base->spaces[i]->get_element_order(e->id));
        else
          for (int son = 0; son < 4; son++)
          {
            int o = get_h_order(base->spaces[i]->get_element_order(e->sons[son]->id));
            if (o > max_o) max_o = o;
          }
        ref_spaces[i]->set_element_order(re->id, max_o + order_inc);
      }

    }
    else
      ref_spaces[i]->copy_orders(base->spaces[i], order_inc);

    dofs += ref_spaces[i]->assign_dofs(dofs);
  }
  
  memcpy(spaces, ref_spaces, sizeof(Space*) * wf->neq);
  memcpy(pss, base->pss, sizeof(PrecalcShapeset*) * wf->neq);
  have_spaces = true;

  LinSystem::assemble(rhsonly);
}


void RefSystem::free_ref_data()
{
  int i, j;
  
  // free reference meshes
  if (ref_meshes != NULL)
  {
    for (i = 0; i < wf->neq; i++)
    {
      for (j = 0; j < i; j++)
        if (ref_meshes[j] == ref_meshes[i])
          break;
      
      if (i == j) delete ref_meshes[i];
    }
    
    delete [] ref_meshes;
    ref_meshes = NULL;
  }
  
  // free reference spaces
  if (ref_spaces != NULL)
  {
    for (i = 0; i < wf->neq; i++)
      delete ref_spaces[i];
    
    delete [] ref_spaces;
    ref_spaces = NULL;
  }
}

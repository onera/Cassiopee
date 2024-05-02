/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "proto.h"

#define one_sixth 0.16666666666666666666

// assumes quad
void compute_face_center(mesh *M, E_Int face)
{
  E_Float *fx = &M->fc[3*face];
  for (E_Int i = 0; i < 3; i++)
    fx[i] = 0.;

  E_Int *pn = &M->NGON[4*face];

  for (E_Int i = 0; i < 4; i++) {
    E_Float *px = &M->xyz[3*pn[i]];
    for (E_Int j = 0; j < 3; j++)
      fx[j] += px[j];
  }

  for (int i = 0; i < 3; i++)
    fx[i] *= 0.25;
}

// assumes hexahedron
void compute_cell_center(mesh *M, E_Int cell)
{
  E_Float *cx = &M->cc[3*cell];
  for (E_Int i = 0; i < 3; i++)
    cx[i] = 0;

  E_Int *pf = &M->NFACE[6*cell];

  for (E_Int i = 0; i < 6; i++) {
    E_Float *fx = &M->fc[3*pf[i]];
    for (E_Int j = 0; j < 3; j++)
      cx[j] += fx[j];
  }

  for (E_Int i = 0; i < 3; i++)
    cx[i] *= one_sixth;
}


void compute_face_centers(mesh *M)
{
  if (M->fc) return;
  M->fc = (E_Float *)XCALLOC(3*M->nfaces, sizeof(E_Float));

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Float *pf = &M->fc[3*i];

    E_Int start = M->xfaces[i];
    E_Int end = M->xfaces[i+1];
    
    for (E_Int j = start; j < end; j++) {
      E_Int pt = M->NGON[j];
      E_Float *px = &M->xyz[3*pt];
      pf[0] += px[0];
      pf[1] += px[1];
      pf[2] += px[2];
    }

    E_Float coeff = 1./(end-start);

    for (E_Int k = 0; k < 3; k++)
      pf[k] *= coeff;
  }
}

void compute_cell_centers(mesh *M)
{
  if (M->cc) return;
  M->cc = (E_Float *)XCALLOC(3*M->ncells, sizeof(E_Float));
  compute_face_centers(M);

  E_Int *po = M->owner;
  E_Int *pn = M->neigh;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Float *fx = &M->fc[3*i];
    
    E_Int own = po[i];
    E_Float *co = &M->cc[3*own];

    for (E_Int j = 0; j < 3; j++)
      co[j] += fx[j];

    E_Int nei = pn[i];
    if (nei < 0) continue;
    E_Float *cn = &M->cc[3*nei];

    for (E_Int j = 0; j < 3; j++)
      cn[j] += fx[j];
  }
  
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Float coeff = 1./(M->xcells[i+1] - M->xcells[i]);
    E_Float *px = &M->cc[3*i];
    for (E_Int j = 0; j < 3; j++)
      px[j] *= coeff;
  }
}

// Zero-based face and cell
E_Int get_neighbour(E_Int cell, E_Int face, mesh *M)
{
  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;
  assert(cell == owner[face] || cell == neigh[face]);
  if (cell == owner[face]) return neigh[face];
  return owner[face];
}

static
void ppatch_free(proc_patch *pp)
{
  XFREE(pp->faces);
  XFREE(pp->gneis);
  XFREE(pp->send_buf_d);
  XFREE(pp->send_buf_i);
  XFREE(pp->recv_buf_d);
  XFREE(pp->recv_buf_i);
}

void mesh_free(mesh *M) {
  XFREE(M->xyz);
  XFREE(M->owner);
  XFREE(M->neigh);
  XFREE(M->NFACE);
  XFREE(M->xcells);
  XFREE(M->NGON);
  XFREE(M->xfaces);
  M->CT.clear();
  M->FT.clear();
  M->PT.clear();
  M->ET.clear();
  XFREE(M->cc);
  XFREE(M->fc);
  XFREE(M->lsqG);
  XFREE(M->lsqGG);
  XFREE(M->lsqH);
  XFREE(M->lsqHH);
  for (E_Int i = 0; i < M->nppatches; i++)
    ppatch_free(&M->ppatches[i]);
  delete [] M->ppatches;
  XFREE(M->req);
  if (M->pnei_coords) {
    for (E_Int i = 0; i < M->nppatches; i++)
      XFREE(M->pnei_coords[i]);
    XFREE(M->pnei_coords);
  }
  if (M->pnei_flds) {
    for (E_Int i = 0; i < M->nppatches; i++)
      XFREE(M->pnei_flds[i]);
    XFREE(M->pnei_flds);
  }
  XFREE(M->pnei_grads);
  XFREE(M->gcells);
  XFREE(M->gfaces);
  XFREE(M->gpoints);
  XFREE(M->ref_data);

  delete M;
}

/*
size_t umap_size(const std::unordered_map<E_Int, E_Int>& umap)
{
  size_t count = 0;
  for (size_t i = 0; i < umap.bucket_count(); i++) {
    size_t bucket_size = umap.bucket_size(i);
    if (bucket_size == 0)
      count++;
    else
      count += bucket_size;
  }
  return count;
}
*/

E_Float mesh_memsize(mesh *M)
{
  size_t memsize = 0;

  E_Int pc = M->predicted_ncells;
  E_Int pf = M->predicted_nfaces;
  E_Int pp = M->predicted_npoints;

  // points
  memsize += pp*3*sizeof(E_Float);

  // cells
  memsize += 6*pc*sizeof(E_Int); // NFACE
  memsize += pc+1*sizeof(E_Int); // xcells
  memsize += 3*pc*sizeof(E_Float); // cc

  // faces
  memsize += 2*pf*sizeof(E_Int); // owner + neigh
  memsize += 4*pf*sizeof(E_Int); // NGON
  memsize += (pf+1)*sizeof(E_Int); // xfaces
  memsize += 3*pf*sizeof(E_Float); // fc

  return memsize/1000000.; // in mega bytes
}

void mesh_save_memory(mesh *M)
{
  M->NGON  = (E_Int *)   XRESIZE(M->NGON,   (4*M->nfaces)  * sizeof(E_Int));
  M->NFACE  = (E_Int *)  XRESIZE(M->NFACE,  (6*M->ncells)  * sizeof(E_Int));
  M->fc     = (E_Float *)XRESIZE(M->fc,     (3*M->nfaces)  * sizeof(E_Float));
  M->cc     = (E_Float *)XRESIZE(M->cc,     (3*M->ncells)  * sizeof(E_Float));
  M->xyz    = (E_Float *)XRESIZE(M->xyz,    (3*M->npoints) * sizeof(E_Float));
  M->owner  = (E_Int *)  XRESIZE(M->owner,  M->nfaces      * sizeof(E_Int));
  M->neigh  = (E_Int *)  XRESIZE(M->neigh,  M->nfaces      * sizeof(E_Int));
  M->xfaces = (E_Int *)  XRESIZE(M->xfaces, (M->nfaces+1)  * sizeof(E_Int));
  M->xcells = (E_Int *)  XRESIZE(M->xcells, (M->ncells+1)  * sizeof(E_Int));

  M->predicted_ncells = M->ncells;
  M->predicted_nfaces = M->nfaces;
  M->predicted_npoints = M->npoints;
}

mesh_leaves mesh_get_leaves(mesh *M, tree *ct, tree *ft, E_Int conformize)
{
  mesh_leaves cM;
  auto &NFACE   = cM.NFACE;
  auto &XCELLS  = cM.XCELLS;
  auto &NGON    = cM.NGON;
  auto &XFACES  = cM.XFACES;
  auto &XYZ     = cM.XYZ;
  auto &ecells  = cM.ecells;
  auto &efaces  = cM.efaces;
  auto &epoints = cM.epoints;
  auto &necells = cM.necells;
  auto &nefaces = cM.nefaces;
  auto &nepoints = cM.nepoints;

  // extract enabled cells
  for (E_Int i = 0; i < M->ncells; i++) {
    if (ct->enabled[i])
      ecells.push_back(i);
  }
  necells = ecells.size();

  // replace faces with leaves
  XCELLS.resize(1, 0);
  XFACES.resize(1, 0);

  nefaces = 0;
  nepoints = 0;

  if (conformize) {
    for (const auto& cell : ecells) {
      E_Int stride = 0;
      for (E_Int j = M->xcells[cell]; j < M->xcells[cell+1]; j++) {
        E_Int face = M->NFACE[j];
        std::vector<E_Int> leaves;
        tree_get_face_leaves(ft, face, leaves);
        stride += leaves.size();
        for (const auto& leaf : leaves)
          NFACE.push_back(leaf);
      }
      XCELLS.push_back(stride);
    }
  } else {
    for (const auto& cell : ecells) {
      for (E_Int j = M->xcells[cell]; j < M->xcells[cell+1]; j++) {
        E_Int face = M->NFACE[j];
        NFACE.push_back(face);
      }
      XCELLS.push_back(M->xcells[cell+1] - M->xcells[cell]);
    }
  }

  for (E_Int i = 0; i < necells; i++)
        XCELLS[i+1] += XCELLS[i];

  for (E_Int i = 0; i < necells; i++) {
    for (E_Int j = XCELLS[i]; j < XCELLS[i+1]; j++) {
      E_Int face = NFACE[j];
      if (efaces.find(face) == efaces.end()) {
        efaces[face] = nefaces++;
        for (E_Int k = M->xfaces[face]; k < M->xfaces[face+1]; k++) {
          E_Int point = M->NGON[k];
          NGON.push_back(point);
        }
        XFACES.push_back(M->xfaces[face+1] - M->xfaces[face]);
      }
    }
  }

  for (E_Int i = 0; i < nefaces; i++)
    XFACES[i+1] += XFACES[i];
  
  for (E_Int i = 0; i < nefaces; i++) {
    for (E_Int j = XFACES[i]; j < XFACES[i+1]; j++) {
      E_Int point = NGON[j];
      if (epoints.find(point) == epoints.end())
        epoints[point] = nepoints++;
    }
  }
  
  XYZ.resize(3*nepoints);
  for (const auto& point : epoints) {
    E_Float *ox = &M->xyz[3*point.first];
    E_Float *nx = &XYZ[3*point.second];
    for (E_Int i = 0; i < 3; i++)
      nx[i] = ox[i];
  }

  return cM;
}

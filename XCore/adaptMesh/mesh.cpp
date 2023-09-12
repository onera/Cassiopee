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
  for (E_Int i = 0; i < M->nppatches; i++) {
    XFREE(M->pnei_coords[i]);
    XFREE(M->pnei_flds[i]);
    //XFREE(M->pnei_grads[i]);
  }
  XFREE(M->pnei_coords);
  XFREE(M->pnei_flds);
  XFREE(M->pnei_grads);
  XFREE(M->gcells);
  XFREE(M->gfaces);
  XFREE(M->gpoints);
  XFREE(M->ref_data);

  delete M;
}
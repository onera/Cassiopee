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
#include "Proto.h"
#include <cassert>

AMesh *init_mesh(K_FLD::FldArrayI &cn, E_Float *px, E_Float *py,
  E_Float *pz, E_Int npts)
{
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int ncells = cn.getNElts();
  E_Int nfaces = cn.getNFaces();

  AMesh *M = new AMesh;

  M->npoints = npts;
  M->nfaces = nfaces;
  M->ncells = ncells;

  M->x = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));
  M->y = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));
  M->z = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));

  memcpy(M->x, px, M->npoints * sizeof(E_Float));
  memcpy(M->y, py, M->npoints * sizeof(E_Float));
  memcpy(M->z, pz, M->npoints * sizeof(E_Float));

  M->indPH = (E_Int *)XMALLOC((M->ncells+1) * sizeof(E_Int));
  M->indPG = (E_Int *)XMALLOC((M->nfaces+1) * sizeof(E_Int));

  memcpy(M->indPH, indPH, (M->ncells+1) * sizeof(E_Int));
  memcpy(M->indPG, indPG, (M->nfaces+1) * sizeof(E_Int));

  M->nface = (E_Int *)XMALLOC(M->indPH[M->ncells] * sizeof(E_Int));
  M->ngon  = (E_Int *)XMALLOC(M->indPG[M->nfaces] * sizeof(E_Int));

  E_Int *ptr = M->nface;

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int *pf = cn.getElt(i, nf, nface, indPH);
    for (E_Int j = 0; j < nf; j++)
      *ptr++ = pf[j]-1;
  }

  ptr = M->ngon;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = cn.getFace(i, np, ngon, indPG);
    for (E_Int j = 0; j < np; j++)
      *ptr++ = pn[j]-1;
  }

  M->ecenter = new std::map<Edge, E_Int>;

  M->CT = new std::unordered_map<E_Int, E_Int>;
  M->FT = new std::unordered_map<E_Int, E_Int>;
  M->PT = new std::unordered_map<E_Int, E_Int>;
 
  return M;
}

E_Int master_face(E_Int face, AMesh *M)
{
  return M->faceTree->children(face) ? face : M->faceTree->parent(face);
}

E_Int master_cell(E_Int cell, AMesh *M)
{
  if (cell == -1) return -1;
  return M->cellTree->children(cell) ? cell : M->cellTree->parent(cell);
}

const char *cell_type_to_string(E_Int type)
{
  switch (type) {
    case HEXA: return "HEXA";
    case TETRA: return "TETRA";
    case PENTA: return "PENTA";
    case PYRA: return "PYRA";
    default:
      return "UNKNOWN";
  }
}

const char *face_type_to_string(E_Int type)
{
  switch (type) {
    case HEXA: return "QUAD";
    case TRI: return "TRI";
    default:
      return "UNKNOWN";
  }
}

void ngon_print(AMesh *M)
{
  puts("");
  for (E_Int i = 0; i < M->nfaces; i++) {
    for (E_Int j = M->indPG[i]; j < M->indPG[i+1]; j++)
      printf(SF_D_ " ", M->ngon[j]);
    puts("");
  }
  puts("");
}

void nface_print(AMesh *M)
{
  puts("");
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->indPH[i]; j < M->indPH[i+1]; j++)
      printf(SF_D_ " ", M->nface[j]);
    puts("");
  }
  puts("");
}

const E_Int normalIn_T[4] = {1, 0, 1, 0};
const E_Int normalIn_Pe[5] = {1, 0, 1, 0, 1};
const E_Int normalIn_Py[5] = {1, 1, 0, 1, 0};
const E_Int normalIn_H[6] = {1, 0, 1, 0, 1, 0};

Edge::Edge()
{}

Edge::Edge(E_Int p0, E_Int p1) :
  p0_(std::min(p0, p1)), p1_(std::max(p0, p1))
{}

void Edge::set(E_Int p0, E_Int p1)
{
  p0_ = std::min(p0, p1);
  p1_ = std::max(p0, p1);
}

bool Edge::operator<(const Edge &e) const
{
  return (p0_ < e.p0_) || (p0_ == e.p0_ && p1_ < e.p1_);
}

void patch_drop(Patch *P)
{
  XFREE(P->pf);
  XFREE(P->pn);
  XFREE(P->sbuf_i);
  XFREE(P->rbuf_i);
  XFREE(P->sbuf_f);
  XFREE(P->rbuf_f);
}

AMesh::AMesh() :
  ncells(-1), nfaces(-1), npoints(-1), nif(-1), nbf(-1), npf(-1),
  x(NULL), y(NULL), z(NULL),
  nface(NULL), indPH(NULL), ngon(NULL), indPG(NULL),
  owner(NULL), neigh(NULL),
  nbc(-1), ptlists(NULL), bcsizes(NULL), bcnames(NULL),
  Tr(-1.0), Tu(-1.0), eps(-1.0), hmin(-1.0), hmax(-1.0), unrefine(-1),
  mode_2D(NULL), ref_data(NULL), ecenter(NULL), ccenter(),
  cellTree(NULL), faceTree(NULL),
  prev_ncells(-1), prev_nfaces(-1), prev_npoints(-1),
  onc(-1), onf(-1), onp(-1),
  pid(-1), npc(-1), nrq(-1), req(NULL),
  gcells(NULL), gfaces(NULL), gpoints(NULL),
  npatches(-1), patches(NULL),
  PT(NULL), FT(NULL), CT(NULL),
  XNEIS(NULL), CADJ(NULL),
  closed_indPG(NULL), closed_ngon(NULL),
  closed_indPH(NULL), closed_nface(NULL)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &npc);
  nrq = 0;
  req = (MPI_Request *)XMALLOC(2*npc * sizeof(MPI_Request));
}

void mesh_drop(AMesh *M)
{
  XFREE(M->x);
  XFREE(M->y);
  XFREE(M->z);

  XFREE(M->nface);
  XFREE(M->indPH);
  XFREE(M->ngon);
  XFREE(M->indPG);

  XFREE(M->owner);
  XFREE(M->neigh);

  for (E_Int i = 0; i < M->nbc; i++) {
    delete M->ptlists[i];
    XFREE(M->bcnames[i]);
  }
  XFREE(M->ptlists);
  XFREE(M->bcnames);
  XFREE(M->bcsizes);

  XFREE(M->mode_2D);
  XFREE(M->ref_data);

  delete M->ecenter;
  
  M->cellTree->drop();
  M->faceTree->drop();
  delete M->cellTree;
  delete M->faceTree;

  XFREE(M->req);
  XFREE(M->gcells);
  XFREE(M->gfaces);
  XFREE(M->gpoints);

  for (E_Int i = 0; i < M->npatches; i++)
    patch_drop(&M->patches[i]);
  XFREE(M->patches);

  delete M->PT;
  delete M->FT;
  delete M->CT;

  XFREE(M->XNEIS);
  XFREE(M->CADJ);
  
  XFREE(M->closed_indPG);
  XFREE(M->closed_ngon);
  XFREE(M->closed_indPH);
  XFREE(M->closed_nface);

  delete M;
}

E_Int get_neighbour(E_Int cell, E_Int face, AMesh *M)
{
  assert(cell == M->owner[face] || cell == M->neigh[face]);
  if (cell == M->owner[face]) return M->neigh[face];
  return M->owner[face];
}

E_Int *get_face(E_Int i, E_Int &np, E_Int *ngon, E_Int *indPG)
{
  np = indPG[i+1] - indPG[i];
  return &ngon[indPG[i]];
}

E_Int *get_cell(E_Int i, E_Int &nf, E_Int *nface, E_Int *indPH)
{
  nf = indPH[i+1] - indPH[i];
  return &nface[indPH[i]];
}

E_Int get_stride(E_Int i, E_Int *indir)
{
  return indir[i+1] - indir[i];
}

E_Int *get_facets(E_Int i, E_Int *cn, E_Int *indir)
{
  return &cn[indir[i]];
}

E_Int get_reorient(E_Int face, E_Int cell, E_Int normalIn, AMesh *M)
{
  assert(M->owner[face] == cell || M->neigh[face] == cell);
  if (M->neigh[face] == cell && normalIn == 1) return 0;
  else if (M->owner[face] == cell && normalIn == 0) return 0;
  else return 1;
}

E_Int is_internal_face(E_Int face, AMesh *M)
{
  return M->neigh[face] != -1;
}

void get_full_cell(E_Int cell, AMesh *M, E_Int &nf, E_Int *pf)
{
  E_Int stride = -1;
  E_Int *FACES = get_cell(cell, stride, M->nface, M->indPH);
  nf = 0;

  for (E_Int i = 0; i < stride; i++) {
    E_Int face = FACES[i];

    E_Int flvl = M->faceTree->level(face);
    E_Int clvl = M->cellTree->level(cell);

    assert(flvl >= clvl);

    if (flvl == clvl) {
      pf[nf++] = face;
    } else {
      assert(master_face(face, M) == face);
      Children *children = M->faceTree->children(face);
      for (E_Int j = 0; j < children->n; j++)
        pf[nf++] = children->pc[j];
    }
  }
}

void reconstruct_parent_quad(E_Int face, AMesh *M, E_Int pn[4])
{
  Children *children = M->faceTree->children(face);
  if (children == NULL) {
    memcpy(pn, &M->ngon[M->indPG[face]], 4*sizeof(E_Int));
  } else {
    if (children->n == 2) {
      E_Int *pn0 = get_facets(children->pc[0], M->ngon, M->indPG);
      E_Int *pn1 = get_facets(children->pc[1], M->ngon, M->indPG);
      E_Int dir = pn0[1] == pn1[0] ? DIRX : DIRY;
      if (dir == DIRX) {
        pn[0] = pn0[0]; pn[1] = pn1[1]; pn[2] = pn1[2]; pn[3] = pn0[3];
      } else {
        pn[0] = pn0[0]; pn[1] = pn0[1]; pn[2] = pn1[2]; pn[3] = pn1[3];
      }
    } else {
      for (E_Int i = 0; i < children->n; i++)
        pn[i] = get_facets(children->pc[i], M->ngon, M->indPG)[i];
    }
  }
}

// Keep refined faces and their children, untouched faces and unrefined
// faces

void update_boundary_faces(AMesh *M)
{
  M->nbf = 0;
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];
    
    // How many faces on current boundary have been refined/unrefined?
    E_Int new_bcsize = 0;
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      E_Int state = M->faceTree->state(face);

      if (state == UNTOUCHED)
        new_bcsize += 1;
      else if (state == REFINED)
        new_bcsize += M->faceTree->children(face)->n;
      else if (state == UNREFINED)
        new_bcsize += 1;
    }

    M->nbf += new_bcsize;

    E_Int *new_ptlist = (E_Int *)XMALLOC(new_bcsize * sizeof(E_Int));
    E_Int *ptr = new_ptlist;

    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      E_Int state = M->faceTree->state(face);
      
      if (state == UNTOUCHED || state == UNREFINED) {
        *ptr++ = face;
      } else if (state == REFINED) {
        Children *children = M->faceTree->children(face);
        for (E_Int k = 0; k < children->n; k++) *ptr++ = children->pc[k];
      }
    }

    XFREE(M->ptlists[i]);

    M->ptlists[i] = new_ptlist;
    M->bcsizes[i] = new_bcsize;
  }

  M->nif = M->nfaces - M->nbf;
}

void update_global_cells_after_ref(AMesh *M)
{
  E_Int nnew_cells = M->ncells - M->prev_ncells;
  
  E_Int first_new_cell;

  MPI_Scan(&nnew_cells, &first_new_cell, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  E_Int gncells;
  MPI_Allreduce(&M->prev_ncells, &gncells, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (nnew_cells > 0) {
    M->gcells = (E_Int *)XRESIZE(M->gcells, M->ncells * sizeof(E_Int));

    for (E_Int i = 0; i < nnew_cells; i++)
      M->gcells[M->prev_ncells + i] = gncells + first_new_cell - nnew_cells + i;
  }
}

void update_global_faces_after_ref(AMesh *M)
{
  E_Int nnew_faces = M->nfaces - M->prev_nfaces;
  E_Int first_new_face;

  MPI_Scan(&nnew_faces, &first_new_face, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  E_Int gnfaces;
  MPI_Allreduce(&M->prev_nfaces, &gnfaces, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (nnew_faces == 0) return;

  M->gfaces = (E_Int *)XRESIZE(M->gfaces, M->nfaces * sizeof(E_Int));

  // Init
  for (E_Int i = M->prev_nfaces; i < M->nfaces; i++)
    M->gfaces[i] = -1;

  // Internal faces first
  E_Int incr = first_new_face - nnew_faces;

  for (E_Int i = 0; i < nnew_faces; i++) {
    E_Int face = M->prev_nfaces + i;
    if (M->neigh[face] == -1) continue;

    M->gfaces[face] = gnfaces + incr++;
    M->FT->insert({M->gfaces[face], face});
  }

  // Boundary faces next
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];
    
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];

      if (M->gfaces[face] == -1) {
        assert(face >= M->prev_nfaces);
        M->gfaces[face] = gnfaces + incr++;
        M->FT->insert({M->gfaces[face], face});
      }
    }
  }

  // Patch faces last as they might be renumbered and we wouldn't want to
  // create holes in global face numbering

  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int face = P->pf[j];

      if (M->gfaces[face] == -1) {
        M->gfaces[face] = gnfaces + incr++;
        M->FT->insert({M->gfaces[face], face});
      }
    }
  }
}

void update_patch_faces_after_ref(AMesh *M)
{
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];
    E_Int new_nf = 0;
    
    for (E_Int j = 0; j < P->nf; j++) {
      E_Int face = P->pf[j];
      E_Int state = M->faceTree->state(face);

      assert(state != NEW);

      if (state == UNTOUCHED)
        new_nf += 1;
      else if (state == REFINED)
        new_nf += M->faceTree->children(face)->n;
      else if (state == UNREFINED)
        new_nf += 1;
    }

    E_Int *new_pf = (E_Int *)XMALLOC(new_nf * sizeof(E_Int));
    E_Int *new_pn = (E_Int *)XMALLOC(new_nf * sizeof(E_Int));
    E_Int *pf = new_pf;
    E_Int *pn = new_pn;

    // Children faces inherit pn of parent, will be updated later

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int face = P->pf[j];
      E_Int nei = P->pn[j];
      E_Int state = M->faceTree->state(face);

      if (state == UNTOUCHED || state == UNREFINED) {
        *pf++ = face;
        *pn++ = nei;
      } else if (state == REFINED) {
        Children *children = M->faceTree->children(face);

        assert(M->faceTree->state(children->pc[0]) == REFINED);
        
        for (E_Int k = 0; k < children->n; k++) {
          *pf++ = children->pc[k];
          *pn++ = nei;
        }
      }
    }

    XFREE(M->patches[i].pf);
    XFREE(M->patches[i].pn);
    M->patches[i].pf = new_pf;
    M->patches[i].pn = new_pn;
    M->patches[i].nf = new_nf;
  }

  if (M->mode_2D) return;

  // Swap first and third children if iso refinement
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    if (M->pid < P->nei) continue;

    for (E_Int j = 0; j < P->nf; ) {
      E_Int face = P->pf[j++];
      E_Int state = M->faceTree->state(face);
      if (state != REFINED) continue;

      assert(M->faceTree->children(face)->n == 4);

      std::swap(P->pf[j], P->pf[j+2]);
      j += 3;
    }
  }
}

void make_dual_graph(AMesh *M)
{
  M->XNEIS = (E_Int *)XCALLOC((M->ncells+1), sizeof(E_Int));
  E_Int *XNEIS = M->XNEIS;

  // Internal faces
  E_Int nif = 0;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    E_Int own = M->owner[i];
    XNEIS[own+1] += 1;
    XNEIS[nei+1] += 1;
    nif += 1;
  }

  // Proc faces
  E_Int npf = 0;

  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    npf += P->nf;

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int face = P->pf[j];
      XNEIS[M->owner[face]+1] += 1;
    }
  }

  for (E_Int i = 0; i < M->ncells; i++) XNEIS[i+1] += XNEIS[i];

  M->CADJ = (E_Int *)XMALLOC((2*nif + npf) * sizeof(E_Int));
  E_Int *CADJ = M->CADJ;

  E_Int *cneis = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  // Internal edges first
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    E_Int own = M->owner[i];

    CADJ[XNEIS[own] + cneis[own]++] = M->gcells[nei];
    CADJ[XNEIS[nei] + cneis[nei]++] = M->gcells[own];
  }

  // Proc edges
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    for (E_Int j = 0; j < P->nf; j++) {
      E_Int face = P->pf[j];
      E_Int own = M->owner[face];
      CADJ[XNEIS[own] + cneis[own]++] = M->patches[i].pn[j];
    }
  }

  XFREE(cneis);
}

static
void compute_cell_center_hexa(E_Int cell, AMesh *M)
{
  E_Float cc[3] = {0., 0., 0.};
  E_Int *pf = get_facets(cell, M->nface, M->indPH);
  E_Int *pn = get_facets(pf[0], M->ngon, M->indPG);
  for (E_Int i = 0; i < 4; i++) {
    cc[0] += M->x[pn[i]];
    cc[1] += M->y[pn[i]];
    cc[2] += M->z[pn[i]];
  }
  pn = get_facets(pf[1], M->ngon, M->indPG);
  for (E_Int i = 0; i < 4; i++) {
    cc[0] += M->x[pn[i]];
    cc[1] += M->y[pn[i]];
    cc[2] += M->z[pn[i]];
  }
  for (E_Int i = 0; i < 3; i++) cc[i] /= 8.0;
  M->ccenter[cell][0] = cc[0];
  M->ccenter[cell][1] = cc[1];
  M->ccenter[cell][2] = cc[2];
}

void compute_ref_cells_centers(AMesh *M, const std::vector<E_Int> &ref_cells)
{
  M->ccenter.clear();

  for (auto cell : ref_cells) {
    switch (M->cellTree->type(cell)) {
      case HEXA:
        compute_cell_center_hexa(cell, M);
        break;
      default:
        fprintf(stderr, "UNIMPLEMENTED");
        assert(0);
        exit(1);
    }
  }
}

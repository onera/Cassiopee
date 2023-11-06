/*    
    Copyright 2013-2023 Onera.

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
#include <stack>

E_Int is_metric_valid(E_Float *M)
{
  E_Float a = M[0];
  E_Float b = M[1];
  E_Float c = M[2];
  E_Float d = M[3];
  E_Float e = M[4];
  E_Float f = M[5];

  E_Float det = a*d*f - (a*e*e + d*c*c + f*b*b) + 2.*b*c*e;

  return (a > 0.) && (a*d > b*b) && (det > 0.);
}

void hessian_to_metric(E_Float *H, mesh *M)
{
  E_Float L[3], v0[3], v1[3], v2[3];
  E_Float L0, L1, L2, *pt;
  E_Int ncells = M->ncells;
  for (E_Int i = 0; i < ncells; i++) {
    pt = &H[6*i];
    eigen(pt, L, v0, v1, v2);
    L0 = sqrt(fabs(L[0])); L1 = sqrt(fabs(L[1])); L2 = sqrt(fabs(L[2]));
    pt[0] = L0*v0[0]*v0[0] + L1*v1[0]*v1[0] + L2*v2[0]*v2[0];
    pt[1] = L0*v0[0]*v0[1] + L1*v1[0]*v1[1] + L2*v2[0]*v2[1];
    pt[2] = L0*v0[0]*v0[2] + L1*v1[0]*v1[2] + L2*v2[0]*v2[2];
    pt[3] = L0*v0[1]*v0[1] + L1*v1[1]*v1[1] + L2*v2[1]*v2[1];
    pt[4] = L0*v0[1]*v0[2] + L1*v1[1]*v1[2] + L2*v2[1]*v2[2];
    pt[5] = L0*v0[2]*v0[2] + L1*v1[2]*v1[2] + L2*v2[2]*v2[2];
  }
}

void compute_principal_vecs(mesh *M, std::vector<pDirs> &Dvec)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int *pf = &M->NFACE[6*i];
    E_Float *p0, *p1;
    auto& dirI = Dvec[i].I;
    auto& dirJ = Dvec[i].J;
    auto& dirK = Dvec[i].K;
    p0 = &M->fc[3*pf[2]];
    p1 = &M->fc[3*pf[3]];
    for (E_Int j = 0; j < 3; j++) dirI[j] = p1[j] - p0[j];
    
    p0 = &M->fc[3*pf[4]];
    p1 = &M->fc[3*pf[5]];
    for (E_Int j = 0; j < 3; j++) dirJ[j] = p1[j] - p0[j];
    
    p0 = &M->fc[3*pf[0]];
    p1 = &M->fc[3*pf[1]];
    for (E_Int j = 0; j < 3; j++) dirK[j] = p1[j] - p0[j];
  }
}

void compute_cells_principal_vecs(mesh *M, const std::vector<E_Int> &cells, std::vector<pDirs> &Dvec)
{
  for (size_t i = 0; i < cells.size(); i++) {
    E_Int cell = cells[i];
    E_Int *pf = &M->NFACE[6*cell];
    E_Float *p0, *p1;
    auto& dirI = Dvec[i].I;
    auto& dirJ = Dvec[i].J;
    auto& dirK = Dvec[i].K;
    p0 = &M->fc[3*pf[2]];
    p1 = &M->fc[3*pf[3]];
    for (E_Int j = 0; j < 3; j++) dirI[j] = p1[j] - p0[j];
    
    p0 = &M->fc[3*pf[4]];
    p1 = &M->fc[3*pf[5]];
    for (E_Int j = 0; j < 3; j++) dirJ[j] = p1[j] - p0[j];
    
    p0 = &M->fc[3*pf[0]];
    p1 = &M->fc[3*pf[1]];
    for (E_Int j = 0; j < 3; j++) dirK[j] = p1[j] - p0[j];
  }
}

#define MINREF -10

static
void compute_ref_data(mesh *M, E_Float *H, const std::vector<pDirs> &Dirs,
  std::vector<std::array<int, 3>> &refData)
{
  compute_face_centers(M);
  if (M->iso_mode == 0) {
    for (E_Int i = 0; i < M->ncells; i++) {
      E_Float *pH = &H[6*i];
      auto &pref = refData[i];
      E_Float L0, L1, L2, dd[3];
      symmat_dot_vec(pH, Dirs[i].I, dd);
      L0 = norm(dd, 3);
      symmat_dot_vec(pH, Dirs[i].J, dd);
      L1 = norm(dd, 3);
      symmat_dot_vec(pH, Dirs[i].K, dd);
      L2 = norm(dd, 3);
      // get max stretch
      E_Float MAX = std::max(std::max(L0, L1), L2);
      if (MAX < M->Tr) continue;	
      E_Float coeff = M->Gmax/MAX;
      if (L0 >= M->Tr) pref[0] = round(coeff*L0);
      if (L1 >= M->Tr) pref[1] = round(coeff*L1);
      if (L2 >= M->Tr) pref[2] = round(coeff*L2);
    }
  } else {
    for (E_Int i = 0; i < M->ncells; i++) {
      E_Float *pH = &H[6*i];
      auto &pref = refData[i];
      E_Float L0, L1, L2, dd[3];
      symmat_dot_vec(pH, Dirs[i].I, dd);
      L0 = norm(dd, 3);
      symmat_dot_vec(pH, Dirs[i].J, dd);
      L1 = norm(dd, 3);
      symmat_dot_vec(pH, Dirs[i].K, dd);
      L2 = norm(dd, 3);
      E_Float MAX = std::max(std::max(L0, L1), L2);
      if (MAX < M->Tr) continue;
      pref[0] = M->Gmax;
      pref[1] = M->Gmax;
      pref[2] = M->Gmax;
    }
  }
}

// Returns 0 if freeze vector applied or not needed
// Returns 1 if fail
static
E_Int apply_freeze_vector(mesh *M, const std::vector<pDirs> &Dirs,
  E_Float *fvec)
{
  if (!fvec) return 0;

  // In order to freeze refinement in direction fvec,
  // all refinement cells should have a dimension that's parallel to fvec
  // Otherwise, inconsistencies might happen
  E_Float TOL = 1e-12;

  // Normalized freeze vector
  E_Float fNorm = norm(fvec, 3);
  if (fabs(fNorm) < TOL) {
    if (M->pid == 0)
      fprintf(stderr, "adaptMesh(): null freeze vector. Not applied.\n");
    return 0;
  }

  E_Float fvec_n[3] = {fvec[0]/fNorm, fvec[1]/fNorm, fvec[2]/fNorm};
  if (M->pid == 0) {
    printf("Normalized freeze vector: %.3f %.3f %.3f\n",
      fvec_n[0], fvec_n[1], fvec_n[2]);
  }

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int *pr = &M->ref_data[3*i];
    if (pr[0] == 0 && pr[1] == 0 && pr[2] == 0) continue;
    E_Float dp;
    
    dp = dot(Dirs[i].I, fvec_n, 3) / norm(Dirs[i].I, 3);
    if (fabs(fabs(dp)-1.0) < TOL) {
      pr[0] = 0;
      continue;
    }

    dp = dot(Dirs[i].J, fvec_n, 3) / norm(Dirs[i].J, 3);
    if (fabs(fabs(dp)-1.0) < TOL) {
      pr[1] = 0;
      continue;
    }

    dp = dot(Dirs[i].K, fvec_n, 3) / norm(Dirs[i].K, 3);
    if (fabs(fabs(dp)-1.0) < TOL) {
      pr[2] = 0;
      continue;
    }
    
    return 1;
  }
  
  return 0;
}

E_Int make_ref_data(mesh *M, E_Float **sols, E_Int csize, E_Float *fvec)
{
  // Allocate final ref data
  M->ref_data = (E_Int *)XMALLOC((3*M->ncells) * sizeof(E_Int));
  for (E_Int i = 0; i < 3*M->ncells; i++) M->ref_data[i] = MINREF;

  // Compute temp data
  std::vector<pDirs> D_(M->ncells);
  compute_principal_vecs(M, D_);

  if (M->sensor == 0) { // metric sensor
    for (E_Int i = 0; i < csize; i++) {
      // compute hessians
      E_Float *H = compute_hessian(M, sols[i]);

      // process hessians
      hessian_to_metric(H, M);
  
      // compute refinement data
      std::vector<std::array<int, 3>> partialRefData(M->ncells);
      compute_ref_data(M, H, D_, partialRefData);

      for (E_Int i = 0; i < M->ncells; i++) {
        E_Int *pr = &M->ref_data[3*i];
        pr[0] = std::max(pr[0], partialRefData[i][0]);
        pr[1] = std::max(pr[1], partialRefData[i][1]);
        pr[2] = std::max(pr[2], partialRefData[i][2]);
      }

      XFREE(H);
    }
  } else if (M->sensor == 1) { // cell sensor
    // Marked cells: consider only the first sol
    E_Float *sol = sols[0];
    for (E_Int i = 0; i < M->ncells; i++) {
      E_Int *pr = &M->ref_data[3*i];
      if (sol[i] == 0) {
        for (E_Int j = 0; j < 3; j++)
          pr[j] = M->Gmax;
      } else {
        for (E_Int j = 0; j < 3; j++)
          pr[j] = 0;
      }
    }
  }

  // Smooth out refinement data
  smooth_ref_data(M);

  // Block refinement in fvec direction if needed
  E_Int ret = apply_freeze_vector(M, D_, fvec);

  return ret;
}

std::vector<E_Int> compute_canon_info(E_Int cell, mesh *M, E_Int *ps)
{
  if (ps[0] != -1)
    return std::vector<E_Int>(0);
  
  E_Int *pn, *pf, reorient, i, local[4], face, i0;
  E_Int nodes[8];
  memset(nodes, -1, 8*sizeof(E_Int));
  
  pf = &M->NFACE[6*cell];

  // bot
  face = pf[0];
  pn = &M->NGON[4*face];
  reorient = get_reorient(face, cell, 1, M);
  order_quad(local, pn, reorient, 0);
  for (i = 0; i < 4; i++) nodes[i] = local[i];
  ps[0] = 0;

  // lft
  face = pf[2];
  pn = &M->NGON[4*face];
  reorient = get_reorient(face, cell, 1, M);
  i0 = get_pos(nodes[0], pn, 4);
  order_quad(local, pn, reorient, i0);
  assert(local[0] == nodes[0]);
  assert(local[1] == nodes[3]);
  nodes[4] = local[3];
  nodes[7] = local[2];
  ps[2] = i0;

  // rgt
  face = pf[3];
  pn = &M->NGON[4*face];
  reorient = get_reorient(face, cell, 0, M);
  i0 = get_pos(nodes[1], pn, 4);
  order_quad(local, pn, reorient, i0);
  assert(local[0] == nodes[1]);
  assert(local[1] == nodes[2]);
  nodes[5] = local[3];
  nodes[6] = local[2];
  ps[3] = i0;

  // top
  face = pf[1];
  pn = &M->NGON[4*face];
  reorient = get_reorient(face, cell, 0, M);
  i0 = get_pos(nodes[4], pn, 4);
  order_quad(local, pn, reorient, i0);
  assert(local[0] == nodes[4]);
  assert(local[1] == nodes[5]);
  assert(local[2] == nodes[6]);
  assert(local[3] == nodes[7]);
  ps[1] = i0;

  // fro
  face = pf[4];
  pn = &M->NGON[4*face];
  reorient = get_reorient(face, cell, 1, M);
  i0 = get_pos(nodes[1], pn, 4);
  order_quad(local, pn, reorient, i0);
  assert(local[0] == nodes[1]);
  assert(local[1] == nodes[0]);
  assert(local[2] == nodes[4]);
  assert(local[3] == nodes[5]);
  ps[4] = i0;

  // bck
  face = pf[5];
  pn = &M->NGON[4*face];
  reorient = get_reorient(face, cell, 0, M);
  i0 = get_pos(nodes[2], pn, 4);
  order_quad(local, pn, reorient, i0);
  assert(local[0] == nodes[2]);
  assert(local[1] == nodes[3]);
  assert(local[2] == nodes[7]);
  assert(local[3] == nodes[6]);
  ps[5] = i0;

  std::vector<E_Int> res(8);
  for (E_Int i = 0; i < 8; i++)
    res[i] = nodes[i];
  return res;
}

static
void bot_to_bot(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int b = cell_dir[1];

  if (!reorient) {
    switch (i0) {
      case 1:
      case 3:
        cell_dir[0] = b; cell_dir[1] = a;
        break;
      default:
        break;
    }
  } else {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = b; cell_dir[1] = a;
        break;
      default:
        break;
    }
  }
}

static
void lft_to_bot(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int b = cell_dir[1];
  E_Int c = cell_dir[2];

  cell_dir[2] = a;

  if (!reorient) {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = b; cell_dir[1] = c;
        break;
      case 1:
      case 3:
        cell_dir[0] = c; cell_dir[1] = b;
        break;
      default:
        break;
    }
  } else {
    switch (i0) {
      case 1:
      case 3:
        cell_dir[0] = b; cell_dir[1] = c;
        break;
      case 0:
      case 2:
        cell_dir[0] = c; cell_dir[1] = b;
        break;
      default:
        break;
    }
  }
}

static
void fro_to_bot(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int b = cell_dir[1];
  E_Int c = cell_dir[2];

  cell_dir[2] = b;

  if (!reorient) {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = a; cell_dir[1] = c;
        break;
      case 1:
      case 3:
        cell_dir[0] = c; cell_dir[1] = a;
        break;
      default:
        break;
    }
  } else {
    switch (i0) {
      case 1:
      case 3:
        cell_dir[0] = a; cell_dir[1] = c;
        break;
      case 0:
      case 2:
        cell_dir[0] = c; cell_dir[1] = a;
        break;
      default:
        break;
    }
  }
}

static
void dir_to_bot(E_Int pos_face, E_Int i0, E_Int reorient, E_Int *ref_dir)
{
  switch (pos_face) {
    case 0:
    case 1:
      bot_to_bot(i0, reorient, ref_dir);
      break;
    case 2:
    case 3:
      lft_to_bot(i0, reorient, ref_dir);
      break;
    case 4:
    case 5:
      fro_to_bot(i0, reorient, ref_dir);
      break;
    default:
      assert(0);
      break;
  }
}

static
void bot_to_lft(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int b = cell_dir[1];
  E_Int c = cell_dir[2];

  cell_dir[0] = c;

  if (!reorient) {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[1] = a; cell_dir[2] = b;
        break;
      case 1:
      case 3:
        cell_dir[1] = b; cell_dir[2] = a;
        break;
      default:
        assert(0);
        break;
    }
  } else {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[1] = b; cell_dir[2] = a;
        break;
      case 1:
      case 3:
        cell_dir[1] = a; cell_dir[2] = b;
        break;
      default:
        assert(0);
        break;
    }
  }
}

static
void lft_to_lft(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int b = cell_dir[1];
  E_Int c = cell_dir[2];

  if (!reorient) {
    switch (i0) {
      case 1:
      case 3:
        cell_dir[1] = c; cell_dir[2] = b;
        break;
      default:
        break;
    }
  } else {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[1] = c; cell_dir[2] = b;
        break;
      default:
        break;
    }
  }
}

static
void fro_to_lft(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int b = cell_dir[1];
  E_Int c = cell_dir[2];

  cell_dir[0] = b;

  if (!reorient) {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[1] = a; cell_dir[2] = c;
        break;
      case 1:
      case 3:
        cell_dir[1] = c; cell_dir[2] = a;
        break;
      default:
        assert(0);
        break;
    }
  } else {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[1] = c; cell_dir[2] = a;
        break;
      case 1:
      case 3:
        cell_dir[1] = a; cell_dir[2] = c;
        break;
      default:
        assert(0);
        break;
    }
  }
}

static
void dir_to_lft(E_Int pos_face, E_Int i0, E_Int reorient, E_Int *ref_dir)
{
  switch (pos_face) {
    case 0:
    case 1:
      bot_to_lft(i0, reorient, ref_dir);
      break;
    case 2:
    case 3:
      lft_to_lft(i0, reorient, ref_dir);
      break;
    case 4:
    case 5:
      fro_to_lft(i0, reorient, ref_dir);
      break;
    default:
      assert(0);
      break;
  }
}

static
void bot_to_fro(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int b = cell_dir[1];
  E_Int c = cell_dir[2];

  cell_dir[1] = c;

  if (!reorient) {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = a; cell_dir[2] = b;
        break;
      case 1:
      case 3:
        cell_dir[0] = b; cell_dir[2] = a;
        break;
      default:
        assert(0);
        break;
    }
  } else {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = b; cell_dir[2] = a;
        break;
      case 1:
      case 3:
        cell_dir[0] = a; cell_dir[2] = b;
        break;
      default:
        assert(0);
        break;
    }
  }
}

static
void lft_to_fro(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int b = cell_dir[1];
  E_Int c = cell_dir[2];

  cell_dir[1] = a;

  if (!reorient) {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = b; cell_dir[2] = c;
        break;
      case 1:
      case 3:
        cell_dir[0] = c; cell_dir[2] = b;
        break;
      default:
        assert(0);
        break;
    }
  } else {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = c; cell_dir[2] = b;
        break;
      case 1:
      case 3:
        cell_dir[0] = b; cell_dir[2] = c;
        break;
      default:
        assert(0);
        break;
    }
  }
}

static
void fro_to_fro(E_Int i0, E_Int reorient, E_Int *cell_dir)
{
  E_Int a = cell_dir[0];
  E_Int c = cell_dir[2];

  if (!reorient) {
    switch (i0) {
      case 1:
      case 3:
        cell_dir[0] = c; cell_dir[2] = a;
        break;
      default:
        break;
    }
  } else {
    switch (i0) {
      case 0:
      case 2:
        cell_dir[0] = c; cell_dir[2] = a;
        break;
      default:
        break;
    }
  }
}

static
void dir_to_fro(E_Int pos_face, E_Int i0, E_Int reorient, E_Int *ref_dir)
{
  switch (pos_face) {
    case 0:
    case 1:
      bot_to_fro(i0, reorient, ref_dir);
      break;
    case 2:
    case 3:
      lft_to_fro(i0, reorient, ref_dir);
      break;
    case 4:
    case 5:
      fro_to_fro(i0, reorient, ref_dir);
      break;
    default:
      assert(0);
      break;
  }
}

void deduce_nei_ref_data(E_Int face_idx, E_Int pos_face_in_nei,
  E_Int start_node, E_Int reorient, E_Int *local_dir)
{
  assert(pos_face_in_nei >= 0 && pos_face_in_nei < 6);
  assert(start_node >= 0 && start_node < 4);

  switch (face_idx) {
    case 0:
    case 1:
      dir_to_bot(pos_face_in_nei, start_node, reorient, local_dir);
      break;
    case 2:
    case 3:
      dir_to_lft(pos_face_in_nei, start_node, reorient, local_dir);
      break;
    case 4:
    case 5:
      dir_to_fro(pos_face_in_nei, start_node, reorient, local_dir);
      break;
    default:
      assert(0);
      break;
  }
}

void fix_conflict(E_Int *ac, E_Int *bc, E_Int *an, E_Int *bn)
{
  E_Int incr[4] = {*ac, *bc, *an, *bn};
  E_Int INCR = 0;
  E_Int i, MIN;
  
  while (1) {
    // one cell is done
    if ((incr[0] == 0 && incr[1] == 0) || (incr[2] == 0 && incr[3] == 0)) break;
    // first direction is done
    else if (incr[0] == 0 && incr[2] == 0) break;
    // second direction is done
    else if (incr[1] == 0 && incr[3] == 0) break;

    MIN = 1000;
    for (i = 0; i < 4; i++) {
      if (incr[i] == 0) continue;
      if (incr[i] < MIN) MIN = incr[i];
    }
    assert(MIN != 1000);

    INCR += MIN;

    // decrement
    for (i = 0; i < 4; i++) {
      incr[i] -= MIN;
      if (incr[i] < 0) incr[i] = 0;
    }
  }

  *ac = std::max(*ac, INCR);
  *bc = std::max(*bc, INCR);
  *an = std::max(*an, INCR);
  *bn = std::max(*bn, INCR);
}

void smooth_ref_data(mesh *M)
{
  std::vector<E_Int> start_nodes(6*M->ncells, -1);

  E_Int max_exchanges = 10;
  E_Int exchange = 0;
  
  E_Int i, j, k, cell, ncells, *pf, *pn, *pr, *ps, *prn, face, ac, bc, gc, an,
    bn, gn, has_changed, reorient_c, reorient_n, nei, fpos, i0, ac_orig,
    bc_orig, gc_orig, An_orig[6][3], an_orig, bn_orig, gn_orig, npfaces,
    *pfaces, *owner;
  
  ncells = M->ncells;
  owner = M->owner;

  // first, exchange topological info (doesn't change during smoothing)
  E_Int **pnei_topo_data = (E_Int **)XCALLOC(M->nppatches, sizeof(E_Int *));

  E_Int dest, l;

  E_Int **send_buf = (E_Int **)XCALLOC(M->nppatches, sizeof(E_Int *));

  for (i = 0; i < M->nppatches; i++) {
    npfaces = M->ppatches[i].nfaces;
    pfaces = M->ppatches[i].faces;
    dest = M->ppatches[i].nei_proc;
  
    // fpos + reorient_n + start_node_nei_of_pface = 3
    pnei_topo_data[i] = (E_Int *)XCALLOC(3*npfaces, sizeof(E_Int));

    send_buf[i] = (E_Int *)XCALLOC(3*npfaces, sizeof(E_Int));

    l = 0;
    for (j = 0; j < npfaces; j++) {
      cell = owner[pfaces[j]];
      pf = &M->NFACE[6*cell];
      ps = &start_nodes[6*cell];
      // fpos
      fpos = get_pos(pfaces[j], pf, 6);
      assert(fpos != -1);
      // reorient_n
      reorient_n = get_reorient(pfaces[j], cell, normalIn[fpos], M);
      // start_node_nei
      compute_canon_info(cell, M, ps);

      send_buf[i][l++] = fpos;
      send_buf[i][l++] = reorient_n;
      send_buf[i][l++] = ps[fpos];
    }
    assert(l == 3*npfaces);

    MPI_Irecv(pnei_topo_data[i], (int)3*npfaces, XMPI_INT, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);
    MPI_Isend(send_buf[i], (int)3*npfaces, XMPI_INT, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);

    assert(M->nreq < 2*M->npc);
  }

  comm_waitall(M);

  // allocate proc ref data 
  E_Int **pnei_ref_data = (E_Int **)XCALLOC(M->nppatches, sizeof(E_Int *));
  for (i = 0; i < M->nppatches; i++) {
    pnei_ref_data[i] = (E_Int *)XCALLOC(3*M->ppatches[i].nfaces, sizeof(E_Int));
  }

  while (++exchange <= max_exchanges) {
    // init refinement stack
    std::stack<E_Int> stk;
    for (i = 0; i < ncells; i++) {
      if (is_cell_to_refine(&M->ref_data[3*i]))
        stk.push(i);
    }

    // smooth locally
    while (!stk.empty()) {
      cell = stk.top();
      stk.pop();

      assert(cell >= 0 && cell < M->ncells);

      ps = &start_nodes[6*cell];
      compute_canon_info(cell, M, ps);

      pf = &M->NFACE[6*cell];

      pr = &M->ref_data[3*cell];
      ac = pr[0];
      bc = pr[1];
      gc = pr[2];

      // deduce ref_data for all neighbours
      for (i = 0; i < 6; i++) {
        face = pf[i];
        reorient_c = get_reorient(face, cell, normalIn[i], M);
        nei = get_neighbour(cell, face, M);
        if (nei == -1) continue;
        assert(nei < M->ncells);
        compute_canon_info(nei, M, &start_nodes[6*nei]);
        pn = &M->NFACE[6*nei];
        fpos = get_pos(face, pn, 6);
        assert(fpos != -1);
        reorient_n = get_reorient(face, nei, normalIn[fpos], M);
        i0 = abs(start_nodes[6*nei+fpos] - ps[i]);
        prn = &M->ref_data[3*nei];
        for (j = 0; j < 3; j++) An_orig[i][j] = prn[j];
        deduce_nei_ref_data(i, fpos, i0, reorient_c != reorient_n, prn);
      }

      do {
        has_changed = 0;

        for (i = 0; i < 6; i++) {
          ac_orig = ac;
          bc_orig = bc;
          gc_orig = gc;

          nei = get_neighbour(cell, pf[i], M);
          if (nei == -1)
            continue;

          prn = &M->ref_data[3*nei];
          an = prn[0];
          bn = prn[1];
          gn = prn[2];

          an_orig = an;
          bn_orig = bn;
          gn_orig = gn;

          // alpha
          if (an > ac + 1) ac = an - 1;
          else if (ac > an + 1) an = ac - 1;
          // beta
          if (bn > bc + 1) bc = bn - 1;
          else if (bc > bn + 1) bn = bc - 1;
          // gamma
          if (gn > gc + 1) gc = gn - 1;
          else if (gc > gn + 1) gn = gc - 1;

          if (i == 0 || i == 1)
            fix_conflict(&ac, &bc, &an, &bn);
          else if (i == 2 || i == 3)
            fix_conflict(&bc, &gc, &bn, &gn);
          else if (i == 4 || i == 5)
            fix_conflict(&ac, &gc, &an, &gn);

          pr[0] = ac; pr[1] = bc; pr[2] = gc;
          prn[0] = an; prn[1] = bn; prn[2] = gn;

          has_changed |=
            !((ac == ac_orig) && (bc == bc_orig) && (gc == gc_orig));
          has_changed |=
            !((an == an_orig) && (bn == bn_orig) && (gn == gn_orig));
        }
      } while (has_changed);

      // reverse
      for (i = 0; i < 6; i++) {
        face = pf[i];
        reorient_c = get_reorient(face, cell, normalIn[i], M);
        nei = get_neighbour(cell, face, M);
        if (nei == -1) continue;
        pn = &M->NFACE[6*nei];
        fpos = get_pos(face, pn, 6);
        assert(fpos != -1);
        reorient_n = get_reorient(face, nei, normalIn[fpos], M);
        i0 = abs(start_nodes[6*nei+fpos] - ps[i]);
        prn = &M->ref_data[3*nei];

        deduce_nei_ref_data(fpos, i, i0, reorient_c != reorient_n, prn);

        for (j = 0; j < 3; j++) {
          if (prn[j] != An_orig[i][j]) {
            stk.push(nei);
            break;
          }
        }
      }
    }

    // exchange proc ref data
    MPI_Barrier(MPI_COMM_WORLD); 
    comm_interface_data_i(M, M->ref_data, 3, pnei_ref_data);

    // address proc cells
    E_Int lstop = 0;
    for (E_Int i = 0; i < M->nppatches; i++) {
      E_Int *nei_ref_data = pnei_ref_data[i];
      E_Int *nei_topo_data = pnei_topo_data[i];
      for (E_Int j = 0; j < M->ppatches[i].nfaces; j++) {
        E_Int face = M->ppatches[i].faces[j];
        E_Int cell = M->owner[face];
        pf = &M->NFACE[6*cell];
        ps = &start_nodes[6*cell];
        compute_canon_info(cell, M, ps);
        k = get_pos(face, pf, 6);
        assert(k != -1);
        reorient_c = get_reorient(face, cell, normalIn[k], M);

        // need fpos, reorient_n, start_node_nei, ref_data_nei
        // reorient_n has to be 0 for nei
        prn = &nei_ref_data[3*j];
        E_Int *pnei_topo = &nei_topo_data[3*j];
        fpos = pnei_topo[0];
        assert(fpos >= 0 && fpos < 6);
        reorient_n = pnei_topo[1];
        i0 = abs(pnei_topo[2] - ps[k]);
        deduce_nei_ref_data(k, fpos, i0, reorient_c == reorient_n, prn);

        pr = &M->ref_data[3*cell];
        ac = pr[0];
        bc = pr[1];
        gc = pr[2];
        ac_orig = ac;
        bc_orig = bc;
        gc_orig = gc;
        an = prn[0];
        bn = prn[1];
        gn = prn[2];

        // alpha
        if (an > ac + 1) ac = an - 1;
        // beta
        if (bn > bc + 1) bc = bn - 1;
        // gamma
        if (gn > gc + 1) gc = gn - 1;

        if (k == 0 || k == 1)
          fix_conflict(&ac, &bc, &an, &bn);
        else if (k == 2 || k == 3)
          fix_conflict(&bc, &gc, &bn, &gn);
        else if (k == 4 || k == 5)
          fix_conflict(&ac, &gc, &an, &gn);

        pr[0] = ac;
        pr[1] = bc;
        pr[2] = gc;

        E_Int pcell_changed =
          (ac != ac_orig) || (bc != bc_orig) || (gc != gc_orig);
        if (pcell_changed)
          lstop = 1;
      }
    }

    // stop condition is nothing to smooth anywhere
    E_Int gstop;
    MPI_Allreduce(&lstop, &gstop, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if (gstop == 0)
      break;
  }
  if (exchange > max_exchanges)
    fprintf(stderr, "smooth_ref_data(): Warning: exceeded max_exchanges");

  for (i = 0; i < M->nppatches; i++) {
    XFREE(pnei_topo_data[i]);
    XFREE(pnei_ref_data[i]);
    XFREE(send_buf[i]);
  }
  XFREE(pnei_topo_data);
  XFREE(pnei_ref_data);
  XFREE(send_buf);
}

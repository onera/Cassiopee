/*    
    Copyright 2013-2025 Onera.

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
#include "metric.h"
#include "Connect/connect.h"
#include "Math/math.h"
#include "String/kstring.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <stack>
#include <set>

#define INTERNAL 0
#define EXTERNAL 1

void K_METRIC::compute_face_center_and_area(
  E_Int id, E_Int stride,
  E_Int *pn, E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa
)
{
  // Init
  fa[0] = fa[1] = fa[2] = 0.0;
  fc[0] = fc[1] = fc[2] = 0.0;

  // Approximate face center
  E_Float fcenter[3] = {0,0,0};

  for (E_Int i = 0; i < stride; i++) {
    E_Int point = pn[i]-1;
    fcenter[0] += x[point];
    fcenter[1] += y[point];
    fcenter[2] += z[point];
  }

  for (E_Int i = 0; i < 3; i++) fcenter[i] /= stride;

  // Sum of triangle area vectors
  E_Float sumN[3] = {0,0,0};
  // Sum of triangle areas
  E_Float sumA = 0;
  // Sum of area-weighted triangle centers
  E_Float sumAc[3] = {0,0,0};

  // Compute area vector and center of stride-2 triangles
  // formed by p0p1p2, p0p2p3, ... p0p(stride-2)p(stride-1)
  E_Int p0 = pn[0]-1;
  for (E_Int i = 1; i < stride-1; i++)
  {
    E_Int p1 = pn[i]-1;
    E_Int p2 = pn[i+1]-1;

    // Triangle center
    E_Float tc[3];
    tc[0] = x[p0] + x[p1] + x[p2];
    tc[1] = y[p0] + y[p1] + y[p2];
    tc[2] = z[p0] + z[p1] + z[p2];

    // Area vector
    E_Float n[3];
    E_Float v10[3] = {x[p1]-x[p0], y[p1]-y[p0], z[p1]-z[p0]};
    E_Float v20[3] = {x[p2]-x[p0], y[p2]-y[p0], z[p2]-z[p0]};
    K_MATH::cross(v10, v20, n);

    // Area
    E_Float a = K_MATH::norm(n, 3);

    for (E_Int j = 0; j < 3; j++) {
      sumN[j] += n[j];
      sumAc[j] += a*tc[j];
    }
    sumA += a;
  }

  // Deal with zero-area faces
  if (sumA < K_MATH::SMALL)
  {
    fprintf(stderr, "compute_face_area_and_center(): "
      "Warning: Face: " SF_D_ " - Area: " SF_F_ " - Tol: %.2e\n", id, sumA, K_MATH::SMALL);
    for (E_Int i = 0; i < 3; i++) {
      fc[i] = fcenter[i];
      fa[i] = 0.0;
    }
  }
  else
  {
    for (E_Int i = 0; i < 3; i++) {
      fc[i] = sumAc[i]/(3.0*sumA);
      fa[i] = 0.5*sumN[i];
    }
  }
}

// Assumes well-oriented mesh
void K_METRIC::compute_cell_center_and_volume(E_Int id, E_Int stride,
  E_Int *pf, E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa,
  E_Int *owner, E_Float &cx, E_Float &cy, E_Float &cz, E_Float &vol)
{
  // Estimate cell centers as average of face centers
  E_Float cEst[3] = {0.0, 0.0, 0.0};

  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i]-1;
    E_Float *fcenter = &fc[3*face];
    
    for (E_Int j = 0; j < 3; j++)
      cEst[j] += fcenter[j];
  }
  for (E_Int i = 0; i < 3; i++) cEst[i] /= stride;

  // Init
  vol = 0.0;
  cx = cy = cz = 0.0;

  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i]-1;
    E_Float *fcenter = &fc[3*face];
    E_Float *farea = &fa[3*face];

    E_Float d[3] = {fcenter[0]-cEst[0], fcenter[1]-cEst[1], fcenter[2]-cEst[2]};
    
    // Tree times the face-pyramid volume
    E_Float pyr3Vol = K_MATH::dot(farea, d, 3);

    // Negate it if current cell is not the owner of face
    if (id != owner[face]-1) pyr3Vol = -pyr3Vol;

    // Face-pyramid center
    E_Float pyrCenter[3];
    for (E_Int j = 0; j < 3; j++)
      pyrCenter[j] = 0.75*fcenter[j] + 0.25*cEst[j];
    
    // Accumulate cell center and volume
    cx += pyr3Vol * pyrCenter[0];
    cy += pyr3Vol * pyrCenter[1];
    cz += pyr3Vol * pyrCenter[2];
    vol += pyr3Vol;
  }

  if (fabs(vol) >= K_MATH::SMALL)
  {
    cx /= vol;
    cy /= vol;
    cz /= vol;
  }
  else
  {
    cx = cEst[0];
    cy = cEst[1];
    cz = cEst[2];
  }
  
  vol *= K_CONST::ONE_THIRD;

  if (vol < 0.0)
    fprintf(stderr, "Warning: cell " SF_D_ " has negative volume %.4e\n", id, vol);
}

// Assumes closed cell
// Does not assume faces pointing outwards: volume might be negative
void K_METRIC::compute_cell_volume(E_Int cell, K_FLD::FldArrayI &cn, E_Float *x, E_Float *y,
  E_Float *z, E_Float &vol, E_Int refIdx)
{
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();

  // Orient the faces coherently
  std::vector<E_Int> NGON;
  std::vector<E_Int> INDPG(1, 0);
  E_Int stride = -1;
  E_Int *pf = cn.getElt(cell, stride, nface, indPH); 
  
  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i]-1;
    E_Int np = -1;
    E_Int *pn = cn.getFace(face, np, ngon, indPG);
    INDPG.push_back(np);
    for (E_Int j = 0; j < np; j++)
      NGON.push_back(pn[j]);
  }

  for (E_Int i = 0; i < stride; i++)
    INDPG[i+1] += INDPG[i];

  // Fix orientation of first face
  std::vector<E_Int> orient(stride);
  orient[refIdx] = 1;
  std::vector<E_Int> neis(NGON.size());
  K_CONNECT::build_face_neighbourhood(NGON, INDPG, neis);
  K_CONNECT::reversi_connex(&NGON[0], &INDPG[0], stride, &neis[0], refIdx, orient);

  // Apply orientation in local NGON
  for (E_Int i = 0; i < stride; i++) {
    if (orient[i] == -1) {
      E_Int start = INDPG[i];
      E_Int np = INDPG[i+1] - start;
      E_Int *pn = &NGON[start];
      std::reverse(pn+1, pn+np);
    }
  }

  // Compute faces area and center
  std::vector<E_Float> faceAreas(3*stride, 0.0);
  std::vector<E_Float> faceCenters(3*stride, 0.0);

  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i]-1;
    E_Float *fa = &faceAreas[3*i];
    E_Float *fc = &faceCenters[3*i];
    E_Int np = INDPG[i+1]-INDPG[i];
    E_Int *pn = &NGON[INDPG[i]];
    K_METRIC::compute_face_center_and_area(face, np, pn, x, y, z, fc, fa);
  }

  // Estimate cell centroid as average of face centers
  E_Float cc[3] = {0,0,0};
  for (E_Int i = 0; i < stride; i++) {
    E_Float *fc = &faceCenters[3*i];
    for (E_Int j = 0; j < 3; j++)
      cc[j] += fc[j];
  }
  for (E_Int i = 0; i < 3; i++)
    cc[i] /= stride;

  // Compute cell volume
  vol = 0.0;

  for (E_Int i = 0; i < stride; i++) {
    E_Float *fa = &faceAreas[3*i];
    E_Float *fc = &faceCenters[3*i];
    
    // Compute 3*face-pyramid volume contribution
    E_Float d[3] = {fc[0]-cc[0], fc[1]-cc[1], fc[2]-cc[2]};
    //E_Float pyr3Vol = K_MATH::dot(fa, fc, 3);
    E_Float pyr3Vol = K_MATH::dot(fa, d, 3);

    vol += pyr3Vol;
  }

  vol *= K_CONST::ONE_THIRD;
}

// Assumes external faces oriented outwards + parent elements computed
// Volume should always be positive
static
void compute_volumes(K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  std::vector<E_Int> &owner, std::vector<E_Int> &neigh, E_Float *vols)
{
  E_Int *nface = cn.getNFace();
  E_Int *ngon = cn.getNGon();
  E_Int *indPH = cn.getIndPH();
  E_Int *indPG = cn.getIndPG();
  E_Int ncells = cn.getNElts();
  E_Int nfaces = cn.getNFaces();

  // Compute face areas and centers
  std::vector<E_Float> fA(3*nfaces, 0);
  std::vector<E_Float> fC(3*nfaces, 0);
  
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = cn.getFace(i, np, ngon, indPG);
    K_METRIC::compute_face_center_and_area(i, np, pn, x, y, z,
      &fC[3*i], &fA[3*i]);
  }

  // Estimate cell centers as average of face centers
  std::vector<E_Float> cC(3*ncells, 0);
  for (E_Int i = 0; i < nfaces; i++) {
    E_Float *fc = &fC[3*i];
    E_Float *cc;
    
    E_Int own = owner[i];
    cc = &cC[3*own];

    for (E_Int j = 0; j < 3; j++) cc[j] += fc[j];

    E_Int nei = neigh[i];
    if (nei == -1) continue;
    
    cc = &cC[3*nei];

    for (E_Int j = 0; j < 3; j++) cc[j] += fc[j];
  }

  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    cn.getElt(i, stride, nface, indPH);
    E_Float *cc = &cC[3*i];
    for (E_Int j = 0; j < 3; j++)
      cc[j] /= stride;
  }

  for (E_Int i = 0; i < ncells; i++) vols[i] = 0.0;

  for (E_Int i = 0; i < nfaces; i++) {
    E_Float *fa = &fA[3*i];
    E_Float *fc = &fC[3*i];
    E_Float pyr3vol, *cc, d[3];

    E_Int own = owner[i];
    cc = &cC[3*own];
    for (E_Int j = 0; j < 3; j++) d[j] = fc[j]-cc[j];
    pyr3vol = K_MATH::dot(fa, d, 3);
    vols[own] += pyr3vol;

    E_Int nei = neigh[i];
    if (nei == -1) continue;

    cc = &cC[3*nei];
    for (E_Int j = 0; j < 3; j++) d[j] = cc[j]-fc[j];
    pyr3vol = K_MATH::dot(fa, d, 3);
    vols[nei] += pyr3vol;
  }

  for (E_Int i = 0; i < ncells; i++) {
    vols[i] *= K_CONST::ONE_THIRD;
    if (vols[i] < 0.0)
      fprintf(stderr, "Warning: cell " SF_D_ " has negative volume %.4e\n", i, vols[i]);
  }
}

// Assumes unsigned faces
// Returns 0 if ok
// Returns 2 if some cells are not closed
E_Int K_METRIC::compute_volumes_ngon(E_Float *x, E_Float *y, E_Float *z,
  K_FLD::FldArrayI &cn, E_Float *vols)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();
  E_Int *nface = cn.getNFace();

  E_Int ret = 0;
  
  // Make unsigned NFACE
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = cn.getElt(i, stride, nface, indPH);
    for (E_Int j = 0; j < stride; j++)
      if (pf[j] < 0) pf[j] = -pf[j];
  }

  // Check mesh
  std::vector<E_Int> is_cell_open(ncells, 0);
  E_Int bad_faces = K_CONNECT::check_overlapping_cells(cn);
  E_Int bad_cells = K_CONNECT::check_open_cells(cn, &is_cell_open[0]);

  if (!bad_cells && !bad_faces)
  {
    K_CONNECT::orient_boundary_ngon(x, y, z, cn);
    E_Int nfaces = cn.getNFaces();
    std::vector<E_Int> owner(nfaces), neigh(nfaces);
    K_CONNECT::build_parent_elements_ngon(cn, &owner[0], &neigh[0]);
    compute_volumes(cn, x, y, z, owner, neigh, vols);
  }
  else
  {
    // Compute closed cells volumes
    for (E_Int i = 0; i < ncells; i++) {
      if (is_cell_open[i]) {
        ret = 2;
        fprintf(stderr,
          "Warning: Cell " SF_D_ " is not closed. Setting its volume to zero\n", i);
        continue;
      }
      K_METRIC::compute_cell_volume(i, cn, x, y, z, vols[i]);
      if (vols[i] < 0.0) vols[i] = -vols[i];
    }
    ret = 2;
  }

  return ret;
}

void K_METRIC::compute_face_centers_and_areas(
  K_FLD::FldArrayI &cn, E_Float *x,
  E_Float *y, E_Float *z, E_Float *fcenters, E_Float *fareas
)
{
  E_Int nfaces = cn.getNFaces();
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = cn.getFace(i, np, ngon, indPG);
    K_METRIC::compute_face_center_and_area(i, np, pn, x, y, z, &fcenters[3*i],
      &fareas[3*i]);
  }
}

void K_METRIC::compute_cell_centers_and_vols(
  K_FLD::FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
  E_Int *owner, E_Int *neigh, E_Float *fcenters, E_Float *fareas,
  E_Float *cx, E_Float *cy, E_Float *cz, E_Float *volumes
)
{
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int ncells = cn.getNElts();
  E_Int nfaces = cn.getNFaces();
  
  E_Float *vols;
  if (volumes) vols = volumes;
  else vols = (E_Float *)malloc(ncells * sizeof(E_Float)); 

  // Estimate cell centers as average of face centers
  std::vector<E_Float> cEst(3*ncells, 0);
  for (E_Int i = 0; i < nfaces; i++) {
    E_Float *fc = &fcenters[3*i];
    E_Float *cE;
    
    E_Int own = owner[i];
    assert(own >= 0 && own < ncells);
    cE = &cEst[3*own];

    for (E_Int j = 0; j < 3; j++) cE[j] += fc[j];

    E_Int nei = neigh[i];
    assert(nei >= -1 && nei < ncells);
    if (nei == -1) continue;
    
    cE = &cEst[3*nei];

    for (E_Int j = 0; j < 3; j++) cE[j] += fc[j];
  }

  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    cn.getElt(i, stride, nface, indPH);
    E_Float *cE = &cEst[3*i];
    for (E_Int j = 0; j < 3; j++)
      cE[j] /= stride;
  }

  memset(vols, 0, ncells*sizeof(E_Float));
  memset(cx,   0, ncells*sizeof(E_Float));
  memset(cy,   0, ncells*sizeof(E_Float));
  memset(cz,   0, ncells*sizeof(E_Float));

  for (E_Int i = 0; i < nfaces; i++) {
    E_Float *fa = &fareas[3*i];
    E_Float *fc = &fcenters[3*i];
    E_Float pyr3vol, pc[3], *cE, d[3];

    E_Int own = owner[i];
    assert(own >= 0 && own < ncells);
    cE = &cEst[3*own];
    for (E_Int j = 0; j < 3; j++) d[j] = fc[j]-cE[j];
    pyr3vol = K_MATH::dot(fa, d, 3);
    for (E_Int j = 0; j < 3; j++) pc[j]= 0.75*fc[j] + 0.25*cE[j];
    cx[own] += pyr3vol*pc[0];
    cy[own] += pyr3vol*pc[1];
    cz[own] += pyr3vol*pc[2];
    vols[own] += pyr3vol;

    E_Int nei = neigh[i];
    assert(nei >= -1 && nei < ncells);
    if (nei == -1) continue;

    cE = &cEst[3*nei];
    for (E_Int j = 0; j < 3; j++) d[j] = cE[j]-fc[j];
    pyr3vol = K_MATH::dot(fa, d, 3);
    for (E_Int j = 0; j < 3; j++) pc[j]= 0.75*fc[j] + 0.25*cE[j];
    cx[nei] += pyr3vol*pc[0];
    cy[nei] += pyr3vol*pc[1];
    cz[nei] += pyr3vol*pc[2];
    vols[nei] += pyr3vol;
  }

  for (E_Int i = 0; i < ncells; i++) {
    E_Float coeff = 1.0/vols[i];
    cx[i] *= coeff;
    cy[i] *= coeff;
    cz[i] *= coeff;
    vols[i] /= 3.0;
    if (vols[i] <= 0.0) {
      fprintf(stderr, "Warning: cell " SF_D_ " has negative volume " SF_F_ "\n",
        i, vols[i]);
    }
    assert(vols[i] > 0.0);
  }

  if (!volumes) free(vols);
}

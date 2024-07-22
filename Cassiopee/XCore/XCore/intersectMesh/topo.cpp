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
#include <stack>

#define INTERNAL 0
#define EXTERNAL 1

#define DSMALL 1e-14

static
Int *mesh_get_face(Int face, Int &stride, IMesh *M)
{
    stride = M->fpoints[face].size();
    return &(M->fpoints[face][0]);
}

static
Int *mesh_get_cell(Int cell, Int &stride, IMesh *M)
{
    stride = M->cfaces[cell].size();
    return &(M->cfaces[cell][0]);
}

static
void Compute_cell_volume(Int cell, IMesh *M, Float *x, Float *y,
  Float *z, Float &vol, Int refIdx)
{
  // Orient the faces coherently
  std::vector<Int> NGON;
  std::vector<Int> INDPG(1, 0);
  Int stride = -1;
  Int *pf = mesh_get_cell(cell, stride, M);

  for (Int i = 0; i < stride; i++) {
    Int face = pf[i];
    Int np = -1;
    Int *pn = mesh_get_face(face, np, M);
    INDPG.push_back(np);
    for (Int j = 0; j < np; j++)
      NGON.push_back(pn[j]);
  }

  for (Int i = 0; i < stride; i++)
    INDPG[i+1] += INDPG[i];

  // Fix orientation of first face
  std::vector<Int> orient(stride);
  orient[refIdx] = 1;
  std::vector<Int> neis(NGON.size());
  K_CONNECT::build_face_neighbourhood(NGON, INDPG, neis);
  K_CONNECT::reversi_connex(&NGON[0], &INDPG[0], stride, &neis[0], refIdx, orient);

  // Apply orientation in local NGON
  for (Int i = 0; i < stride; i++) {
    if (orient[i] == -1) {
      Int start = INDPG[i];
      Int np = INDPG[i+1] - start;
      Int *pn = &NGON[start];
      std::reverse(pn+1, pn+np);
    }
  }

  // Compute faces area and center
  std::vector<Float> faceAreas(3*stride, 0.0);
  std::vector<Float> faceCenters(3*stride, 0.0);

  for (Int i = 0; i < stride; i++) {
    Int face = pf[i];
    Float *fa = &faceAreas[3*i];
    Float *fc = &faceCenters[3*i];
    Int np = INDPG[i+1]-INDPG[i];
    Int *pn = &NGON[INDPG[i]];
    for (Int j = 0; j < np; j++) pn[j] += 1;
    K_METRIC::compute_face_center_and_area(face, np, pn, x, y, z, fc, fa);
  }
  
  // Estimate cell centroid as average of face centers
  Float cc[3] = {0,0,0};
  for (Int i = 0; i < stride; i++) {
    Float *fc = &faceCenters[3*i];
    for (Int j = 0; j < 3; j++)
      cc[j] += fc[j];
  }
  for (Int i = 0; i < 3; i++)
    cc[i] /= stride;

  // Compute cell volume
  vol = 0.0;

  for (Int i = 0; i < stride; i++) {
    Float *fa = &faceAreas[3*i];
    Float *fc = &faceCenters[3*i];

    // Compute 3*face-pyramid volume contribution
    Float d[3] = {fc[0]-cc[0], fc[1]-cc[1], fc[2]-cc[2]};
    //Float pyr3Vol = K_MATH::dot(fa, fc, 3);
    Float pyr3Vol = K_MATH::dot(fa, d, 3);

    vol += pyr3Vol;
  }

  vol *= K_MATH::ONE_THIRD;
}

static
void Flag_and_get_external_faces(IMesh *M, std::vector<Int> &fflags,
  std::vector<Int> &efaces)
{
  std::vector<Int> face_count(M->nf, 0);

  // Loop through the elements and increment face_count
  for (Int i = 0; i < M->nc; i++) {
    Int stride = -1;
    Int *pf = mesh_get_cell(i, stride, M);
    for (Int j = 0; j < stride; ++j)
      face_count[pf[j]]++;
  }

  // External faces are those with a count equal to 1
  fflags.resize(M->nf);
  for (Int i = 0; i < M->nf; i++) {
    if (face_count[i] == 1) {
      fflags[i] = EXTERNAL;
      efaces.push_back(i);
    } else {
      fflags[i] = INTERNAL;
    }
  }
}

static
Int _Orient_boundary
(
  IMesh *M,
  Float *x, Float *y, Float *z,
  Int ncells,
  Int *efadj, Int *efxadj, Int nefaces,
  Int *fneis, Int *efaces, std::vector<Int> &forient,
  const std::vector<Int> &cflags, const std::vector<Int> &fflags,
  Int *cells)
{
  // Look for a cell whose volume is "definitely well computed"
  Float cvol = 0.0;
  Int seed = -1;
  Int refPG = -1;
  Int refIdx = -1;

  while (++seed < ncells) {
    if (cflags[seed] != EXTERNAL) continue;

    Int cid = (cells != NULL) ? cells[seed] : seed;

    Int stride = -1;
    Int *pf = mesh_get_cell(cid, stride, M);
    refPG = -1;
    Int local_idx = -1;
    for (Int j = 0; j < stride; j++) {
      Int face = pf[j];
      if (fflags[face] == EXTERNAL) {
        refPG = face;
        local_idx = j;
        break;
      }
    }

    if (refPG == -1) {
      fprintf(stderr, "orient_boundary(): couldn't find an external face within external cell " SF_D_ "\n", cid);
      return 1;
    }

    // Look for index of refPG in efaces (0-based)
    refIdx = -1;
    for (Int i = 0; i < nefaces; i++) {
      if (efaces[i] == refPG) {
        refIdx = i;
        break;
      }
    }

    if (refIdx == -1) {
      fprintf(stderr, "orient_boundary(): couldn't find reference face " SF_D_ " in external faces list\n", refPG);
      return 1;
    }

    // Set orientation of refPG to +1.
    // Reorient seed's faces based on orientation of refPG.
    // Compute cvol, the volume of seed.
    // If cvol > 0, orientation of all faces including refPG, is outwards
    // Otherwise, set orientation of refPG to -1.

    Compute_cell_volume(cid, M, x, y, z, cvol, local_idx);

    if (fabs(cvol) < DSMALL) continue;

    // set reference orientation of refPG and exit
    forient[refIdx] = (cvol > 0.0) ? 1 : -1;

    break;
  }

   if (seed >= ncells) {
    fprintf(stderr, "orient_boundary_ngon(): couldn't find reference polyhedron\n");
    assert(0);
    return 1;
  }

  // propagate
  K_CONNECT::reversi_connex(efadj, efxadj, nefaces, fneis, refIdx, forient);

  return 0;
}

static
void Extract_nface_of_kept_pgs(IMesh *M, const std::vector<bool> &kept_pgs,
  std::vector<Int> &NFACE, std::vector<Int> &xadj,
  std::vector<Int> &cells)
{
  Int ncells = M->nc;

  NFACE.clear();
  xadj.resize(1, 0);
  cells.clear();

  for (Int i = 0; i < ncells; i++) {
    Int stride = -1;
    Int *pf = mesh_get_cell(i, stride, M);
    bool keep = false;
    for (Int j = 0; j < stride && !keep; j++)
      keep = kept_pgs[pf[j]];
    if (keep) {
      cells.push_back(i);
      xadj.push_back(stride);
      for (Int j = 0; j < stride; j++)
        NFACE.push_back(pf[j]);
    }
  }

  for (size_t i = 0; i < xadj.size(); i++)
    xadj[i+1] += xadj[i];
}

static
void Flag_marked_external_cells(IMesh *M, const std::vector<Int> &cells,
  const std::vector<Int> &fflags, std::vector<Int> &cflags)
{
  // External cells are those with at least one external face
  cflags.resize(cells.size(), INTERNAL);
  for (size_t i = 0; i < cells.size(); i++) {
    Int cell = cells[i];
    Int stride = -1;
    Int *pf = mesh_get_cell(cell, stride, M);
    for (Int j = 0; j < stride; j++) {
      Int face = pf[j];
      if (fflags[face] == EXTERNAL) {
        cflags[i] = EXTERNAL;
        break;
      }
    }
  }
}

static
void Flag_all_external_cells(IMesh *M, const std::vector<Int> &fflags,
  std::vector<Int> &cflags)
{
  Int ncells = M->nc;

  // External cells are those with at least one external face
  cflags.resize(ncells, INTERNAL);
  for (Int i = 0; i < ncells; i++) {
    Int stride = -1;
    Int *pf = mesh_get_cell(i, stride, M);
    for (Int j = 0; j < stride; j++) {
      Int face = pf[j];
      if (fflags[face] == EXTERNAL) {
        cflags[i] = EXTERNAL;
        break;
      }
    }
  }
}

Int mesh_orient_boundary(IMesh *M)
{
  Int nfaces = M->nf;
  Int ncells = M->nc;
  Float *x = &(M->x[0]);
  Float *y = &(M->y[0]);
  Float *z = &(M->z[0]);

  // flag external cells and faces
  std::vector<Int> fflags, efaces;
  Flag_and_get_external_faces(M, fflags, efaces);

  // extract external faces connectivity
  std::vector<Int> fadj;
  std::vector<Int> xadj(1, 0);
  for (Int i = 0; i < nfaces; i++) {
    if (fflags[i] == EXTERNAL) {
      Int stride = -1;;
      Int *pn = mesh_get_face(i, stride, M);
      xadj.push_back(stride);
      for (Int j = 0; j < stride; j++)
        fadj.push_back(pn[j]);
    }
  }

  Int nefaces = (Int)efaces.size();

  for (Int i = 0; i < nefaces; i++)
    xadj[i+1] += xadj[i];

  // build skin neighbourhood
  std::vector<Int> fneighbours;
  K_CONNECT::build_face_neighbourhood(fadj, xadj, fneighbours);

  // color the faces by connex part
  std::vector<Int> colors(xadj.size()-1);
  Int nconnex = K_CONNECT::colorConnexParts(&fneighbours[0], &xadj[0],
    nefaces, &colors[0]);

  printf("orient_boundary(): connex parts: " SF_D_ "\n", nconnex);

  assert(efaces.size() == xadj.size()-1);
  std::vector<Int> forient(nefaces, 0);
  std::vector<Int> cflags;
  Int ret = 0;
  if (nconnex > 1) {
    // extract nconnex nface-ngon for separate orientation
    for (Int color = 0; color < nconnex; color++) {
      std::vector<bool> keep_pgs(nfaces, false);
      for (Int i = 0; i < nefaces; i++) {
        keep_pgs[efaces[i]] = (colors[i] == color);
      }
      // extract nface corresponding to kept faces
      std::vector<Int> NFACE, cxadj(1, 0), cells;
      Extract_nface_of_kept_pgs(M, keep_pgs, NFACE, cxadj, cells);

      std::vector<Int> cflags;
      Flag_marked_external_cells(M, cells, fflags, cflags);

      ret |= _Orient_boundary(M, x, y, z, (Int)cells.size(),
        &fadj[0], &xadj[0], nefaces, &fneighbours[0], &efaces[0], forient,
        cflags, fflags, &cells[0]);
    }
  } else {
    std::vector<Int> cflags;
    Flag_all_external_cells(M, fflags, cflags);
    ret = _Orient_boundary(M, x, y, z, ncells, &fadj[0], &xadj[0],
      nefaces, &fneighbours[0], &efaces[0], forient, cflags, fflags, NULL);
  }

  // Apply orientation
  Int nrev = 0;
  for (Int i = 0; i < nefaces; i++) {
    if (forient[i] == -1) {
      Int face = efaces[i]; // 0-based
      Int stride = -1;
      Int *pn = mesh_get_face(face, stride, M);
      std::reverse(pn+1, pn+stride);
      nrev++;
    }
  }
  printf("orient_boundary(): reversed " SF_D_ " faces\n", nrev);

  return ret;
}

static
void Build_cell_neighbourhood(IMesh *M, std::vector<Int>& neighbours,
  std::vector<Int> &xadj)
{
  Int nfaces = M->nf;
  Int ncells = M->nc;

  xadj.resize(ncells+1);
  xadj[0] = 0;

  // TODO(Imad): this is a safe resize
  Int *ptr = &xadj[0]+1;
  for (Int i = 0; i < ncells; i++) {
    Int stride = -1;
    mesh_get_cell(i, stride, M);
    *ptr++ = stride;
  }

  for (Int i = 0; i < ncells; i++) xadj[i+1] += xadj[i];

  Int sz = xadj[ncells];
  neighbours.resize(sz, -1);

  std::vector<Int> neigh(nfaces, -1);

  Int count = 0;
  while (count++ != 2) {
    for (Int i = 0; i < ncells; i++) {
      Int stride = -1;
      Int *pf = mesh_get_cell(i, stride, M);
      Int *pn = &neighbours[xadj[i]];
      for (Int j = 0; j < stride; j++) {
        Int face = pf[j];
        Int &nei = neigh[face];
        Int &Kn = pn[j];
        if (nei != -1 && nei != i)
          Kn = nei;
        neigh[face] = i;
      }
    }
  }
}

Int mesh_build_own_nei(IMesh *M)
{
  Int ncells = M->nc;
  Int *owner = &(M->owner[0]);
  Int *neigh = &(M->neigh[0]);

  std::vector<Int> neighbours, xadj;
  Build_cell_neighbourhood(M, neighbours, xadj);

  std::vector<Int> exPH(ncells, -1);
  for (Int i = 0; i < ncells; i++) {
    Int stride = -1;
    Int *pf = mesh_get_cell(i, stride, M);
    Int *pn = &neighbours[xadj[i]];
    for (Int j = 0; j < stride; j++) {
      if (pn[j] == -1) {
        owner[pf[j]] = i;
        exPH[i] = pf[j]+1;
        break;
      }
    }
  }

  // look for first external cell
  std::vector<Int> processed(ncells, 0);

  Int nconnex = 0;

  Int seed = 0;

  std::stack<Int> cpool;

   while (1) {
    while ((seed < ncells) && (processed[seed] || exPH[seed] == -1))
      ++seed;

    if (seed >= ncells)
      break;

    nconnex++;

    cpool.push(seed);

    while (!cpool.empty()) {
      Int cell = cpool.top();
      cpool.pop();

      if (processed[cell])
        continue;

      processed[cell] = 1;

      // build faces neighbourhood based on shared edges' nodes order
      Int stride = -1;
      Int *pf = mesh_get_cell(cell, stride, M);
      std::vector<Int> oids;
      std::vector<Int> orient(stride, 1);
      std::vector<Int> pgs;
      std::vector<Int> xpgs(1, 0);
      for (Int i = 0; i < stride; i++) {
        Int face = pf[i];
        Int np = -1;
        Int *pn = mesh_get_face(face, np, M);
        for (Int j = 0; j < np; j++)
          pgs.push_back(pn[j]);
        xpgs.push_back(np);
        oids.push_back(face+1);
      }

      for (Int i = 0; i < stride; i++)
        xpgs[i+1] += xpgs[i];

      std::vector<Int> PGneighbours(pgs.size());
      K_CONNECT::build_face_neighbourhood(pgs, xpgs, PGneighbours);

      Int revers = 0;

      // reference face is the external face
      Int PGref = exPH[cell];
      // face can be negative
      if (PGref < 0) {
        revers = 1;
        PGref = -PGref;
      }

      // find reference face index in oids
      Int iref = -1;
      for (size_t i = 0; i < oids.size(); i++) {
        if (PGref == oids[i]) {
          iref = i;
          break;
        }
      }
      assert(iref != -1);

      // set orientation of face if prescribed
      if (revers)
        orient[iref] = -1;

      // all connected faces must follow the orientation of the reference face
      K_CONNECT::reversi_connex(&pgs[0], &xpgs[0], stride, &PGneighbours[0],
        iref, orient);

       // set the owner and neighbour of the faces
      Int *pn = &neighbours[xadj[cell]];
      for (Int i = 0; i < stride; i++) {
        Int face = pf[i];
        Int nei = pn[i];
        assert(nei < ncells && nei >= -1);

        owner[face] = cell;
        neigh[face] = nei;

        if (nei == -1)
          continue;

        // set the reference face for neighbour
        exPH[nei] = -(face+1);

        if (orient[i] == -1) {
          std::swap(owner[face], neigh[face]);
          exPH[nei] = face+1;
        }

        if (!processed[nei])
          cpool.push(nei);
      }
    }
  }

  //printf("build_parent_elements(): connex parts: " SF_D_ "\n", nconnex);

  return 0;
}

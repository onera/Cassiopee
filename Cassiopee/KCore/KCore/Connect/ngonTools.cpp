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
#include "connect.h"
#include "String/kstring.h"
#include <map>
#include <algorithm>
#include <stack>
#include "Metric/metric.h"
#include <unordered_map>

#define EXTERNAL 1
#define INTERNAL 0
#define DSMALL 1e-18

struct edge {
  E_Int p0_;
  E_Int p1_;

  edge()
  {}

  edge(E_Int p0, E_Int p1)
  :
  p0_(std::min(p0, p1)),
  p1_(std::max(p0, p1))
  {}

  void set(E_Int p0, E_Int p1)
  {
    p0_ = std::min(p0, p1);
    p1_ = std::max(p0, p1);
  }

  bool operator<(const edge &a) const
  {
    return (p0_ < a.p0_) || (p0_ == a.p0_ && p1_ < a.p1_);
  }
};

// Return 0 if closed, 1 if open
static
E_Int check_open_cell(E_Int cell, K_FLD::FldArrayI &cn)
{
  // Closed cell: edge count = 2 for all edges
  E_Int stride = -1;
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();
  E_Int *pf = cn.getElt(cell, stride, nface, indPH);
  edge E;
  std::map<edge, E_Int> edgeCount;
  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i]-1;
    E_Int np = -1;
    E_Int *pn = cn.getFace(face, np, ngon, indPG);
    for (E_Int j = 0; j < np; j++) {
      E_Int p0 = pn[j];
      E_Int p1 = pn[(j+1)%np];
      E.set(p0, p1);
      edgeCount[E]++;
    }
  }

  for (auto& ec : edgeCount) {
    if (ec.second != 2) {
      fprintf(stderr, "Warning: Cell " SF_D_ " is not closed.\n", cell); 
      return 1;
    }
  }

  return 0;
}

// Returns 0 if all cells are closed, returns 1 otherwise
E_Int K_CONNECT::check_open_cells(K_FLD::FldArrayI& cn, E_Int *is_cell_open)
{
  E_Int ncells = cn.getNElts();
  
  if (is_cell_open) {
    // Up to the caller to allocate is_cell_open
    for (E_Int i = 0; i < ncells; i++)
      is_cell_open[i] = check_open_cell(i, cn);

    for (E_Int i = 0; i < ncells; i++) {
      if (is_cell_open[i])
        return 1;
    }
  } else {
    for (E_Int i = 0; i < ncells; i++) {
      if (check_open_cell(i, cn))
        return 1;
    }
  }
  
  return 0;
}

E_Int K_CONNECT::check_overlapping_cells(K_FLD::FldArrayI &cn)
{
  E_Int *nface = cn.getNFace();
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();

  E_Int ret = 0;

  std::unordered_map<E_Int, std::vector<E_Int>> f2e;

  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = cn.getElt(i, stride, nface, indPH);
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j];
      f2e[face].push_back(i);
    }
  }

  for (const auto &FE : f2e) {
    if (FE.second.size() != 1 && FE.second.size() != 2) {
      fprintf(stderr, "Warning: Face " SF_D_ " belongs to more than two cells\n",
        FE.first);
      ret = 1;
    }
  }

  return ret;
}

void K_CONNECT::build_face_neighbourhood
(
  std::vector<E_Int> &pgs,
  std::vector<E_Int> &xpgs,
  std::vector<E_Int> &neighbour
)
{
  neighbour.resize(pgs.size(), -1);

  std::map<edge, std::pair<std::pair<E_Int, E_Int>, std::pair<E_Int, E_Int>>> EM;

  size_t nf = xpgs.size() - 1;

  for (size_t i = 0; i < nf; i++) {
    E_Int start = xpgs[i];
    E_Int end = xpgs[i+1];
    E_Int stride = end - start;

    E_Int *pn = &pgs[start];
    for (E_Int j = 0; j < stride; j++) {
      E_Int n0 = pn[j];
      E_Int n1 = pn[(j+1)%stride];

      edge E(n0, n1);

      auto search = EM.find(E);

      if (search == EM.end()) {
        // first time encoutering this edge
        EM[E].first = std::make_pair(i, j);
        EM[E].second = std::make_pair(-1, -1);
      } else {
        if (search->second.second.first == -1)
          search->second.second = std::make_pair(i, j);
        else
          search->second.second = std::make_pair(E_IDX_NONE, E_IDX_NONE);
      }
    }
  }

  for (auto &elem : EM) {
    E_Int pg0 = elem.second.first.first;
    E_Int n0 = elem.second.first.second;
    E_Int pg1 = elem.second.second.first;
    E_Int n1 = elem.second.second.second;

    if (pg1 == -1 || pg1 == E_IDX_NONE)
      continue;

    E_Int s0 = xpgs[pg0];
    E_Int s1 = xpgs[pg1];

    neighbour[s0 + n0] = pg1;
    neighbour[s1 + n1] = pg0;
  }

  // handle non-manifoldness
  std::map<edge, E_Int> edge_to_count;
  edge E;
  for (size_t i = 0; i < nf; i++) {
    E_Int start = xpgs[i];
    E_Int end = xpgs[i+1];
    E_Int stride = end-start;
    E_Int *pn = &pgs[start];
    for (E_Int j = 0; j < stride; j++) {
      E_Int ni = pn[j];
      E_Int nj = pn[(j+1)%stride];
      E.set(ni, nj);
      auto it = edge_to_count.find(E);
      if (it == edge_to_count.end())
        edge_to_count.insert(std::make_pair(E, 1));
      else
        it->second++;
    }
  }

  for (size_t i = 0; i < nf; i++) {
    E_Int start = xpgs[i];
    E_Int end = xpgs[i+1];
    E_Int stride = end-start;
    E_Int *pn = &pgs[start];
    E_Int *pk = &neighbour[start];
    for (E_Int j = 0; j < stride; j++) {
      E_Int ni = pn[j];
      E_Int nj = pn[(j+1)%stride];
      E.set(ni, nj);
      if (edge_to_count[E] != 2)
        pk[j] = -1;
    }
  }
}

static
E_Int get_orientation(E_Int *pn, E_Int stride, E_Int ni, E_Int nj,
  E_Int *same_orient)
{
  *same_orient = 0;
  for (E_Int i = 0; i < stride; i++) {
    if (pn[i] == ni && pn[(i+1)%stride] == nj) {
      *same_orient = 1;
      return 0;
    }
    if (pn[i] == nj && pn[(i+1)%stride] == ni) {
      *same_orient = 0;
      return 0;
    }
  }
  return -1;
}

static
void get_boundary(E_Int *pn0, E_Int s0, E_Int *pn1, E_Int s1, E_Int *m,
  E_Int *n)
{
  for (E_Int i = 0; i < s0; i++) {
    E_Int n00 = pn0[i];
    E_Int n01 = pn0[(i+1)%s0];
    for (E_Int j = 0; j < s1; j++) {
      E_Int n10 = pn1[j];
      E_Int n11 = pn1[(j+1)%s1];
      if ((n00 == n10 || n00 == n11) && (n01 == n10 || n01 == n11)) {
        *m = i;
        *n = j;
        return;
      }
    }
  }
}

// Deduce orientation of connected faces based on orientation of seed face
void K_CONNECT::reversi_connex(E_Int *pgs, E_Int *xpgs, E_Int npgs,
  E_Int *neighbours, E_Int kseed, std::vector<E_Int> &orient)
{
  assert(kseed < npgs);
  std::vector<E_Int> cpool;
  cpool.push_back(kseed);
  
  std::vector<E_Int> processed(npgs, 0);

  while (!cpool.empty()) {
    E_Int K = cpool.back();
    cpool.pop_back();

    processed[K] = 1;

    E_Int *pf = &pgs[xpgs[K]];
    E_Int stride = xpgs[K+1]-xpgs[K];

    for (E_Int i = xpgs[K]; i < xpgs[K+1]; i++) {
      E_Int nei = neighbours[i];
      if (nei == -1 || processed[nei])
        continue;

      // get the shared edge between face K and face nei
      E_Int k, l;
      k = l = -1;
      E_Int *pnn = &pgs[xpgs[nei]];
      E_Int sn = xpgs[nei+1] - xpgs[nei];
      get_boundary(pf, stride, pnn, sn, &k, &l);

      E_Int ni = pf[k];
      E_Int nj = pf[(k+1)%stride];

      E_Int reverse = 2;
      get_orientation(pnn, sn, ni, nj, &reverse);

      if (orient[K] == -1) reverse = !reverse;
      if (reverse) orient[nei] = -1;

      cpool.push_back(nei);
    }
  }
}

static
E_Int _orient_boundary
(
  K_FLD::FldArrayI &cn,
  E_Float *x, E_Float *y, E_Float *z,
  E_Int ncells,
  E_Int *efadj, E_Int *efxadj, E_Int nefaces,
  E_Int *fneis, E_Int *efaces, std::vector<E_Int> &forient,
  const std::vector<E_Int> &cflags, const std::vector<E_Int> &fflags,
  E_Int *cells)
{
  // Look for a cell whose volume is "definitely well computed"
  E_Float cvol = 0.0;
  E_Int seed = -1;
  E_Int refPG = -1;
  E_Int refIdx = -1;

  while (++seed < ncells) { 
    if (cflags[seed] != EXTERNAL) continue;

    E_Int cid = (cells != NULL) ? cells[seed] : seed;

    if (check_open_cell(cid, cn)) {
      fprintf(stderr, "_orient_boundary(): non-closed cell found. Aborting.\n");
      return 1;
    }
    
    E_Int stride = -1;
    E_Int *pf = cn.getElt(cid, stride, cn.getNFace(), cn.getIndPH());
    refPG = -1;
    E_Int local_idx = -1;
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j]-1;
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
    for (E_Int i = 0; i < nefaces; i++) {
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

    K_METRIC::compute_cell_volume(cid, cn, x, y, z, cvol, local_idx);

    if (fabs(cvol) < DSMALL) continue;
    
    // set reference orientation of refPG and exit
    forient[refIdx] = (cvol > 0.0) ? 1 : -1;
    
    break;
  }
  
  if (seed >= ncells) {
    assert(0);
    fprintf(stderr, "orient_boundary_ngon(): couldn't find reference polyhedron\n");
    return 1;
  }

  // propagate
  K_CONNECT::reversi_connex(efadj, efxadj, nefaces, fneis, refIdx, forient);

  return 0;
}

static
void extract_nface_of_kept_pgs(K_FLD::FldArrayI &cn, const std::vector<bool> &kept_pgs,
  std::vector<E_Int> &NFACE, std::vector<E_Int> &xadj, std::vector<E_Int> &cells)
{
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int ncells = cn.getNElts();

  NFACE.clear();
  xadj.resize(1, 0);
  cells.clear();

  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = cn.getElt(i, stride, nface, indPH);
    bool keep = false;
    for (E_Int j = 0; j < stride && !keep; j++)
      keep = kept_pgs[pf[j]-1];
    if (keep) {
      cells.push_back(i);
      xadj.push_back(stride);
      for (E_Int j = 0; j < stride; j++)
        NFACE.push_back(pf[j]);
    }
  }

  for (size_t i = 0; i < xadj.size(); i++)
    xadj[i+1] += xadj[i];
}

static
void flag_and_get_external_faces(K_FLD::FldArrayI &cn, std::vector<E_Int> &fflags,
  std::vector<E_Int> &efaces)
{
  E_Int nfaces = cn.getNFaces();
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();
  E_Int *nface = cn.getNFace();

  std::vector<E_Int> face_count(nfaces, 0);
  
  // Loop through the elements and increment face_count
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = cn.getElt(i, stride, nface, indPH);
    for (E_Int j = 0; j < stride; ++j)
      face_count[pf[j]-1]++;
  }
  
  // External faces are those with a count equal to 1
  fflags.resize(nfaces);
  for (E_Int i = 0; i < nfaces; i++) {
    if (face_count[i] == 1) {
      fflags[i] = EXTERNAL;
      efaces.push_back(i);
    } else {
      fflags[i] = INTERNAL;
    }
  }
}

static
void flag_marked_external_cells(K_FLD::FldArrayI &cn, const std::vector<E_Int> &cells,
  const std::vector<E_Int> &fflags, std::vector<E_Int> &cflags)
{
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();

  // External cells are those with at least one external face
  cflags.resize(cells.size(), INTERNAL);
  for (size_t i = 0; i < cells.size(); i++) {
    E_Int cell = cells[i];
    E_Int stride = -1;
    E_Int *pf = cn.getElt(cell, stride, nface, indPH);
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j]-1;
      if (fflags[face] == EXTERNAL) {
        cflags[i] = EXTERNAL;
        break;
      }
    }
  }
}

static
void flag_all_external_cells(K_FLD::FldArrayI &cn, const std::vector<E_Int> &fflags,
  std::vector<E_Int> &cflags)
{
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH(); 
  E_Int ncells = cn.getNElts();

  // External cells are those with at least one external face
  cflags.resize(ncells, INTERNAL);
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = cn.getElt(i, stride, nface, indPH);
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j]-1;
      if (fflags[face] == EXTERNAL) {
        cflags[i] = EXTERNAL;
        break;
      }
    }
  }
}

E_Int K_CONNECT::orient_boundary_ngon(E_Float *x, E_Float *y, E_Float *z,
  K_FLD::FldArrayI &cn)
{
  E_Int nfaces = cn.getNFaces();
  E_Int *indPG = cn.getIndPG();
  E_Int *ngon = cn.getNGon();
  E_Int ncells = cn.getNElts();
  //E_Int *indPH = cn.getIndPH();
  //E_Int *nface = cn.getNFace();
  
  // flag external cells and faces
  std::vector<E_Int> fflags, efaces;
  flag_and_get_external_faces(cn, fflags, efaces);

  // extract external faces connectivity
  std::vector<E_Int> fadj;
  std::vector<E_Int> xadj(1, 0);
  for (E_Int i = 0; i < nfaces; i++) {
    if (fflags[i] == EXTERNAL) {
      E_Int stride = -1;
      E_Int *pn = cn.getFace(i, stride, ngon, indPG);
      xadj.push_back(stride);
      for (E_Int j = 0; j < stride; j++)
        fadj.push_back(pn[j]);
    }
  }

  E_Int nefaces = (E_Int)efaces.size();

  for (E_Int i = 0; i < nefaces; i++)
    xadj[i+1] += xadj[i];

  // build skin neighbourhood
  std::vector<E_Int> fneighbours;
  build_face_neighbourhood(fadj, xadj, fneighbours);

  // color the faces by connex part
  std::vector<E_Int> colors(xadj.size()-1);
  E_Int nconnex = K_CONNECT::colorConnexParts(&fneighbours[0], &xadj[0],
    nefaces, &colors[0]);
  
  //printf("orient_boundary(): connex parts: " SF_D_ "\n", nconnex);

  assert(efaces.size() == xadj.size()-1);

  std::vector<E_Int> forient(nefaces, 0);
  std::vector<E_Int> cflags;
  E_Int ret = 0;
  if (nconnex > 1) {
    // extract nconnex nface-ngon for separate orientation
    for (E_Int color = 0; color < nconnex; color++) {
      std::vector<bool> keep_pgs(nfaces, false);
      for (E_Int i = 0; i < nefaces; i++) {
        keep_pgs[efaces[i]] = (colors[i] == color);
      }
      // extract nface corresponding to kept faces
      std::vector<E_Int> NFACE, cxadj(1, 0), cells;
      extract_nface_of_kept_pgs(cn, keep_pgs, NFACE, cxadj, cells);

      std::vector<E_Int> cflags;
      flag_marked_external_cells(cn, cells, fflags, cflags);

      ret |= _orient_boundary(cn, x, y, z, (E_Int)cells.size(), 
        &fadj[0], &xadj[0], nefaces, &fneighbours[0], &efaces[0], forient, cflags,
        fflags, &cells[0]);
    }
  } else {
    std::vector<E_Int> cflags;
    flag_all_external_cells(cn, fflags, cflags);
    ret = _orient_boundary(cn, x, y, z, ncells, &fadj[0], &xadj[0],
      nefaces, &fneighbours[0], &efaces[0], forient, cflags, fflags, NULL);
  }
  
  // Apply orientation
  //E_Int nrev = 0;
  for (E_Int i = 0; i < nefaces; i++) {
    if (forient[i] == -1) {
      E_Int face = efaces[i]; // 0-based
      E_Int stride = -1;
      E_Int *pn = cn.getFace(face, stride, ngon, indPG);
      std::reverse(pn+1, pn+stride);
      //nrev++;
    }
  }
  //printf("Orient_boundary(): flipped " SF_D_ " faces\n", nrev);

  return ret;
}

static
void build_cell_neighbourhood(K_FLD::FldArrayI &cn, std::vector<E_Int>& neighbours,
  std::vector<E_Int> &xadj)
{

  E_Int nfaces = cn.getNFaces();
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();
  E_Int *nface = cn.getNFace();

  xadj.resize(ncells+1);
  xadj[0] = 0;

  // TODO(Imad): this is a safe resize
  E_Int *ptr = &xadj[0]+1;
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    cn.getElt(i, stride, nface, indPH);
    *ptr++ = stride;
  }

  for (E_Int i = 0; i < ncells; i++) xadj[i+1] += xadj[i];

  E_Int sz = xadj[ncells];
  neighbours.resize(sz, -1);

  std::vector<E_Int> neigh(nfaces, -1);
  
  E_Int count = 0;
  while (count++ != 2) {
    for (E_Int i = 0; i < ncells; i++) {
      E_Int stride = -1;
      E_Int *pf = cn.getElt(i, stride, nface, indPH);
      E_Int *pn = &neighbours[xadj[i]];
      for (E_Int j = 0; j < stride; j++) {
        E_Int face = pf[j]-1;
        E_Int &nei = neigh[face];
        E_Int &Kn = pn[j];
        if (nei != -1 && nei != i)
          Kn = nei;
        neigh[face] = i;
      }
    }
  }
}

// Assumes that the external faces follow the same orientation (all outwards/inwards)
// Arrays owner and neigh should be allocated by caller
// Returns zero-based parent elements
E_Int K_CONNECT::build_parent_elements_ngon(K_FLD::FldArrayI &cn, E_Int *owner,
  E_Int *neigh)
{
  E_Int ncells = cn.getNElts();
  E_Int *nface = cn.getNFace();
  E_Int *ngon = cn.getNGon();
  E_Int *indPH = cn.getIndPH();
  E_Int *indPG = cn.getIndPG();
  
  std::vector<E_Int> neighbours, xadj;
  build_cell_neighbourhood(cn, neighbours, xadj);

  std::vector<E_Int> exPH(ncells, -1);
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = cn.getElt(i, stride, nface, indPH);
    E_Int *pn = &neighbours[xadj[i]];
    for (E_Int j = 0; j < stride; j++) {
      if (pn[j] == -1) {
        owner[pf[j]-1] = i;
        exPH[i] = pf[j];
        break;
      }
    }
  }

  // look for first external cell
  std::vector<E_Int> processed(ncells, 0);

  E_Int nconnex = 0;
  
  E_Int seed = 0;

  std::stack<E_Int> cpool;

  while (1) {
    while ((seed < ncells) && (processed[seed] || exPH[seed] == -1))
      ++seed;

    if (seed >= ncells)
      break;
    
    nconnex++;

    cpool.push(seed);

    while (!cpool.empty()) {
      E_Int cell = cpool.top();
      cpool.pop();

      if (processed[cell])
        continue;

      processed[cell] = 1;

      // build faces neighbourhood based on shared edges' nodes order
      E_Int stride = -1;
      E_Int *pf = cn.getElt(cell, stride, nface, indPH);
      std::vector<E_Int> oids;
      std::vector<E_Int> orient(stride, 1);
      std::vector<E_Int> pgs;
      std::vector<E_Int> xpgs(1, 0);
      for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i];
        E_Int np = -1;
        E_Int *pn = cn.getFace(face-1, np, ngon, indPG);
        for (E_Int j = 0; j < np; j++)
          pgs.push_back(pn[j]);
        xpgs.push_back(np);
        oids.push_back(face);
      }

      for (E_Int i = 0; i < stride; i++)
        xpgs[i+1] += xpgs[i];

      std::vector<E_Int> PGneighbours(pgs.size());
      build_face_neighbourhood(pgs, xpgs, PGneighbours);

      E_Int revers = 0;

      // reference face is the external face
      E_Int PGref = exPH[cell];
      // face can be negative
      if (PGref < 0) {
        revers = 1;
        PGref = -PGref;
      }

      // find reference face index in oids
      E_Int iref = -1;
      for (size_t i = 0; i < oids.size(); i++) {
        if (PGref == oids[i]) {
          iref = i;
          break;
        }
      }

      if (iref == -1) {
        break;
        return 1;
      }

      assert(iref != -1);

      // set orientation of face if prescribed
      if (revers)
        orient[iref] = -1;

      // all connected faces must follow the orientation of the reference face
      reversi_connex(&pgs[0], &xpgs[0], stride, &PGneighbours[0], iref, orient);

      // set the owner and neighbour of the faces
      E_Int *pn = &neighbours[xadj[cell]];
      for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i];
        E_Int nei = pn[i];
        assert(nei < ncells && nei >= -1);

        owner[face-1] = cell;
        neigh[face-1] = nei;

        if (nei == -1)
          continue;

        // set the reference face for neighbour
        exPH[nei] = -face;

        if (orient[i] == -1) {
          std::swap(owner[face-1], neigh[face-1]);
          exPH[nei] = face;
        }

        if (!processed[nei])
          cpool.push(nei);
      }
    }
  }

  //printf("build_parent_elements(): connex parts: " SF_D_ "\n", nconnex);
  
  return 0;
}

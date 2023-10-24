#include "metric.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <stack>

#define INTERNAL 0
#define EXTERNAL 1
#define ONE_THIRD 0.3333333333333333

E_Float dot(E_Float *a, E_Float *b, E_Int n)
{
  E_Float res = 0;
  for (E_Int i = 0; i < n; i++)
    res += a[i]*b[i];
  return res;
}

void cross(E_Float a[3], E_Float b[3], E_Float c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

E_Float norm(E_Float *a, E_Int n)
{
  return sqrt(dot(a, a, n));
}

static
void compute_face_center_and_area(E_Float *x, E_Float *y, E_Float *z, E_Int *pn,
  E_Int stride, E_Float tol, E_Float fc[3], E_Float fa[3])
{
  fa[0] = fa[1] = fa[2] = 0.0;
  fc[0] = fc[1] = fc[2] = 0.0;

  // Approximate face center
    E_Float fcenter[3] = {0,0,0};
    for (E_Int j = 0; j < stride; j++) {
      E_Int point = pn[j]-1;
      fcenter[0] += x[point];
      fcenter[1] += y[point];
      fcenter[2] += z[point];
    }
    fcenter[0] /= stride;
    fcenter[1] /= stride;
    fcenter[2] /= stride;

    // Sum of triangle area vectors
    E_Float sumN[3] = {0,0,0};
    // Sum of triangle areas
    E_Float sumA = 0;
    // Area-weighted sum of triangle centers
    E_Float sumAc[3] = {0,0,0};

    // Compute area vector and center of stride-2 triangles
    // formed by p0p1p2, p0p2p3, ... p0p(stride-2)p(stride-1)
    E_Int p0 = pn[0]-1;
    for (E_Int j = 1; j < stride-1; j++) {
      E_Int p1 = pn[j]-1;
      E_Int p2 = pn[j+1]-1;

      // Triangle center
      E_Float tc[3];
      tc[0] = x[p0] + x[p1] + x[p2];
      tc[1] = y[p0] + y[p1] + y[p2];
      tc[2] = z[p0] + z[p1] + z[p2];

      // Area vector
      E_Float n[3];
      E_Float v10[3] = {x[p1]-x[p0], y[p1]-y[p0], z[p1]-z[p0]};
      E_Float v20[3] = {x[p2]-x[p0], y[p2]-y[p0], z[p2]-z[p0]};
      cross(v10, v20, n);

      // Area
      E_Float a = norm(n, 3);

      for (E_Int k = 0; k < 3; k++) {
        sumN[k] += n[k];
        sumAc[k] += a*tc[k];
      }
      sumA += a;
    }

    // Deal with zero-area faces
    if (sumA < tol) {
      for (E_Int j = 0; j < 3; j++) {
        fc[j] = fcenter[j];
        fa[j] = 0.0;
      }
    } else {
      for (E_Int j = 0; j < 3; j++) {
        fc[j] = sumAc[j]/(3.0*sumA);
        fa[j] = 0.5*sumN[j];
      }
    }
}


// Assumes consistent face orientation
static
E_Float compute_cell_volume(E_Float *x, E_Float *y, E_Float *z, K_FLD::FldArrayI &cn,
  E_Int cell, E_Float tol)
{
  E_Int *nface = cn.getNFace();
  E_Int *ngon = cn.getNGon();
  E_Int *indPH = cn.getIndPH();
  E_Int *indPG = cn.getIndPG();

  E_Int stride = -1;
  E_Int *pf = cn.getElt(cell, stride, nface, indPH);
  E_Float vol = 0.;

  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i]-1;
    E_Int np = -1;
    E_Int *pn = cn.getFace(face, np, ngon, indPG);

    // face center and area
    E_Float fc[3];
    E_Float fa[3];
    compute_face_center_and_area(x, y, z, pn, np, tol, fc, fa);

    E_Float pyr3Vol = dot(fa, fc, 3);
    
    vol += pyr3Vol;
  }

  vol /= 3.0;

  return vol;
}


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

static
void build_face_neighbourhood
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

static
void reversi_connex(E_Int *pgs, E_Int *xpgs, E_Int npgs, E_Int *neighbours, E_Int kseed,
  std::vector<E_Int> &orient)
{
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
      if (processed[nei])
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
  E_Int *nface, E_Int *indPH, E_Int ncells,
  E_Int *efadj, E_Int *efxadj, E_Int nefaces,
  E_Int *fneis, E_Int *efaces, std::vector<E_Int> &forient,
  const std::vector<E_Int> &cflags, const std::vector<E_Int> &fflags,
  E_Int *cells, E_Float tol)
{
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();

  // Look for a cell whose volume is "definitely well computed"
  E_Float cvol = 0.0;
  E_Int seed = -1;
  E_Int refPG = -1;
  E_Int refIdx = -1;
  while (++seed < ncells) { 
    if (cflags[seed] != EXTERNAL) continue;
    
    E_Int stride = -1;
    E_Int *pf = cn.getElt(seed, stride, nface, indPH);
    refPG = -1;
    E_Int idx = -1;
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j]-1;
      if (fflags[face] == EXTERNAL) {
        refPG = face;
        idx = j;
        break;
      }
    }

    if (refPG == -1) {
      fprintf(stderr, "orient_boundary_ngon(): couldn't find an external polygon within external polyhedron %d\n", seed);
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
      fprintf(stderr, "orient_boundary_ngon(): couldn't find reference polygon %d in external faces list\n", refPG);
      return 1;
    }

    // Set orientation of refPG to +1. 
    // Reorient seed's faces based on orientation of refPG.
    // Compute cvol, the volume of seed.
    // If cvol > 0, orientation of all faces including refPG, is outwards
    // Otherwise, set orientation of refPG to -1.

    std::vector<E_Int> NGON;
    std::vector<E_Int> INDPG(1, 0);
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

    std::vector<E_Int> orient(stride);
    orient[idx] = 1;
    std::vector<E_Int> neis(NGON.size());
    build_face_neighbourhood(NGON, INDPG, neis);
    reversi_connex(&NGON[0], &INDPG[0], stride, &neis[0], idx, orient);
    
    // apply orientation directly in ngon
    for (E_Int i = 0; i < stride; i++) {
      if (orient[i] == -1) {
        E_Int face = pf[i]-1;
        E_Int np = -1;
        E_Int *pn = cn.getFace(face, np, ngon, indPG);
        std::reverse(pn+1, pn+np);
      }
    }

    // compute cvol
    if (cells)
      cvol = compute_cell_volume(x, y, z, cn, cells[seed], tol);
    else
      cvol = compute_cell_volume(x, y, z, cn, seed, tol);

    if (cvol == 0.0) continue;
    
    // set reference orientation of refPG and exit
    forient[refIdx] = (cvol > 0.0) ? 1 : -1;
    break;
  }
  
  if (seed >= ncells) {
    fprintf(stderr, "orient_boundary_ngon(): couldn't find reference polyhedron\n");
    return 1;
  }

  // propagate
  reversi_connex(efadj, efxadj, nefaces, fneis, refIdx, forient);

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

  // External cells are those with at least on external face
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

  // External cells are those with at least on external face
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

void K_METRIC::orient_boundary_ngon(E_Float *x, E_Float *y, E_Float *z,
  K_FLD::FldArrayI &cn, E_Float tol)
{
  E_Int nfaces = cn.getNFaces();
  E_Int *indPG = cn.getIndPG();
  E_Int *ngon = cn.getNGon();
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();
  E_Int *nface = cn.getNFace();
  
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
  E_Int nconnex = K_METRIC::colorConnexParts(&fneighbours[0], &xadj[0],
    nefaces, &colors[0]);
  
  //printf("orient_boundary(): connex parts: %d\n", nconnex);

  assert(efaces.size() == xadj.size()-1);

  std::vector<E_Int> forient(nefaces, 0);
  std::vector<E_Int> cflags;
  if (nconnex > 1) {
    // extract nconnex nface-ngon for separate orientation
    for (E_Int color = 0; color < nconnex; color++) {
      std::vector<bool> keep_pgs(nfaces, false);
      for (E_Int i = 0; i < nefaces; i++)
        keep_pgs[efaces[i]] = (colors[i] == color);
      // extract nface corresponding to kept faces
      std::vector<E_Int> NFACE, cxadj(1, 0), cells;
      extract_nface_of_kept_pgs(cn, keep_pgs, NFACE, cxadj, cells);

      std::vector<E_Int> cflags;
      flag_marked_external_cells(cn, cells, fflags, cflags);

      _orient_boundary(cn, x, y, z, &NFACE[0], &cxadj[0], (E_Int)cells.size(), 
        &fadj[0], &xadj[0], nefaces, &fneighbours[0], &efaces[0], forient, cflags,
        fflags, &cells[0], tol);
    }
  } else {
    std::vector<E_Int> cflags;
    flag_all_external_cells(cn, fflags, cflags);

    _orient_boundary(cn, x, y, z, nface, indPH, ncells, &fadj[0], &xadj[0], nefaces,
      &fneighbours[0], &efaces[0], forient, cflags, fflags, NULL, tol);
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
  //printf("orient_boundary(): reversed %d faces\n", nrev);
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

// Assumes external faces have been properly oriented outwards
void K_METRIC::build_parent_elements_ngon(K_FLD::FldArrayI &cn, std::vector<E_Int> &owner,
  std::vector<E_Int> &neigh)
{
  E_Int nfaces = cn.getNFaces();
  E_Int ncells = cn.getNElts();
  E_Int *nface = cn.getNFace();
  E_Int *ngon = cn.getNGon();
  E_Int *indPH = cn.getIndPH();
  E_Int *indPG = cn.getIndPG();
  
  std::vector<E_Int> neighbours, xadj;
  build_cell_neighbourhood(cn, neighbours, xadj);

  owner.resize(nfaces, -1);
  neigh.resize(nfaces, -1);

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

  while (1) {
    while ((seed < ncells) && (processed[seed] || exPH[seed] == -1))
      seed++;

    if (seed >= ncells)
      break;
    
    nconnex++;
    
    std::vector<E_Int> cpool;

    cpool.push_back(seed);

    while (!cpool.empty()) {
      E_Int cell = cpool.back();
      cpool.pop_back();

      if (processed[cell])
        continue;

      processed[cell] = 1;

      // build faces neighbourhood based on shared edges' nodes order
      E_Int stride = -1;
      E_Int *pf = cn.getElt(cell, stride, nface, indPH);
      std::vector<E_Int> oids;
      std::vector<E_Int> orient(stride, 1);
      std::vector<E_Int> pgs;
      std::vector<E_Int> xpgs;
      xpgs.push_back(0);
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

      // set orientation of face if prescribed
      if (revers)
        orient[iref] = -1;

      // all connected faces must follow the orientation of the reference face
      reversi_connex(&pgs[0], &xpgs[0], stride, &PGneighbours[0], iref, orient);

      // set the owner and neighbour of the faces
      E_Int *pn = &neighbours[xadj[cell]];
      for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i]-1;
        E_Int nei = pn[i];

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
          cpool.push_back(nei);
      }
    }
  }

  //printf("build_parent_elements(): connex parts: %d\n", nconnex);
}


// Assumes unsigned faces.
E_Int K_METRIC::compute_volumes_ngon(E_Float *x, E_Float *y, E_Float *z,
  K_FLD::FldArrayI &cn, E_Float *vols, E_Float tol)
{

  E_Int nfaces = cn.getNFaces();
  E_Int *indPG = cn.getIndPG();
  E_Int *ngon = cn.getNGon();
  E_Int ncells = cn.getNElts();
  
  orient_boundary_ngon(x, y, z, cn, tol);
  
  std::vector<E_Int> owner, neigh;
  build_parent_elements_ngon(cn, owner, neigh);

  // make face center and area vectors
  std::vector<E_Float> faceCenters(3*nfaces, 0);
  std::vector<E_Float> faceAreas(3*nfaces, 0);

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int stride = -1;
    E_Int *pn = cn.getFace(i, stride, ngon, indPG);
    E_Float *fc = &faceCenters[3*i];
    E_Float *fa = &faceAreas[3*i];

    compute_face_center_and_area(x, y, z, pn, stride, tol, fc, fa);
  }

  // init
  for (E_Int i = 0; i < ncells; i++)
    vols[i] = 0.0;

  for (E_Int i = 0; i < nfaces; i++) {
    E_Float *fc = &faceCenters[3*i];
    E_Float *fa = &faceAreas[3*i];

    E_Int own = owner[i];

    // Calculate 3*face-pyramid volume
    E_Float pyr3Vol = dot(fa, fc, 3);
    
    // Accumulate face-pyramid volume
    vols[own] += pyr3Vol;

    // Same for second parent element
    E_Int nei = neigh[i];
    if (nei == -1) continue;
    
    // Accumulate face-pyramid volume
    vols[nei] -= pyr3Vol;
  }

  E_Float vmin, vmax, vsum;
  E_Int cmin = 0, cmax = 0;
  
  vols[0] *= ONE_THIRD;
  vmin = vmax = vsum = vols[0];
  
  for (E_Int i = 1; i < ncells; i++) {
    vols[i] *= ONE_THIRD;
    if (vols[i] > vmax) {
      vmax = vols[i];
      cmax = i;
    }
    if (vols[i] < vmin) {
      vmin = vols[i];
      cmin = i;
    }
    vsum += vols[i];
  }

  printf("Min volume: %.6e at cell %d\n", vmin, cmin);
  printf("Max volume: %.6e at cell %d\n", vmax, cmax);
  printf("Total volume: %.6e\n", vsum);
  
  return 0;
}

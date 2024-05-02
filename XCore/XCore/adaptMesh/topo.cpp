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
#include <stdio.h>
#include "proto.h"
#include <stack>

const E_Int normalIn[6] = {1,0,1,0,1,0};

// Zero based face and cell
E_Int get_reorient(E_Int face, E_Int cell, E_Int normalIn, mesh *M)
{
  assert(M->owner[face] == cell || M->neigh[face] == cell);
  if (M->neigh[face] == cell && normalIn == 1) return 0;
  else if (M->owner[face] == cell && normalIn == 0) return 0;
  else return 1;
}

static
void build_face_neighbourhood
(
  std::vector<E_Int> &pgs,
  std::vector<E_Int> &xpgs,
  std::vector<E_Int> &neighbour
)
{
  for (size_t i = 0; i < neighbour.size(); i++)
    neighbour[i] = -1;

  std::map<edge,
           std::pair<std::pair<E_Int, E_Int>, std::pair<E_Int, E_Int>>> EM;

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
          search->second.second = std::make_pair(IDX_NONE, IDX_NONE);
      }
    }
  }

  for (auto &elem : EM) {
    E_Int pg0 = elem.second.first.first;
    E_Int n0 = elem.second.first.second;
    E_Int pg1 = elem.second.second.first;
    E_Int n1 = elem.second.second.second;

    if (pg1 == -1 || pg1 == IDX_NONE)
      continue;

    E_Int s0 = xpgs[pg0];
    E_Int s1 = xpgs[pg1];
    
    neighbour[s0 + n0] = pg1;
    neighbour[s1 + n1] = pg0;
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
void colouring(const std::vector<E_Int>& neighbours, const std::vector<E_Int>& xadj,
  std::vector<E_Int>& colours)
{
  E_Int K, Kseed(0), colour(0), nelems(xadj.size()-1);
  std::vector<E_Int> kpool;

  colours.resize(nelems, IDX_NONE);

  while (1) {
    while ((Kseed < nelems) && (colours[Kseed] != IDX_NONE))
      Kseed++;

    if (Kseed >= nelems)
      return;

    kpool.push_back(Kseed);

    while (!kpool.empty()) {
      K = kpool.back();
      kpool.pop_back();

      if (colours[K] != IDX_NONE)
        continue;

      colours[K] = colour;

      E_Int start = xadj[K];
      E_Int end = xadj[K+1];
      E_Int stride = end - start;
      const E_Int *pn = &neighbours[start];
      for (E_Int i = 0; i < stride; i++) {
        E_Int Kn = pn[i];
        if ((Kn != IDX_NONE) && (colours[Kn] == IDX_NONE))
          kpool.push_back(Kn);
      }
    }

    colour++;
  }
}


static
void orient_boundary(mesh *M)
{
  std::vector<E_Int> fadj;
  std::vector<E_Int> xadj(1, 0);
  std::vector<E_Int> efaces;

  // gather all the external faces
  // external faces are faces that belong to only one cell
  std::vector<E_Int> fcount(M->nfaces, 0);
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->xcells[i]; j < M->xcells[i+1]; j++) {
      E_Int face = M->NFACE[j];
      fcount[face] += 1;
    }
  }

  for (E_Int i = 0; i < M->nfaces; i++) {
    assert(fcount[i] == 1 || fcount[i] == 2);
    if (fcount[i] == 1) {
      efaces.push_back(i);
      xadj.push_back(4);
      for (E_Int j = M->xfaces[i]; j < M->xfaces[i+1]; j++)
        fadj.push_back(M->NGON[j]);
    }
  }

  E_Int nefaces = (E_Int)efaces.size();

  for (E_Int i = 0; i < nefaces; i++)
    xadj[i+1] += xadj[i];

  std::vector<E_Int> neighbours(fadj.size(), -1);
  build_face_neighbourhood(fadj, xadj, neighbours);
  // Neighbour i of facej is the face shared edge(pn[i], pn[i+1])

  // remove non-manifoldness
  {
    std::map<edge, E_Int> edge_count;
    edge E;
    for (E_Int i = 0; i < nefaces; i++) {
      E_Int start = xadj[i];
      E_Int end = xadj[i+1];
      E_Int stride = end - start;
      E_Int *pn = &fadj[start];

      for (E_Int j = 0; j < stride; j++) {
        E_Int n0 = pn[j];
        E_Int n1 = pn[(j+1)%stride];
        E.set(n0, n1);
        auto search = edge_count.find(E);
        if (search == edge_count.end())
          edge_count.insert(std::make_pair(E, 1));
        else
          search->second++;
      }
    }

    for (E_Int i = 0; i < nefaces; i++) {
      E_Int start = xadj[i];
      E_Int end = xadj[i+1];
      E_Int stride = end - start;

      E_Int *pn = &fadj[start];
      E_Int *pk = &neighbours[start];

      for (E_Int n = 0; n < stride; n++) {
        E_Int n0 = pn[n];
        E_Int n1 = pn[(n+1)%stride];
        E.set(n0, n1);

        if (edge_count[E] != 2)
          pk[n] = IDX_NONE;
      }
    }
  }

  //E_Int nb_connex = 1;
  std::vector<E_Int> colours;
  colouring(neighbours, xadj, colours);
  //nb_connex = 1 + *std::max_element(colours.begin(), colours.end());

  // faces belonging to the same connex part should have the same colour
  //std::vector<E_Int> colours(M->nfaces);
  //colour_elements(neighbours, colours);

  std::vector<E_Int> processed(nefaces, 0);
  std::vector<E_Int> orient(nefaces, 0);

  // As many iterations as connex bits
  while (1) {
    // look for first unprocessed face
    E_Int seed = 0;
    for (; seed < nefaces; seed++)
      if (!processed[seed])
        break;
    
    if (seed >= nefaces) // all done
      break;
    
    std::stack<E_Int> fpool; // face pool
    fpool.push(seed);

    // is the seed well oriented?
    // TODO(Imad): improve this check
    compute_cell_centers(M);
    E_Int gf = efaces[seed];
    E_Int *pn = &M->NGON[M->xfaces[gf]];
    E_Int n0 = pn[0];
    E_Int n1 = pn[1];
    E_Int n2 = pn[2];
    E_Float *X0 = &M->xyz[3*n0];
    E_Float *X1 = &M->xyz[3*n1];
    E_Float *X2 = &M->xyz[3*n2];
    E_Float n[3]; // face normal
    E_Float e1[3], e2[3]; // edges 10 & 20
    for (E_Int i = 0; i < 3; i++) {
      e1[i] = X1[i] - X0[i];
      e2[i] = X2[i] - X0[i];
    }
    cross(e1, e2, n);
    E_Float *fx = &M->fc[3*gf];
    E_Int own = M->owner[gf];
    assert(own >= 0 && M->neigh[gf] == -1);
    E_Float *cx = &M->cc[3*own];
    E_Float d[3]; // vector from face center to cell center
    for (E_Int i = 0; i < 3; i++)
      d[i] = cx[i] - fx[i];
    // face should be oriented outwards: dot product should be negative
    if (dot(d, n, 3) > 0.)
      orient[seed] = -1;

    while (!fpool.empty()) {
      E_Int face = fpool.top();
      fpool.pop();

      if (processed[face]) {
        continue;
      }
      
      processed[face] = 1;

      E_Int start = xadj[face];
      E_Int end = xadj[face+1];
      E_Int stride = end - start;
      E_Int *pf = &fadj[start];
      
      for (E_Int i = xadj[face]; i < xadj[face+1]; i++) {
        E_Int nei = neighbours[i];
        if (nei == IDX_NONE)
          continue;
        if (processed[nei] == 1)
          continue;
        

        E_Int nstart = xadj[nei];
        E_Int nend = xadj[nei+1];
        E_Int nstride = nend - nstart;
        E_Int *pn = &fadj[nstart];

        // find the shared edge
        E_Int k, l;
        k = l = -1;
        get_boundary(pf, stride, pn, nstride, &k, &l);
        assert(k != -1 && l != -1);

        E_Int ni = pf[k];
        E_Int nj = pf[(k+1)%stride];

        E_Int reverse = 2;
        E_Int ret = get_orientation(pn, nstride, ni, nj, &reverse);
        assert(ret != -1);

        if (orient[face] == -1) reverse = !reverse;
        if (reverse) orient[nei] = -1;
        fpool.push(nei);
      }
    }
  }

  for (E_Int i = 0; i < nefaces; i++) {
    if (orient[i] == -1) {
      E_Int face = efaces[i];
      E_Int start = M->xfaces[face];
      E_Int stride = M->xfaces[face+1] - start;
      E_Int *pn = &M->NGON[start];
      std::reverse(pn+1, pn+stride);
      //std::swap(pn[1], pn[3]);
    }
  }
}

void topo_init_mesh(mesh *M)
{
  //reorient_skin(M);
  orient_boundary(M);
  build_own_nei(M);
  reorder_hexa(M);
}

void reorder_hexa(mesh *M)
{
  E_Int common[4];
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int *pf = &M->NFACE[6*i];
    E_Int bot = pf[0];
    E_Int *pn = &M->NGON[4*bot];
    
    E_Int map[4];
    for (E_Int j = 0; j < 4; j++) map[j] = pn[j];
    E_Int reorient = get_reorient(bot, i, normalIn[0], M);
    if (reorient) std::swap(map[1], map[3]);
    
    E_Int top, lft, rgt, fro, bck;
    top = -1;
    lft = -1;
    rgt = -1;
    fro = -1;
    bck = -1;
    
    for (E_Int j = 1; j < 6; j++) {
      E_Int face = pf[j];

      for (E_Int k = 0; k < 4; k++)
        common[k] = 0;

      pn = &M->NGON[4*face];

      for (E_Int k = 0; k < 4; k++) {
        E_Int point = pn[k];
        // look for point in map
        for (E_Int l = 0; l < 4; l++) {
          if (map[l] == point)
            common[l] = 1;
        }
      }
      if (common[0] && common[3]) lft = face;
      else if (common[0] && common[1]) fro = face;
      else if (common[1] && common[2]) rgt = face;
      else if (common[3] && common[2]) bck = face;
      else top = face;
    }
    pf[1] = top;
    pf[2] = lft;
    pf[3] = rgt;
    pf[4] = fro;
    pf[5] = bck;
  }
}

/*
void reorient_skin(mesh *M)
{
  compute_cell_centers(M);
  
  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = neigh[i];

    if (nei != -1)
      continue;

    // compute face normal
    E_Int start = M->xfaces[i];
    E_Int *pn = &M->NGON[start];

    E_Float e0[3], e1[3];

    E_Float *X0 = &M->xyz[3*pn[0]];
    E_Float *X1 = &M->xyz[3*pn[1]];
    E_Float *X2 = &M->xyz[3*pn[2]];
    
    for (E_Int j = 0; j < 3; j++) {
      e0[j] = X1[j] - X0[j];
      e1[j] = X2[j] - X0[j];
    }
    
    E_Float n[3];
    cross(e0, e1, n);

    // get owner's cell center
    E_Int own = owner[i];
    E_Float *cx = &M->cc[3*own];
    E_Float *fx = &M->fc[3*i];

    // compute vector linking face center to cell center
    E_Float d[3];
    for (E_Int j = 0; j < 3; j++)
      d[j] = cx[j] - fx[j];

    // normal should be pointing outwards: negative dot product
    if (dot(d, n, 3) > 0.) {
      //std::reverse(M->NGON+start+1, M->NGON+start+stride);
      std::swap(pn[1], pn[3]);
    }
  }
}
*/

static
void reversi_connex(std::vector<E_Int> &pgs, std::vector<E_Int> &xpgs,
  std::vector<E_Int> &neighbours, E_Int kseed, std::vector<E_Int> &orient)
{
  std::vector<E_Int> cpool;
  cpool.push_back(kseed);

  std::vector<E_Int> processed(xpgs.size()-1, 0);

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
void build_cell_neighbourhood(mesh *M, std::vector<E_Int>& neighbours)
{
  neighbours.resize(M->xcells[M->ncells], -1);

  std::vector<E_Int> neigh(M->nfaces, -1);
  E_Int count = 0;

  while (count++ != 2) {
    for (E_Int i = 0; i < M->ncells; i++) {
      E_Int start = M->xcells[i];
      E_Int end = M->xcells[i+1];
      E_Int stride = end - start;
      E_Int *pf = &M->NFACE[start];
      for (E_Int j = 0; j < stride; j++) {
        E_Int face = pf[j];
        E_Int &nei = neigh[face];
        E_Int &Kn = neighbours[6*i+j];
        if (nei != -1 && nei != i)
          Kn = nei;
        neigh[face] = i;
      }
    }
  }
}

void build_own_nei(mesh *M)
{
  std::vector<E_Int> neighbours;
  build_cell_neighbourhood(M, neighbours);

  // assumes external faces have been properly oriented outwards
  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  for (E_Int i = 0; i < M->nfaces; i++) {
    owner[i] = -1;
    neigh[i] = -1;
  }

  std::vector<E_Int> exPH(M->ncells, -1);
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->xcells[i]; j < M->xcells[i+1]; j++) {
      if (neighbours[j] == -1) {
        owner[M->NFACE[j]] = i;
        exPH[i] = M->NFACE[j]+1;
        break;
      }
    }
  }

  // look for first external cell
  std::vector<E_Int> processed(M->ncells, 0);

  E_Int nconnex = 0;
  
  E_Int seed = 0;

  while (1) {
    while ((seed < M->ncells) && (processed[seed] || exPH[seed] == -1))
      seed++;

    if (seed >= M->ncells)
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
      E_Int start = M->xcells[cell];
      E_Int end = M->xcells[cell+1];
      E_Int stride = end - start;
      std::vector<E_Int> oids;
      std::vector<E_Int> orient(stride, 1);

      std::vector<E_Int> pgs;
      std::vector<E_Int> xpgs;
      xpgs.push_back(0);
      for (E_Int i = start; i < end; i++) {
        E_Int face = M->NFACE[i];
        for (E_Int j = M->xfaces[face]; j < M->xfaces[face+1]; j++)
          pgs.push_back(M->NGON[j]);
        xpgs.push_back(M->xfaces[face+1]-M->xfaces[face]);
        oids.push_back(face+1);
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
      reversi_connex(pgs, xpgs, PGneighbours, iref, orient);

      // set the owner and neighbour of the faces
      E_Int *pn = &neighbours[start];
      E_Int *pf = &M->NFACE[start];
      for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i];
        E_Int nei = pn[i];

        M->owner[face] = cell;
        M->neigh[face] = nei;

        if (nei == -1)
          continue;

        // set the reference face for neighbour
        exPH[nei] = -(face+1);

        if (orient[i] == -1) {
          std::swap(M->owner[face], M->neigh[face]);
          exPH[nei] = face+1;
        }

        if (!processed[nei])
          cpool.push_back(nei);
      }
    }
  }
}

E_Int get_pos(E_Int e, E_Int *pn, E_Int size)
{
  for (E_Int i = 0; i < size; i++) {
    if (pn[i] == e)
      return i;
  }
  return -1;
}

void right_shift(E_Int *pn, E_Int pos, E_Int size)
{
  E_Int tmp[10];
  assert(size <= 10);
  for (E_Int i = 0; i < size; i++) tmp[i] = pn[i];
  for (E_Int i = 0; i < size; i++)
    pn[i] = tmp[(i+pos)%size];
}

void order_quad(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0)
{
  for (E_Int i = 0; i < 4; i++) local[i] = pn[i];
  right_shift(local, i0, 4);
  if (reorient) std::swap(local[1], local[3]);
}

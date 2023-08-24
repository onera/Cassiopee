#include <stdio.h>
#include "proto.h"
//#include <algorithm>

const E_Int normalIn[6] = {1,0,1,0,1,0};

// Zero based face and cell
E_Int get_reorient(E_Int face, E_Int cell, E_Int normalIn, mesh *M)
{
  if (M->neigh[face] == cell && normalIn == 1) return 0;
  else if (M->owner[face] == cell && normalIn == 0) return 0;
  else return 1;
}

void topo_init_mesh(mesh *M)
{
  reorient_skin(M);
  build_own_nei(M);
  reorder_hexa(M);
}

void reorder_hexa(mesh *M)
{
  E_Int map[4];
  E_Int common[4];
  E_Int face, *pf, *pn;
  E_Int bot, top, lft, rgt, fro, bck;
  for (E_Int i = 0; i < M->ncells; i++) {
    pf = &M->NFACE[6*i];
    bot = pf[0];
    pn = &M->NGON[4*bot];
    for (E_Int j = 0; j < 4; j++) map[j] = pn[j];
    E_Int reorient = get_reorient(bot, i, normalIn[0], M);
    if (reorient) std::swap(map[1], map[3]);
    top = lft = rgt = fro = bck = -1;
    for (E_Int j = 1; j < 6; j++) {
      for (E_Int k = 0; k < 4; k++)
        common[k] = 0;
      face = pf[j];
      pn = &M->NGON[4*face];
      for (E_Int k = 0; k < 4; k++) {
        for (E_Int l = 0; l < 4; l++) {
          if (map[k] == pn[l])
            common[k] = 1;
        }
      }
      if (common[0] && common[3]) lft = face;
      else if (common[0] && common[1]) fro = face;
      else if (common[1] && common[2]) rgt = face;
      else if (common[3] && common[2]) bck = face;
      else top = face;
    }
    assert(!(top == -1 || lft == -1 || rgt == -1 || fro == -1 || bck == -1));
    pf[1] = top;
    pf[2] = lft;
    pf[3] = rgt;
    pf[4] = fro;
    pf[5] = bck;
  }
}

void reorient_skin(mesh *M)
{
  compute_cell_centers(M);
  
  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = neigh[i];

    assert(nei < M->ncells);
    
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
      E_Int stride = M->xfaces[i+1] - start;
      std::reverse(M->NGON+start+1, M->NGON+start+stride);
      //std::swap(pn[1], pn[3]);
    }
  }
}

static
void build_face_neighbourhood
(
  std::vector<E_Int> &pgs,
  std::vector<E_Int> &xpgs,
  mesh *M,
  std::vector<E_Int> &neighbour
)
{
  for (size_t i = 0; i < neighbour.size(); i++)
    neighbour[i] = -1;

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
      } else {
        search->second.second = std::make_pair(i, j);
      }
    }
  }

  for (auto &elem : EM) {
    E_Int pg0 = elem.second.first.first;
    E_Int n0 = elem.second.first.second;
    E_Int pg1 = elem.second.second.first;
    E_Int n1 = elem.second.second.second;

    E_Int s0 = xpgs[pg0];
    E_Int s1 = xpgs[pg1];
    
    neighbour[s0 + n0] = pg1;
    neighbour[s1 + n1] = pg0;
  }
}

static
void get_boundary(E_Int *pn0, E_Int s0, E_Int *pn1, E_Int s1, E_Int *m, E_Int *n)
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
E_Int get_orientation(E_Int *pn, E_Int stride, E_Int ni, E_Int nj, E_Int *same_orient)
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
void reversi_connex
(
  std::vector<E_Int> &pgs,
  std::vector<E_Int> &xpgs,
  std::vector<E_Int> &neighbours,  
  E_Int kseed,
  std::vector<E_Int> &orient
)
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
      assert(k != -1);
      assert(l != -1);

      E_Int ni = pf[k];
      E_Int nj = pf[(k+1)%stride];

      E_Int reverse = 2;
      E_Int ret = get_orientation(pnn, sn, ni, nj, &reverse);
      assert(ret != -1);

      if (orient[K] == -1) {
        if (reverse == 0) reverse = 1;
        else if (reverse == 1) reverse = 0;
        else assert(0);
      }

      if (reverse)
        orient[nei] = -1;

      cpool.push_back(nei);
    }
  }
}

static
void build_cell_neighbourhood(mesh *M, std::vector<E_Int>& neighbours)
{
  neighbours.resize(M->xcells[M->ncells]);
  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->xcells[i]; j < M->xcells[i+1]; j++) {
      E_Int face = M->NFACE[j];
      assert(face >= 0 && face < M->nfaces);
      E_Int own = owner[face];
      E_Int nei = neigh[face];
      assert(own < M->ncells);
      assert(nei < M->ncells);

      if (i == own) neighbours[j] = nei;
      else neighbours[j] = own;
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
  E_Int seed = 0;
  while ((seed < M->ncells) && (processed[seed] || exPH[seed] == 0))
    seed++;

  assert(seed < M->ncells);
  
  std::vector<E_Int> cpool;

  cpool.push_back(seed);

  E_Int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  while (!cpool.empty()) {
    E_Int cell = cpool.back();
    cpool.pop_back();

    assert(cell < M->ncells);

    if (processed[cell])
      continue;

    processed[cell] = 1;

    // build faces neighbourhood based on shared edges' nodes order
    E_Int start = M->xcells[cell];
    E_Int end = M->xcells[cell+1];
    E_Int stride = end - start;
    assert(stride == 6);
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

    // Note(Imad): REMOVE THESE ASSERTIONS
    assert(pgs.size() == 24);
    assert(xpgs.size() == 7);

    std::vector<E_Int> PGneighbours(pgs.size());
    build_face_neighbourhood(pgs, xpgs, M, PGneighbours);

    E_Int revers = 0;

    // reference face is the external face
    E_Int PGref = exPH[cell];
    assert(PGref != 0);
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

    assert(orient.size() == oids.size());

    assert(iref != -1);

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

      assert(nei < M->ncells);

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

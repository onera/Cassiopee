#include "Proto.h"
#include <cassert>
#include <stack>

#define MINREF -10

E_Int Get_pos(E_Int e, E_Int *pn, E_Int size)
{
  for (E_Int i = 0; i < size; i++) {
    if (pn[i] == e)
      return i;
  }
  return -1;
}

void Right_shift(E_Int *pn, E_Int pos, E_Int size)
{
  E_Int tmp[10];
  assert(size <= 10);
  for (E_Int i = 0; i < size; i++) tmp[i] = pn[i];
  for (E_Int i = 0; i < size; i++)
    pn[i] = tmp[(i+pos)%size];
}

static
void hessian_to_metric(E_Float *H, E_Int ncells)
{
  E_Float L[3], v0[3], v1[3], v2[3];
  E_Float L0, L1, L2, *pt;
  for (E_Int i = 0; i < ncells; i++) {
    pt = &H[6*i];
    K_LINEAR::sym3mat_eigen(pt, L, v0, v1, v2);
    L0 = sqrt(fabs(L[0])); L1 = sqrt(fabs(L[1])); L2 = sqrt(fabs(L[2]));
    pt[0] = L0*v0[0]*v0[0] + L1*v1[0]*v1[0] + L2*v2[0]*v2[0];
    pt[1] = L0*v0[0]*v0[1] + L1*v1[0]*v1[1] + L2*v2[0]*v2[1];
    pt[2] = L0*v0[0]*v0[2] + L1*v1[0]*v1[2] + L2*v2[0]*v2[2];
    pt[3] = L0*v0[1]*v0[1] + L1*v1[1]*v1[1] + L2*v2[1]*v2[1];
    pt[4] = L0*v0[1]*v0[2] + L1*v1[1]*v1[2] + L2*v2[1]*v2[2];
    pt[5] = L0*v0[2]*v0[2] + L1*v1[2]*v1[2] + L2*v2[2]*v2[2];
  }
}

static
void make_pdirs_pyra(E_Int cell, AMesh *M, pDirs &Dirs)
{
  // Setup bot
  E_Int *pf = get_facets(cell, M->nface, M->indPH);
  E_Int face = pf[0];
  E_Int *pn = get_facets(face, M->ngon, M->indPG);
  E_Int NODES[5] = {pn[0], pn[1], pn[2], pn[3], -1,};
  E_Int reorient = get_reorient(face, cell, normalIn_Py[0], M);
  if (reorient) std::swap(NODES[1], NODES[3]);

  E_Int lft = pf[1];
  pn = &M->ngon[M->indPG[lft]];

  for (E_Int i = 0; i < 3; i++) {
    E_Int pt = pn[i];

    E_Int found = 0;
    // pt must not be in n
    for (E_Int j = 0; j < 4; j++) {
      if (pt == NODES[j]) {
        found = 1;
        break;
      }
    }

    if (!found) {
      NODES[4] = pt;
      break;
    }
  }

  assert(NODES[4] != NODES[0]);
  assert(NODES[4] != NODES[1]);
  assert(NODES[4] != NODES[2]);
  assert(NODES[4] != NODES[3]);

  E_Float *f0, *f1;

  // Face 0 (bot) <-> n4
  f0 = &M->fc[3*pf[0]];
  Dirs.I[0] = M->x[NODES[4]] - f0[0];
  Dirs.I[1] = M->y[NODES[4]] - f0[1];
  Dirs.I[2] = M->z[NODES[4]] - f0[2];

  // LFT - RGT
  f0 = &M->fc[3*pf[1]];
  f1 = &M->fc[3*pf[2]];
  for (E_Int i = 0; i < 3; i++) Dirs.J[i] = f1[i] - f0[i];

  // FRO - BCK
  f0 = &M->fc[3*pf[3]];
  f1 = &M->fc[3*pf[4]];
  for (E_Int i = 0; i < 3; i++) Dirs.K[i] = f1[i] - f0[i];
}

static
void make_pdirs_penta(E_Int cell, AMesh *M, pDirs &Dirs)
{
  // Setup bot
  E_Int *pf = get_facets(cell, M->nface, M->indPH);
  E_Int face = pf[0];
  E_Int *pn = get_facets(face, M->ngon, M->indPG);
  E_Int NODES[6] = {pn[0], pn[1], pn[2], -1, -1, -1};
  E_Int reorient = get_reorient(face, cell, normalIn_Pe[0], M);
  if (reorient) std::swap(NODES[1], NODES[2]);

  // Deduce left and right
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  E_Int i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[2], M);
  E_Int local[4];
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[2]);
  NODES[5] = local[2];
  NODES[3] = local[3];

  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[3], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[3] == NODES[3]);
  NODES[4] = local[2];

  E_Float *fc, *fcc;
  E_Float ec[3];

  // First dir: Edge(0,3) - BCK
  ec[0] = 0.5*(M->x[NODES[0]]+M->x[NODES[3]]);
  ec[1] = 0.5*(M->y[NODES[0]]+M->y[NODES[3]]);
  ec[2] = 0.5*(M->z[NODES[0]]+M->z[NODES[3]]);
  fc = &M->fc[3*pf[4]];
  Dirs.I[0] = ec[0] - fc[0];
  Dirs.I[1] = ec[1] - fc[1];
  Dirs.I[2] = ec[2] - fc[2];
  
  // Second dir: Edge(2,5) - RGT
  ec[0] = 0.5*(M->x[NODES[2]]+M->x[NODES[5]]);
  ec[1] = 0.5*(M->y[NODES[2]]+M->y[NODES[5]]);
  ec[2] = 0.5*(M->z[NODES[2]]+M->z[NODES[5]]);
  fc = &M->fc[3*pf[3]];
  Dirs.J[0] = ec[0] - fc[0];
  Dirs.J[1] = ec[1] - fc[1];
  Dirs.J[2] = ec[2] - fc[2];

  // Third dir: Edge(1,4) - LFT
  ec[0] = 0.5*(M->x[NODES[1]]+M->x[NODES[4]]);
  ec[1] = 0.5*(M->y[NODES[1]]+M->y[NODES[4]]);
  ec[2] = 0.5*(M->z[NODES[1]]+M->z[NODES[4]]);
  fc = &M->fc[3*pf[2]];
  Dirs.K[0] = ec[0] - fc[0];
  Dirs.K[1] = ec[1] - fc[1];
  Dirs.K[2] = ec[2] - fc[2];

  // Fourth dir: BOT - TOP
  fc =  &M->fc[3*pf[0]];
  fcc = &M->fc[3*pf[1]];
  Dirs.L[0] = fc[0] - fcc[0];
  Dirs.L[1] = fc[1] - fcc[1];
  Dirs.L[2] = fc[2] - fcc[2];
}

static
void make_pdirs_tetra(E_Int cell, AMesh *M, pDirs &Dirs)
{
  // Make nodes
  E_Int *pf = &M->nface[M->indPH[cell]];
  E_Int bot = pf[0];
  E_Int *pn = &M->ngon[M->indPG[bot]];
  E_Int n[4] = {pn[0], pn[1], pn[2], -1};
  if (get_reorient(bot, cell, normalIn_T[0], M)) std::swap(n[1], n[2]);

  E_Int lft = pf[1];
  pn = &M->ngon[M->indPG[lft]];

  for (E_Int i = 0; i < 3; i++) {
    E_Int pt = pn[i];

    E_Int found = 0;
    // pt must not be in n
    for (E_Int j = 0; j < 3; j++) {
      if (pt == n[j]) {
        found = 1;
        break;
      }
    }

    if (!found) {
      n[3] = pt;
      break;
    }
  }

  assert(n[3] != n[0]);
  assert(n[3] != n[1]);
  assert(n[3] != n[2]);

  E_Float *fc;

  // Face 0 (bot) <-> n3
  fc = &M->fc[3*bot];
  Dirs.I[0] = M->x[n[3]] - fc[0];
  Dirs.I[1] = M->y[n[3]] - fc[1];
  Dirs.I[2] = M->z[n[3]] - fc[2];

  // Face 1 (lft) <-> n1
  fc = &M->fc[3*lft];
  Dirs.J[0] = M->x[n[1]] - fc[0];
  Dirs.J[1] = M->y[n[1]] - fc[1];
  Dirs.J[2] = M->z[n[1]] - fc[2];

  // Face 2 (rgt) <-> n2
  fc = &M->fc[3*pf[2]];
  Dirs.K[0] = M->x[n[2]] - fc[0];
  Dirs.K[1] = M->y[n[2]] - fc[1];
  Dirs.K[2] = M->z[n[2]] - fc[2];

  // Face 3 (bck) <-> n0
  fc = &M->fc[3*pf[3]];
  Dirs.L[0] = M->x[n[0]] - fc[0];
  Dirs.L[1] = M->y[n[0]] - fc[1];
  Dirs.L[2] = M->z[n[0]] - fc[2];
}

static
void make_pdirs_hexa(E_Int cell, AMesh *M, pDirs &Dirs)
{
  E_Int *pf = &M->nface[M->indPH[cell]];
  E_Float *f0, *f1;

  f0 = &M->fc[3*pf[2]];
  f1 = &M->fc[3*pf[3]];
  for (E_Int i = 0; i < 3; i++) Dirs.I[i] = f1[i] - f0[i];

  f0 = &M->fc[3*pf[4]];
  f1 = &M->fc[3*pf[5]];
  for (E_Int i = 0; i < 3; i++) Dirs.J[i] = f1[i] - f0[i];

  f0 = &M->fc[3*pf[0]];
  f1 = &M->fc[3*pf[1]];
  for (E_Int i = 0; i < 3; i++) Dirs.K[i] = f1[i] - f0[i];
}

static
void compute_principal_vecs(AMesh *M, std::vector<pDirs> &Dirs)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    switch (M->cellTree[i]->type) {
      case TETRA:
        make_pdirs_tetra(i, M, Dirs[i]);
        break;
      case PENTA:
        make_pdirs_penta(i, M, Dirs[i]);
        break;
      case PYRA:
        make_pdirs_pyra(i, M, Dirs[i]);
        break;
      case HEXA:
        make_pdirs_hexa(i, M, Dirs[i]);
        break;
      default:
        assert(0);
        break;
    }
  }
}

// We refine if metric length of one dir is bigger than ref threshold
// We unrefine if metric length of all dirs is smaller than unref threshold
static
void make_ref_data_hexa(AMesh *M, E_Float *pH, const pDirs &Dirs,
  E_Int *pref)
{
  E_Float L0, L1, L2, dd[3];

  K_MATH::sym3mat_dot_vec(pH, Dirs.I, dd);
  L0 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.J, dd);
  L1 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.K, dd);
  L2 = K_MATH::norm(dd, 3);

  //printf("%f %f %f\n", L0, L1, L2);

  if      (L0 >= M->ref_Tr) { pref[0] = 1; }
  else if (L1 >= M->ref_Tr) { pref[1] = 1; }
  else if (L2 >= M->ref_Tr) { pref[2] = 1; }

  if (pref[0] || pref[1] || pref[2]) {
    pref[0] = pref[1] = pref[2] = 1;
    return;
  }

  if (L0 < M->unref_Tr && L1 < M->unref_Tr && L2 < M->unref_Tr) {
    pref[0] = -1;
    pref[1] = -1;
    pref[2] = -1;
  }
}

static
void make_ref_data_tetra(AMesh *M, E_Float *pH, const pDirs &Dirs,
  E_Int *pref)
{
  E_Float L0, L1, L2, L3, dd[3];

  K_MATH::sym3mat_dot_vec(pH, Dirs.I, dd);
  L0 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.J, dd);
  L1 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.K, dd);
  L2 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.L, dd);
  L3 = K_MATH::norm(dd, 3);

  if (L0 >= M->ref_Tr || L1 >= M->ref_Tr ||
      L2 >= M->ref_Tr || L3 >= M->ref_Tr) {
    pref[0] = 1; pref[1] = 1; pref[2] = 1;
    return;
  }

  if (L0 <= M->unref_Tr && L1 <= M->unref_Tr &&
      L2 <= M->unref_Tr && L3 <= M->unref_Tr) {
    pref[0] = -1;
    pref[1] = -1;
    pref[2] = -1;
  }
}

static
void make_ref_data_penta(AMesh *M, E_Float *pH, const pDirs &Dirs,
  E_Int *pref)
{
  E_Float L0, L1, L2, L3, dd[3];

  K_MATH::sym3mat_dot_vec(pH, Dirs.I, dd);
  L0 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.J, dd);
  L1 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.K, dd);
  L2 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.L, dd);
  L3 = K_MATH::norm(dd, 3);

  if (L0 >= M->ref_Tr || L1 >= M->ref_Tr ||
      L2 >= M->ref_Tr || L3 >= M->ref_Tr) {
    pref[0] = 1; pref[1] = 1; pref[2] = 1;
    return;
  }

  if (L0 <= M->unref_Tr && L1 <= M->unref_Tr &&
      L2 <= M->unref_Tr && L3 <= M->unref_Tr) {
    pref[0] = -1;
    pref[1] = -1;
    pref[2] = -1;
  }
}

static
void make_ref_data_pyra(AMesh *M, E_Float *pH, const pDirs &Dirs,
  E_Int *pref)
{
  E_Float L0, L1, L2, dd[3];

  K_MATH::sym3mat_dot_vec(pH, Dirs.I, dd);
  L0 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.J, dd);
  L1 = K_MATH::norm(dd, 3);
  K_MATH::sym3mat_dot_vec(pH, Dirs.K, dd);
  L2 = K_MATH::norm(dd, 3);

  //printf("%f %f %f\n", L0, L1, L2);

  if      (L0 >= M->ref_Tr) { pref[0] = 1; }
  else if (L1 >= M->ref_Tr) { pref[1] = 1; }
  else if (L2 >= M->ref_Tr) { pref[2] = 1; }

  if (pref[0] || pref[1] || pref[2]) {
    pref[0] = pref[1] = pref[2] = 1;
    return;
  }

  if (L0 < M->unref_Tr && L1 < M->unref_Tr && L2 < M->unref_Tr) {
    pref[0] = -1;
    pref[1] = -1;
    pref[2] = -1;
  }
}

static
void make_ref_data(AMesh *M, E_Float *H, const std::vector<pDirs> &Dirs,
  std::vector<E_Int> &partialRefData)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int type = M->cellTree[i]->type;
    E_Float *pH = &H[6*i];
    E_Int *pref = &partialRefData[3*i];
    switch (type) {
      case HEXA:
        make_ref_data_hexa(M, pH, Dirs[i], pref);
        break;
      case TETRA:
        make_ref_data_tetra(M, pH, Dirs[i], pref);
        break;
      case PENTA:
        make_ref_data_penta(M, pH, Dirs[i], pref);
        break;
      case PYRA:
        make_ref_data_pyra(M, pH, Dirs[i], pref);
        break;
      default:
        assert(0);
        break;
    }
  }
}

static
void get_neighbours(E_Int cell, E_Int face, AMesh *M, std::vector<E_Int> &neis)
{
  E_Int nei = get_neighbour(cell, face, M);
  if (nei == -1) return;

  neis.push_back(nei);

  assert(M->cellTree[cell]->level == 0);
  assert(M->cellTree[nei]->level == 0);

  if (M->cellTree[cell]->level >= M->cellTree[nei]->level) return;
  puts("never here");

  E_Int nchildren = M->faceTree[face]->nchildren;
  E_Int *children = M->faceTree[face]->children;
  assert(nchildren == 3);
  for (E_Int i = 0; i < nchildren; i++) {
    E_Int small_face = children[i];
    E_Int nei = get_neighbour(cell, small_face, M);
    if (nei != -1) neis.push_back(nei);
  }
}

static
void smooth_ref_data(AMesh *M)
{
  std::stack<E_Int> stk;
  for (E_Int i = 0; i < M->ncells; i++) {
    if (M->ref_data[3*i] != 0)
      stk.push(i);
  }

  while (!stk.empty()) {
    E_Int cell = stk.top();
    stk.pop();

    std::vector<E_Int> neis;
    
    for (E_Int i = M->indPH[cell]; i < M->indPH[cell+1]; i++) {
      E_Int face = M->nface[i];
      get_neighbours(cell, face, M, neis);
    }

    E_Int incr_cell = M->ref_data[3*cell];
    incr_cell += M->cellTree[cell]->level;
    
    E_Int *pr = &M->ref_data[3*cell];

    for (size_t i = 0; i < neis.size(); i++) {
      E_Int nei = neis[i];
      E_Int incr_nei = M->ref_data[3*nei];
      incr_nei += M->cellTree[nei]->level;
      if (abs(incr_nei - incr_cell) <= 1) continue;

      E_Int cell_to_mod = incr_cell > incr_nei ? nei : cell;

      E_Int *pr_mod = &M->ref_data[3*cell_to_mod];
      if (pr_mod[0] >= 0) {
        pr_mod[0] += 1;
        pr_mod[1] += 1;
        pr_mod[2] += 1;
        stk.push(cell_to_mod);
      } else {
        E_Int parent = M->cellTree[cell_to_mod]->parent;
        E_Int *children = M->cellTree[parent]->children;
        E_Int nchildren = M->cellTree[parent]->nchildren;
        
        E_Int *prp = &M->ref_data[3*parent];
        prp[0] += 1; prp[1] += 1; prp[2] += 1;
        stk.push(parent);
        
        for (E_Int j = 0; j < nchildren; j++) {
          E_Int *prc = &M->ref_data[3*children[j]];
          prc[0] += 1; prc[1] += 1; prc[2] += 1;
          stk.push(children[j]);
        }
      }
    }
  }
}


void compute_ref_data(AMesh *M, E_Float **fields, E_Int nfields)
{
  // Init ref data
  M->ref_data = (int *)XMALLOC(3*M->ncells * sizeof(int));
  for (E_Int i = 0; i < 3*M->ncells; i++) M->ref_data[i] = MINREF;

  // Make cells principal directions
  std::vector<pDirs> Dirs(M->ncells);
  compute_principal_vecs(M, Dirs);

  for (E_Int i = 0; i < nfields; i++) {
    E_Float *H = compute_hessian(M, fields[i]);
    hessian_to_metric(H, M->ncells);

    std::vector<int> partialRefData(3*M->ncells, 0);
    make_ref_data(M, H, Dirs, partialRefData);

    for (E_Int j = 0; j < M->ncells; j++) {
      E_Int *pr = &M->ref_data[3*j];
      E_Int *ppr = &partialRefData[3*j];
      for (E_Int k = 0; k < 3; k++)
        pr[k] = std::max(pr[k], ppr[k]);
    }

    XFREE(H);
  }

  // Fix unrefinement data
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int *pr = &M->ref_data[3*i];
    if (pr[0] >= 0) continue;

    // For a cell to be unrefined, all its siblings must be tagged for
    // unrefinement
    E_Int parent = M->cellTree[i]->parent;
    if (parent == i) {
      pr[0] = pr[1] = pr[2] = 0;
      continue;
    }

    E_Int nchildren = M->cellTree[parent]->nchildren;
    E_Int *children = M->cellTree[parent]->children;
    assert(children);
    
    for (E_Int j = 0; j < nchildren; j++) {
      E_Int child = children[j];
      E_Int *prc = &M->ref_data[3*child];
      if (prc[0] < 0) prc[0] = prc[1] = prc[2] = 0;
    }
  }

  // Smooth out refinement data
  //smooth_ref_data(M);
}

#include "Proto.h"
#include "../common/mem.h"
#include <cassert>

#define MINREF -10

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
void make_pdirs_penta(E_Int cell, AMesh *M, pDirs &Dirs)
{}

static
void make_pdirs_pyra(E_Int cell, AMesh *M, pDirs &Dirs)
{}

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

static
void smooth_ref_data(AMesh *M)
{}

void compute_ref_data(AMesh *M, E_Float **fields, E_Int nfields)
{
  // Init ref data
  M->ref_data = (short *)XMALLOC(3*M->ncells * sizeof(short));
  for (E_Int i = 0; i < 3*M->ncells; i++) M->ref_data[i] = MINREF;

  // Make cells principal directions
  std::vector<pDirs> Dirs(M->ncells);
  compute_principal_vecs(M, Dirs);

  for (E_Int i = 0; i < nfields; i++) {
    E_Float *H = compute_hessian(M, fields[i]);
    hessian_to_metric(H, M->ncells);

    /*
    std::vector<std::array<short, 3>> partialRefData(M->ncells);
    make_ref_data(M, H, D_, partialRefData);

    if (i > 0) {
      for (E_Int j = 0; j < M->ncells; j++) {
        E_Int *pr = &M->ref_data[3*j];
        for (E_Int k = 0; k < 3; k++)
          pr[k] = std::max(pr[k], partialRefData[k]);
      } 
    }
    */

    XFREE(H);
  }

  // Smooth out refinement data
  smooth_ref_data(M);
}

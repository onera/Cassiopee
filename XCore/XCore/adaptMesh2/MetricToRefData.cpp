#include "Proto.h"
#include <stack>

#define MINREF -10

static
void make_ref_data_hexa(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{
  E_Float L0, L1, L2, dd[3];
  //E_Float l[3];


  K_MATH::sym3mat_dot_vec(pM, Dirs.I, dd);
  L0 = K_MATH::norm(dd, 3);

  K_MATH::sym3mat_dot_vec(pM, Dirs.J, dd);
  L1 = K_MATH::norm(dd, 3);

  K_MATH::sym3mat_dot_vec(pM, Dirs.K, dd);
  L2 = K_MATH::norm(dd, 3);

  M->ref_data[cell] = 0;

  if (L0 >= M->Tr || L1 >= M->Tr || L2 >= M->Tr) {
    M->ref_data[cell] = 1;
  } else if (L0 < M->Tu && L1 < M->Tu && L2 < M->Tu) {
    if (M->cellTree->level(cell) != 0)
      M->ref_data[cell] = -1;
  }
}

static
void smooth_ref_data(AMesh *M)
{
  //printf("init min ref: %d\n", *std::min_element(M->ref_data, M->ref_data+M->ncells));
  //printf("init max ref: %d\n", *std::max_element(M->ref_data, M->ref_data+M->ncells));

  // Test 2 strategies:
  // 1 - first eliminate all impossible unrefinement, then smooth
  // 2 - smooth, then eliminate impossible unrefinement

  // Note(Imad): I think it should be recursive

  E_Int iter = 0;
  E_Int changed = 1;

  puts("Smoothing ref data");

  while (changed) {
    printf("    Iter %d\n", iter++);

    // Eliminate impossible refinement

    std::vector<E_Int> tmp_ref_data(M->ncells, 0);

    for (E_Int i = 0; i < M->ncells; i++) {
      if (M->ref_data[i] > 0)
        tmp_ref_data[i] = M->ref_data[i];
    }

    for (E_Int i = 0; i < M->ncells; i++) {
      Children *children = M->cellTree->children(i);

      if (children == NULL) continue;

      E_Int skip = 0;

      for (E_Int j = 1; j < children->n; j++) {
        E_Int child = children->pc[j];

        if (M->cellTree->children(child)) {
          skip = 1;
          break;
        }
      }

      if (skip) continue;

      // This is a leaf node
      // To unrefine, all siblings must have negative ref data
      for (E_Int j = 0; j < children->n; j++) {
        E_Int child = children->pc[j];
        if (M->ref_data[child] >= 0) {
          skip = 1;
          break;
        }
      }

      if (skip) {
        continue;
      }

      for (E_Int j = 0; j < children->n; j++) {
        E_Int child = children->pc[j];
        assert(M->ref_data[child] < 0);

        tmp_ref_data[child] = M->ref_data[child];
      }
    }
    
    for (E_Int i = 0; i < M->ncells; i++) M->ref_data[i] = tmp_ref_data[i];

    E_Int cells_to_agglo = 0;
    for (E_Int i = 0; i < M->ncells; i++) {
      if (M->ref_data[i] < 0) cells_to_agglo += 1;
    }

    assert(cells_to_agglo % 8 == 0);


    std::stack<E_Int> stk;
    for (E_Int i = 0; i < M->ncells; i++) {
      if (M->ref_data[i] != 0) stk.push(i);
    }

    for (E_Int i = 0; i < M->ncells; i++) {
      M->ref_data[i] += M->cellTree->level(i);
    }

    changed = 0;

    while (!stk.empty()) {
      E_Int cell = stk.top();
      stk.pop();

      E_Int nf = -1;
      E_Int pf[24];
      get_full_cell(cell, M, nf, pf);

      std::vector<E_Int> neis(nf);
      for (E_Int i = 0; i < nf; i++) {
        neis[i] = get_neighbour(cell, pf[i], M);
      }

      E_Int incr_cell = M->ref_data[cell];
      //incr_cell += M->cellTree->level(cell);

      for (E_Int i = 0; i < nf; i++) {
        E_Int nei = neis[i];
        if (nei == -1) continue;

        E_Int incr_nei = M->ref_data[nei];
        //incr_nei += M->cellTree->level(nei);

        E_Int diff = abs(incr_nei - incr_cell);

        if (diff <= 1) continue;

        changed = 1;

        E_Int cell_to_mod = incr_cell > incr_nei ? nei : cell;

        E_Int &rmod = M->ref_data[cell_to_mod];

        rmod += 1;
        
        //rmod += diff - 1;

        stk.push(cell_to_mod);
      }
    }

    for (E_Int i = 0; i < M->ncells; i++) {
      M->ref_data[i] -= M->cellTree->level(i);
    }

    if (changed == 0) break;
    
    //printf("Min ref: %d\n", *std::min_element(M->ref_data, M->ref_data+M->ncells));
    //printf("Max ref: %d\n", *std::max_element(M->ref_data, M->ref_data+M->ncells));
  }

  E_Int cells_to_agglo = 0;
  for (E_Int i = 0; i < M->ncells; i++) {
    if (M->ref_data[i] < 0) cells_to_agglo += 1;
  }

  assert(cells_to_agglo % 8 == 0);
  

  // Note(Imad): should we?
  // Clip ref_data
  for (E_Int i = 0; i < M->ncells; i++) {
    if (M->ref_data[i] > 0) M->ref_data[i] = 1;
    else if (M->ref_data[i] < 0) M->ref_data[i] = -1;
  }

  puts("    Done.");
  
  // We still okay?
  /*for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    E_Int own = M->owner[i];

    E_Int oval = M->cellTree->level(own) + M->ref_data[own];
    E_Int nval = M->cellTree->level(nei) + M->ref_data[nei];

    assert(abs(oval-nval) <= 1);
  }*/
}

static
void make_ref_data_tetra(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

static
void make_ref_data_penta(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

static
void make_ref_data_pyra(E_Int cell, AMesh *M, E_Float *pM, const pDirs &Dirs)
{}

static
void compute_ref_data(AMesh *M, E_Float *metric, const std::vector<pDirs> &Dirs)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int type = M->cellTree->type(i);
    E_Float *pM = &metric[6*i];
    switch (type) {
      case HEXA:
        make_ref_data_hexa(i, M, pM, Dirs[i]);
        break;
      case TETRA:
        make_ref_data_tetra(i, M, pM, Dirs[i]);
        break;
      case PENTA:
        make_ref_data_penta(i, M, pM, Dirs[i]);
        break;
      case PYRA:
        make_ref_data_pyra(i, M, pM, Dirs[i]);
        break;
      default:
        assert(0);
        break;
    }
  }
}

static
void make_pdirs_tetra(E_Int cell, AMesh *M, pDirs &Dirs)
{}

static
void make_pdirs_penta(E_Int cell, AMesh *M, pDirs &Dirs)
{}

static
void make_pdirs_pyra(E_Int cell, AMesh *M, pDirs &Dirs)
{}

static
void reconstruct_parent_quad(E_Int face, AMesh *M, E_Int pn[4])
{
  Children *children = M->faceTree->children(face);
  if (children == NULL) {
    memcpy(pn, &M->ngon[M->indPG[face]], 4*sizeof(E_Int));
  } else {
    pn[0] = get_facets(children->pc[0], M->ngon, M->indPG)[0];
    pn[1] = get_facets(children->pc[1], M->ngon, M->indPG)[1];
    pn[2] = get_facets(children->pc[2], M->ngon, M->indPG)[2];
    pn[3] = get_facets(children->pc[3], M->ngon, M->indPG)[3];
  }
}

static
void make_pdirs_hexa(E_Int cell, AMesh *M, pDirs &Dirs)
{
  E_Int *pf = &M->nface[M->indPH[cell]];

  E_Int p0[4], p1[4];
  E_Float f0[3], f1[3];

  // LFT && RGT
  reconstruct_parent_quad(pf[2], M, p0);
  reconstruct_parent_quad(pf[3], M, p1);
  f0[0] = f0[1] = f0[2] = 0.0;
  f1[0] = f1[1] = f1[2] = 0.0;

  for (E_Int i = 0; i < 4; i++) {
    f0[0] += M->x[p0[i]];
    f0[1] += M->y[p0[i]];
    f0[2] += M->z[p0[i]];
    f1[0] += M->x[p1[i]];
    f1[1] += M->y[p1[i]];
    f1[2] += M->z[p1[i]];
  }

  for (E_Int i = 0; i < 3; i++) Dirs.I[i] = 0.25*(f1[i] - f0[i]);

  // FRO && BCK
  reconstruct_parent_quad(pf[4], M, p0);
  reconstruct_parent_quad(pf[5], M, p1);
  f0[0] = f0[1] = f0[2] = 0.0;
  f1[0] = f1[1] = f1[2] = 0.0;

  for (E_Int i = 0; i < 4; i++) {
    f0[0] += M->x[p0[i]];
    f0[1] += M->y[p0[i]];
    f0[2] += M->z[p0[i]];
    f1[0] += M->x[p1[i]];
    f1[1] += M->y[p1[i]];
    f1[2] += M->z[p1[i]];
  }

  for (E_Int i = 0; i < 3; i++) Dirs.J[i] = 0.25*(f1[i] - f0[i]);

  // LFT && RGT
  reconstruct_parent_quad(pf[0], M, p0);
  reconstruct_parent_quad(pf[1], M, p1);
  f0[0] = f0[1] = f0[2] = 0.0;
  f1[0] = f1[1] = f1[2] = 0.0;

  for (E_Int i = 0; i < 4; i++) {
    f0[0] += M->x[p0[i]];
    f0[1] += M->y[p0[i]];
    f0[2] += M->z[p0[i]];
    f1[0] += M->x[p1[i]];
    f1[1] += M->y[p1[i]];
    f1[2] += M->z[p1[i]];
  }

  for (E_Int i = 0; i < 3; i++) Dirs.K[i] = 0.25*(f1[i] - f0[i]);
}


static
void compute_principal_vecs(AMesh *M, std::vector<pDirs> &Dirs)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    switch (M->cellTree->type(i)) {
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

PyObject *K_XCORE::_metricToRefData(PyObject *self, PyObject *args)
{
  PyObject *ARRAY, *METRIC, *AMESH;

  if (!PYPARSETUPLE_(args, OOO_, &ARRAY, &METRIC, &AMESH)) {
    RAISE("Wrong input.");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  // Parse field
  E_Float *metric = NULL;
  E_Int nfld, size;
  K_NUMPY::getFromNumpyArray(METRIC, metric, size, nfld, true);
  //assert(ret == 1 && nfld == 1 && size == M->ncells*6);

  // Init ref data
  M->ref_data = (int *)XRESIZE(M->ref_data, M->ncells*sizeof(int));
  for (E_Int i = 0; i < M->ncells; i++) M->ref_data[i] = 0;

  // Compute cells principal directions
  std::vector<pDirs> Dirs(M->ncells);
  compute_principal_vecs(M, Dirs);

  // Compute ref data
  compute_ref_data(M, metric, Dirs);

  // Smooth ref data
  smooth_ref_data(M);

  return Py_None;
}

#include "proto.h"
#include <stack>

#define TAG 1.0
#define TOL 1e-12

static
E_Int FEQ(E_Float a, E_Float b)
{
  return fabs(a-b) < TOL;
}

static
void get_higher_level_neighbours(E_Int cell, E_Int face, mesh *M, std::vector<E_Int> &neis)
{
  auto& FT = M->ftree;
  E_Int nchild = tree_get_nchildren(&FT, face);
  E_Int *children = tree_get_children(&FT, face);
  for (E_Int i = 0; i < nchild; i++) {
    E_Int fchild = children[i];
    E_Int nei = M->owner[face] == cell ? M->neigh[fchild] : M->owner[fchild];
    neis.push_back(nei);
  }
}

PyObject *K_XCORE::adaptMeshDir(PyObject *self, PyObject *args)
{
  PyObject *MESH, *arr, *FIELD;
  if (!PYPARSETUPLE_(args, OOO_, &MESH, &arr, &FIELD)) {
    PyErr_SetString(PyExc_ValueError, "adaptMeshDir: wrong input.");
    return NULL;
  }

  E_Int ret;

  // Unpack adaptMesh
  if (!PyCapsule_IsValid(MESH, "adaptMesh")) {
    PyErr_SetString(PyExc_TypeError, "adaptMeshDir: bad adaptMesh hook.");
    return NULL;
  }
  mesh *M = (mesh *)PyCapsule_GetPointer(MESH, "adaptMesh");

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  ret = K_ARRAY::getFromArray3(arr, varString, f, ni, nj, nk, cn, eltType);
  
  if (ret <= 0) {
    PyErr_SetString(PyExc_TypeError, "adaptMeshDir(): only for NGons.");
    return NULL;
  }

  if (ret == 1) { 
    PyErr_SetString(PyExc_TypeError, "adaptMeshDir(): only for NGons."); 
    RELEASESHAREDS(arr, f);
    return NULL; 
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(ret, arr, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "adaptMeshDir(): can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  E_Float *x = f->begin(posx);
  E_Float *y = f->begin(posy);
  E_Float *z = f->begin(posz);

  // Check field
  E_Float *fld;
  E_Int field_size = -1, field_nfld = -1;
  ret = K_NUMPY::getFromNumpyArray(FIELD, fld, field_size, field_nfld, true);
  
  if (ret != 1) {
    PyErr_SetString(PyExc_TypeError, "adaptMeshDir: field should be an array.");
    return NULL;
  }

  if (field_size != cn->getNElts() || field_nfld != 1) {
    assert(0);
    PyErr_SetString(PyExc_TypeError, "adaptMeshDir: bad field size.");
    return NULL;
  }

  // Compute hessians on array
  E_Int nefaces = cn->getNFaces();
  E_Int *owner = (E_Int *)XMALLOC(nefaces * sizeof(E_Int));
  E_Int *neigh = (E_Int *)XMALLOC(nefaces * sizeof(E_Int));
  K_CONNECT::build_parent_elements_ngon(*cn, owner, neigh);

  E_Float *fcenters = (E_Float *)XMALLOC(3*nefaces * sizeof(E_Float));
  E_Float *fareas = (E_Float *)XMALLOC(3*nefaces * sizeof(E_Float));
  compute_face_centers_and_areas(*cn, x, y, z, fcenters, fareas);

  E_Int necells = cn->getNElts();
  E_Float *ccenters = (E_Float *)XMALLOC(3*necells * sizeof(E_Float));
  compute_cell_centers_and_vols(*cn, x, y, z, owner, neigh, fcenters,
    fareas, ccenters, NULL);

  std::vector<E_Float *> Gs(1), flds;
  for (auto& G : Gs)
    G = (E_Float *)XMALLOC(3*necells * sizeof(E_Float));
  flds.push_back(fld);
  compute_gradients_ngon(*cn, x, y, z, owner, neigh, ccenters, flds, Gs);
 
  std::vector<E_Float *> Hs(1);
  for (auto &H : Hs)
    H = (E_Float *)XMALLOC(6*necells * sizeof(E_Float));
  compute_hessians_ngon(*cn, x, y, z, owner, neigh, ccenters, Gs, flds,
    Hs);

  // Transform hessians to metrics
  for (auto& H : Hs) hessian_to_metric(H, necells);

  // Compute principal directions (HEXA)
  const auto &CT = M->ctree;
  std::vector<pDirs> Dirs(necells);
  for (E_Int i = 0; i < necells; i++) {
    E_Int cell = CT.l2g[i];
    assert(cell >= 0 && cell < M->ncells);
    E_Int *pf = &M->NFACE[6*cell];
    E_Float *p0, *p1;
    auto &dirI = Dirs[i].I;
    auto &dirJ = Dirs[i].J;
    auto &dirK = Dirs[i].K;

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

  // Compute raw ref data
  M->ref_data = (E_Int *)XRESIZE(M->ref_data, 3*M->ncells*sizeof(E_Int));

  // TODO(Imad): hardcoded params
  M->Tr = 0.5;
  M->Gmax = 2;

  for (const auto &H : Hs) {
    // TODO(Imad): take max metric
    for (E_Int i = 0; i < necells; i++) {
      E_Float *pH = &H[6*i];
      E_Int cell = CT.l2g[i];
      E_Int *pr = &M->ref_data[3*cell];
      E_Float L0, L1, L2, dd[3];
      symmat_dot_vec(pH, Dirs[i].I, dd);
      L0 = norm(dd, 3);
      symmat_dot_vec(pH, Dirs[i].J, dd);
      L1 = norm(dd, 3);
      symmat_dot_vec(pH, Dirs[i].K, dd);
      L2 = norm(dd, 3);

      // Get max stretch
      E_Float MAX = std::max(std::max(L0, L1), L2);
      if (MAX < M->Tr) continue;
      E_Float coeff = M->Gmax/MAX;
      if (L0 >= M->Tr) pr[0] = round(coeff*L0);
      if (L1 >= M->Tr) pr[1] = round(coeff*L1);
      if (L2 >= M->Tr) pr[2] = round(coeff*L2);
    }
  }

  // Smooth ref data
  smooth_ref_data_seq(M);


  for (auto& H : Hs) XFREE(H);
  for (auto& G : Gs) XFREE(G);

  Py_INCREF(Py_None);
  return Py_None;
}

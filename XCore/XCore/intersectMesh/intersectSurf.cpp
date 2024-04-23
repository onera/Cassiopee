#include "proto.h"

PyObject *K_XCORE::intersectSurf(PyObject *self, PyObject *args)
{
  PyObject *MASTER, *SLAVE, *PATCH, *TAG;
  
  if (!PYPARSETUPLE_(args, OOOO_, &MASTER, &SLAVE, &PATCH, &TAG)) {
    RAISE("Bad input.");
    return NULL;
  }

  // Check master mesh
  E_Int ret;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *fm;
  K_FLD::FldArrayI *cnm;
  char *varString, *eltType;
  ret = K_ARRAY::getFromArray3(MASTER, varString, fm, ni, nj, nk, cnm, eltType);

  if (ret <= 0) {
    RAISE("Bad master mesh");
    return NULL;
  }

  if (ret == 1) {
    RAISE("Master mesh is not NGon.");
    RELEASESHAREDS(MASTER, fm);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RAISE("Bad master point coordinates.");
    return NULL;
  }

  posx++; posy++; posz++;

  E_Float *Xm = fm->begin(posx);
  E_Float *Ym = fm->begin(posy);
  E_Float *Zm = fm->begin(posz);
  E_Int npm = fm->getSize();

  // Check slave mesh
  K_FLD::FldArrayF *fs;
  K_FLD::FldArrayI *cns;
  ret = K_ARRAY::getFromArray3(SLAVE, varString, fs, ni, nj, nk, cns, eltType);

  if (ret <= 0) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RAISE("Bad slave mesh");
    return NULL;
  }

  if (ret != 2) {
    RAISE("Slave mesh is not NGon.");
    RELEASESHAREDU(MASTER, fm, cnm);
    RELEASESHAREDB(ret, SLAVE, fs, cns);
    return NULL;
  }

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RELEASESHAREDU(SLAVE, fs, cns);
    RAISE("Bad slave point coordinates.");
    return NULL;
  }

  posx++; posy++; posz++;

  E_Float *Xs = fs->begin(posx);
  E_Float *Ys = fs->begin(posy);
  E_Float *Zs = fs->begin(posz);
  E_Int nps = fs->getSize();

  // Check intersection patch
  E_Int *mpatch = NULL;
  E_Int mpatch_size = -1;
  ret = K_NUMPY::getFromNumpyArray(PATCH, mpatch, mpatch_size, true);
  if (ret != 1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RELEASESHAREDU(SLAVE, fs, cns);
    RAISE("Bad master patch.");
    return NULL;
  }

  printf("Master patch: %d faces\n", mpatch_size);

  // Check slave point tags
  E_Float *tag = NULL;
  E_Int tag_size = -1;
  ret = K_NUMPY::getFromNumpyArray(TAG, tag, tag_size, true);
  if (ret != 1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RELEASESHAREDU(SLAVE, fs, cns);
    RAISE("Bad slave points tag.");
    return NULL;
  }

  // Init and orient master/slave meshes
  Mesh *M = mesh_init(*cnm, Xm, Ym, Zm, npm);
  Mesh *S = mesh_init(*cns, Xs, Ys, Zs, nps);

  // Extract Mf and Sf, the planar surfaces to intersect
  // TODO(Imad): quasi-planar surfaces
  std::vector<E_Int> spatch;
  for (E_Int i = 0; i < S->nf; i++) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(i, np, S);
    assert(np == 4 || np == 3);

    int keep = 1;

    for (E_Int j = 0; j < np; j++) {
      E_Int point = pn[j];
      if (tag[point] == 0) {
        keep = 0;
        break;
      }
    }

    if (keep) spatch.push_back(i);
  }

  mesh_patch_intersect(M, S, mpatch, mpatch_size, &spatch[0], spatch.size());

  return Py_None;
}
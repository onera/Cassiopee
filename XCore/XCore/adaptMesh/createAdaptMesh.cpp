/*    
    Copyright 2013-2023 Onera.

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

PyObject *K_XCORE::createAdaptMesh(PyObject *self, PyObject *args)
{
  PyObject *arr;
  if (!PYPARSETUPLE_(args, O_, &arr)) {
    PyErr_SetString(PyExc_ValueError, "adaptMeshSeq: wrong input.");
    return NULL;
  }

  // Check input
  E_Int ret;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *f;
  K_FLD::FldArrayI *cn;
  char *varString, *eltType;
  ret = K_ARRAY::getFromArray3(arr, varString, f, ni, nj, nk, cn, eltType);
  
  if (ret <= 0) {
    PyErr_SetString(PyExc_TypeError, "createAdaptMesh: bad mesh.");
    return NULL;
  }

  if (ret == 1) {
    PyErr_SetString(PyExc_TypeError, "createAdaptMesh: mesh is not NGon.");
    RELEASESHAREDS(arr, f);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  
  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "createAdaptMesh: bad coordinates.");
    return NULL;
  }

  posx++; posy++; posz++;

  // Check mesh validity
  ret = K_CONNECT::check_overlapping_cells(*cn);
  if (ret != 0) {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "createAdaptMesh: cells should not overlap.");
    return NULL;
  }

  ret = K_CONNECT::check_open_cells(*cn, NULL);
  if (ret != 0) {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "createAdaptMesh: cells should be closed.");
    return NULL;
  }

  // Apply consistent orientation on mesh faces
  E_Float *px = f->begin(posx);
  E_Float *py = f->begin(posy);
  E_Float *pz = f->begin(posz);
  ret = K_CONNECT::orient_boundary_ngon(px, py, pz, *cn);
  if (ret != 0) {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError,
      "createAdaptMesh: failed boundary orientation.");
    return NULL;
  }

  // Init mesh
  mesh *M = new mesh();
  M->nfaces = cn->getNFaces();
  M->owner = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  M->neigh = (E_Int *)XMALLOC(M->nfaces * sizeof(E_Int));
  
  // Orient the internal faces and build parent elements
  K_CONNECT::build_parent_elements_ngon(*cn, M->owner, M->neigh);

  E_Int *ngon = cn->getNGon();
  E_Int *nface = cn->getNFace();
  E_Int *indPG = cn->getIndPG();
  E_Int *indPH = cn->getIndPH();
  M->ncells = cn->getNElts();
  M->npoints = f->getSize();
  
  M->xcells = (E_Int *)XMALLOC((M->ncells+1) * sizeof(E_Int));
  for (E_Int i = 0; i < M->ncells+1; i++) M->xcells[i] = indPH[i];

  M->xfaces = (E_Int *)XMALLOC((M->nfaces+1) * sizeof(E_Int));
  for (E_Int i = 0; i < M->nfaces+1; i++) M->xfaces[i] = indPG[i];

  M->NFACE = (E_Int *)XMALLOC(M->xcells[M->ncells] * sizeof(E_Int));
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int *pf = cn->getElt(i, nf, nface, indPH);
    E_Int *ptr = &M->NFACE[M->xcells[i]];
    for (E_Int j = 0; j < nf; j++)
      ptr[j] = pf[j]-1;
  }
  
  M->NGON = (E_Int *)XMALLOC(M->xfaces[M->nfaces] * sizeof(E_Int));
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = cn->getFace(i, np, ngon, indPG);
    E_Int *ptr = &M->NGON[M->xfaces[i]];
    for (E_Int j = 0; j < np; j++)
      ptr[j] = pn[j]-1;
  }
 
  M->ctree.setSizeAndStride(M->ncells, 8);
  M->ftree.setSizeAndStride(M->nfaces, 4);
  M->xyz = (E_Float *)XMALLOC(3*M->npoints * sizeof(E_Float));
  for (E_Int i = 0; i < M->npoints; i++) {
    E_Float *ptr = &M->xyz[3*i];
    ptr[0] = px[i];
    ptr[1] = py[i];
    ptr[2] = pz[i];
  }

  // Reorder hexa to follow canon config
  reorder_hexa(M);
  compute_cell_centers(M);

  // Make hook
  PyObject *hook = PyCapsule_New((void *)M, "adaptMesh", NULL);

  RELEASESHAREDU(arr, f, cn);

  return hook;
}

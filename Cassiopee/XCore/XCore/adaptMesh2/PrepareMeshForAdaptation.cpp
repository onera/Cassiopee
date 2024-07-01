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
#include "Proto.h"

PyObject *K_XCORE::_prepareMeshForAdaptation(PyObject *self, PyObject *args)
{
  PyObject *ARRAY;
  if (!PYPARSETUPLE_(args, O_, &ARRAY)) {
    RAISE("Wrong input.");
    return NULL;
  }

  // Check input
  E_Int ret;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *f;
  K_FLD::FldArrayI *cn;
  char *varString, *eltType;
  ret = K_ARRAY::getFromArray3(ARRAY, varString, f, ni, nj, nk, cn, eltType);

  if (ret <= 0) {
    RAISE("Bad mesh.");
    return NULL;
  }

  if (ret == 1) {
    RAISE("Mesh is not NGon.");
    RELEASESHAREDS(ARRAY, f);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(ARRAY, f, cn);
    RAISE("Bad coordinates.");
    return NULL;
  }

  posx++; posy++; posz++;

  E_Float *px = f->begin(posx);
  E_Float *py = f->begin(posy);
  E_Float *pz = f->begin(posz);

  E_Int nfaces = cn->getNFaces();
  npy_intp dims[2];
  dims[0] = (npy_intp)nfaces;
  dims[1] = 1;
  PyArrayObject *OWN = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  PyArrayObject *NEI = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *owner = (E_Int *)PyArray_DATA(OWN);
  E_Int *neigh = (E_Int *)PyArray_DATA(NEI);

  // Make parent elements
  ret = K_CONNECT::orient_boundary_ngon(px, py, pz, *cn);
  assert(ret == 0);

  K_CONNECT::build_parent_elements_ngon(*cn, owner, neigh);

  // Renumber faces: internal faces first
  // Boundary faces are grouped together
    
  PyObject *out = Py_BuildValue("[OO]", (PyObject *)OWN, (PyObject *)NEI);
  Py_DECREF(OWN);
  Py_DECREF(NEI);
  return out;

































}

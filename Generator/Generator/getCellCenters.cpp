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

#include "generator.h"

/* Returns cell centers of a well-oriented NGon */

PyObject *K_GENERATOR::getCellCenters(PyObject *self, PyObject *args)
{
  PyObject* arr, *PE;

  if (!PYPARSETUPLE_(args, OO_, &arr, &PE)) return NULL;

  E_Int ret;

  // Check NGon Connectivity and coordinates
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  ret = K_ARRAY::getFromArray3(arr, varString, f, ni, nj, nk, cn, eltType);

  if (ret <= 0) {
    PyErr_SetString(PyExc_TypeError, "Bad mesh.");
    return NULL;
  }

  if (ret == 1) {
    PyErr_SetString(PyExc_TypeError, "Only for NGons.");
    RELEASESHAREDS(arr, f);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(ret, arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "Coordinates not found.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float *x = f->begin(posx);
  E_Float *y = f->begin(posy);
  E_Float *z = f->begin(posz);

  // Parent Elements
  E_Int nfaces = cn->getNFaces();
  E_Int nf, nfld;
  E_Int *pe;
  ret = K_NUMPY::getFromNumpyArray(PE, pe, nf, nfld, true);
  if (ret != 1 || nf != nfaces || nfld != 2) {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "Bad PE.");
    return NULL;
  }

  E_Int *owner = pe;
  E_Int *neigh = pe + nfaces;

  // Face centers and areas
  E_Float *fcenters = (E_Float *)malloc(3*nfaces * sizeof(E_Float));
  E_Float *fareas = (E_Float *)malloc(3*nfaces * sizeof(E_Float));

  K_METRIC::compute_face_centers_and_areas(*cn, x, y, z, fcenters, fareas);

  // Cell centers
  E_Int ncells = cn->getNElts();
  PyObject *Cx = K_NUMPY::buildNumpyArray(ncells, 1, 0, 1);
  PyObject *Cy = K_NUMPY::buildNumpyArray(ncells, 1, 0, 1);
  PyObject *Cz = K_NUMPY::buildNumpyArray(ncells, 1, 0, 1);
  E_Float *cx = K_NUMPY::getNumpyPtrF(Cx);
  E_Float *cy = K_NUMPY::getNumpyPtrF(Cy);
  E_Float *cz = K_NUMPY::getNumpyPtrF(Cz);

  K_METRIC::compute_cell_centers_and_vols(*cn, x, y, z, owner, neigh, fcenters,
    fareas, cx, cy, cz, NULL);

  free(fcenters);
  free(fareas);

  PyObject *out = PyList_New(0);
  PyList_Append(out, Cx);
  PyList_Append(out, Cy);
  PyList_Append(out, Cz);
  Py_DECREF(Cx);
  Py_DECREF(Cy);
  Py_DECREF(Cz);

  return out;
}

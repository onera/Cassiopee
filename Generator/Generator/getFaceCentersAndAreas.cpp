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

/* Returns face areas and centers of a well-oriented NGon */

PyObject *K_GENERATOR::getFaceCentersAndAreas(PyObject *self, PyObject *args)
{
  PyObject* arr;

  if (!PYPARSETUPLE_(args, O_, &arr)) return NULL;

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
  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDB(ret, arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "Coordinates not found.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float *x = f->begin(posx);
  E_Float *y = f->begin(posy);
  E_Float *z = f->begin(posz);

  // Face centers and areas
  E_Int nfaces = cn->getNFaces();
  PyObject *fcenters = K_NUMPY::buildNumpyArray(3*nfaces, 1, 0, 0);
  PyObject *fareas = K_NUMPY::buildNumpyArray(3*nfaces, 1, 0, 0);
  E_Float *fc = K_NUMPY::getNumpyPtrF(fcenters);
  E_Float *fa = K_NUMPY::getNumpyPtrF(fareas);

  K_METRIC::compute_face_centers_and_areas(*cn, x, y, z, fc, fa);

  PyObject *out = PyList_New(0);
  PyList_Append(out, fcenters);
  PyList_Append(out, fareas);
  Py_DECREF(fc);
  Py_DECREF(fa);

  return out;
}

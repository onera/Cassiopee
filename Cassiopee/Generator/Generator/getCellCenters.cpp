/*    
    Copyright 2013-2025 Onera.

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

PyObject *K_GENERATOR::getCellCenters(PyObject *self, PyObject *args)
{
  PyObject *ARR, *FC, *FA, *OWN, *NEI;

  if (!PYPARSETUPLE_(args, OOOO_ O_, &ARR, &FC, &FA, &OWN, &NEI)) return NULL;

  E_Int ret;

  // Check NGon Connectivity and coordinates
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  ret = K_ARRAY::getFromArray3(ARR, varString, f, ni, nj, nk, cn, eltType);

  if (ret <= 0) {
    PyErr_SetString(PyExc_TypeError, "Bad mesh.");
    return NULL;
  }

  if (ret == 1)
  {
    PyErr_SetString(PyExc_TypeError, "Only for NGons.");
    RELEASESHAREDS(ARR, f);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(ret, ARR, f, cn);
    PyErr_SetString(PyExc_ValueError, "Coordinates not found.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float *x = f->begin(posx);
  E_Float *y = f->begin(posy);
  E_Float *z = f->begin(posz);

  // Face centers and areas
  E_Int nfaces = cn->getNFaces();
  E_Float *fc, *fa;
  E_Int size, nfld;
  ret = K_NUMPY::getFromNumpyArray(FC, fc, size, nfld, true);
  assert(ret == 1 && size == nfaces*3 && nfld == 1);
  ret = K_NUMPY::getFromNumpyArray(FA, fa, size, nfld, true);
  assert(ret == 1 && size == nfaces*3 && nfld == 1);

  // Parent elements
  E_Int *owner, *neigh;
  if (OWN != Py_None && NEI != Py_None)
  {
    ret = K_NUMPY::getFromNumpyArray(OWN, owner, size, nfld, true);
    assert(ret == 1 && size == nfaces && nfld == 1);
    ret = K_NUMPY::getFromNumpyArray(NEI, neigh, size, nfld, true);
    assert(ret == 1 && size == nfaces && nfld == 1);
  }
  else if (OWN == Py_None && NEI == Py_None)
  {
    K_CONNECT::orient_boundary_ngon(x, y, z, *cn);
    owner = (E_Int *)malloc(nfaces * sizeof(E_Int));
    assert(owner);
    neigh = (E_Int *)malloc(nfaces * sizeof(E_Int));
    assert(neigh);
    K_CONNECT::build_parent_elements_ngon(*cn, &owner[0], &neigh[0]);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "getCellCenters: wrong arguments.");
    RELEASESHAREDS(ARR, f);
    return NULL;
  }
  
  // Cell centers
  E_Int ncells = cn->getNElts();
  PyObject *Cx = K_NUMPY::buildNumpyArray(ncells, 1, 0, 1);
  PyObject *Cy = K_NUMPY::buildNumpyArray(ncells, 1, 0, 1);
  PyObject *Cz = K_NUMPY::buildNumpyArray(ncells, 1, 0, 1);
  E_Float *cx = K_NUMPY::getNumpyPtrF(Cx);
  E_Float *cy = K_NUMPY::getNumpyPtrF(Cy);
  E_Float *cz = K_NUMPY::getNumpyPtrF(Cz);

  K_METRIC::compute_cell_centers_and_vols(*cn, x, y, z, owner, neigh, fc, fa,
    cx, cy, cz, NULL);

  PyObject *out = PyList_New(0);
  PyList_Append(out, Cx);
  PyList_Append(out, Cy);
  PyList_Append(out, Cz);
  Py_DECREF(Cx);
  Py_DECREF(Cy);
  Py_DECREF(Cz);

  if (OWN == Py_None && NEI == Py_None)
  {
    free(owner);
    free(neigh);
  }
  
  RELEASESHAREDU(ARR, f, cn);

  return out;
}

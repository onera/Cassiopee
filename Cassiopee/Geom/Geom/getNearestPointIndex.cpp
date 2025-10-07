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

# include "geom.h"
using namespace K_FLD;
using namespace std;
using namespace K_SEARCH;

//=============================================================================
/* Return the nearest point index.
   IN: array: un array
   IN: liste pts sous forme [(x,y,z)]
   OUT: liste [(ind,distance2)] (indices du point proche + distance au carre)
 */
//=============================================================================
PyObject* K_GEOM::getNearestPointIndex(PyObject* self, PyObject* args)
{
  PyObject* array, *listPts;
  if ( !PYPARSETUPLE_(args, OO_, &array, &listPts) ) return NULL;

  // Recuperation de la liste des points
  E_Int npts = PyList_Size(listPts);
  for (E_Int i = 0; i < npts;i++)
  {
    PyObject* tpli = PyList_GetItem(listPts, i);
    if (PyTuple_Check(tpli) == 0)
    {

      PyErr_SetString(PyExc_TypeError,
                      "getNearestPointIndex: each element of the list must be (x,y,z).");
      return NULL;
    }
    E_Int dim = PyTuple_Size(tpli);
    // verifie que chq element de la liste est un triplet (x,y,z)
    if ( dim != 3 )
    {
      PyErr_SetString(PyExc_TypeError,
                      "getNearestPointIndex: 3 coordinates must be found in each point of the list.");
      return NULL;
    }
  }

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f,
                                     ni, nj, nk, cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, "getNearestPointIndex: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getNearestPointIndex: array must be contains coordinates.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }
  posx++; posy++; posz++;

  // Recherche par KdTree
  ArrayAccessor<FldArrayF> coordAcc(*f, posx, posy, posz);
  E_Float* coordx = f->begin(posx);
  E_Float* coordy = f->begin(posy);
  E_Float* coordz = f->begin(posz);

  KdTree<FldArrayF> kdt(coordAcc);
  E_Float pt[3];

  PyObject* tpl = PyList_New(0);

  for (E_Int i = 0; i < npts; i++)
  {
    PyObject* tpli = PyList_GetItem(listPts, i);
    PyObject* cx = PyTuple_GetItem(tpli, 0);
    PyObject* cy = PyTuple_GetItem(tpli, 1);
    PyObject* cz = PyTuple_GetItem(tpli, 2);
    E_Float xx = PyFloat_AsDouble(cx);
    E_Float yy = PyFloat_AsDouble(cy);
    E_Float zz = PyFloat_AsDouble(cz);

    pt[0] = xx; pt[1] = yy; pt[2] = zz;
    E_Int ind = kdt.getClosest(pt);
    E_Float dx = coordx[ind]-xx;
    E_Float dy = coordy[ind]-yy;
    E_Float dz = coordz[ind]-zz;
    E_Float dist2 = dx*dx+dy*dy+dz*dz;
    PyObject* out = Py_BuildValue("(l,d)", ind, dist2);
    PyList_Append(tpl, out);
  }

  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}

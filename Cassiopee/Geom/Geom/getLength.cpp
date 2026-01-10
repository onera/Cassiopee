/*
    Copyright 2013-2026 ONERA.

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

// Information on geometries

# include <string.h>
# include "geom.h"
using namespace K_FUNC;
using namespace K_FLD;
using namespace std;
using namespace K_CONST;

// ============================================================================
/* Get the length of a 1D array defining a mesh */
// ============================================================================
PyObject* K_GEOM::getLength(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  E_Int posx, posy, posz;
  E_Float length;
  E_Float dx, dy, dz;

  length = 0.;
  if (res == 1)
  {
    if (jm != 1 || km != 1)
      printf("Warning: getLength: only line j=1, k=1 is taken into account.\n");

    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "getLength: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);

    for (E_Int i = 1; i < im; i++)
    {
      E_Int i1 = i-1;
      dx = xt[i] - xt[i1];
      dy = yt[i] - yt[i1];
      dz = zt[i] - zt[i1];
      length += sqrt(dx*dx+dy*dy+dz*dz);
    }
    RELEASESHAREDS(array, f);
  }
  else if (res == 2)
  {
    if (strcmp(eltType, "BAR") == 0)
    {
      posx = K_ARRAY::isCoordinateXPresent(varString);
      posy = K_ARRAY::isCoordinateYPresent(varString);
      posz = K_ARRAY::isCoordinateZPresent(varString);
      if (posx == -1 || posy == -1 || posz == -1)
      {
        RELEASESHAREDU(array, f, cn);
        PyErr_SetString(PyExc_TypeError,
                        "getLength: coordinates not found in array.");
        return NULL;
      }
      posx++; posy++; posz++;

      E_Float* xt = f->begin(posx);
      E_Float* yt = f->begin(posy);
      E_Float* zt = f->begin(posz);

      E_Int ind1, ind2;
      for (E_Int i = 0; i < cn->getSize(); i++)
      {
        ind1 = (*cn)(i,1)-1;
        ind2 = (*cn)(i,2)-1;
        dx = xt[ind2] - xt[ind1];
        dy = yt[ind2] - yt[ind1];
        dz = zt[ind2] - zt[ind1];
        length += sqrt(dx*dx + dy*dy + dz*dz);
      }
      RELEASESHAREDU(array, f, cn);
    }
    else
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "getLength: not a valid element type.");
      return NULL;
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLength: invalid array.");
    return NULL;
  }
#ifdef E_DOUBLEREAL
  return Py_BuildValue("d", length);
#else
  return Py_BuildValue("f", length);
#endif
}

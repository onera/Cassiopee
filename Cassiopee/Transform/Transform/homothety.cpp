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

# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
// Homothety array2/3 - in place
//=============================================================================
PyObject* K_TRANSFORM::homothety(PyObject* self, PyObject* args)
{
  E_Float xc, yc, zc, alpha;
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_ TRRR_ R_,
                    &array, &xc, &yc, &zc, &alpha))
      return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl,
                                     cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "homothety: invalid array.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "homothety: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  E_Int npts = f->getSize();
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

 #pragma omp parallel default(shared)
  {
    E_Float rx, ry, rz;
#pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      rx = xt[ind]-xc;
      ry = yt[ind]-yc;
      rz = zt[ind]-zc;
      xt[ind] = xc + alpha*rx;
      yt[ind] = yc + alpha*ry;
      zt[ind] = zc + alpha*rz;
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}

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
/* Symetry */
/* make a symetric mesh considering the plane passing */
/* by (x0,y0,z0) and of director vectors (v1x,v1y,v1z) */
/* and (v2x,v2y,v2z). */
// ============================================================================
PyObject* K_TRANSFORM::symetrize(PyObject* self, PyObject* args)
{
  E_Float x0, y0, z0;
  E_Float v1x, v1y, v1z;
  E_Float v2x, v2y, v2z;
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ TRRR_,
                    &array, &x0, &y0, &z0, &v1x, &v1y, &v1z, &v2x, &v2y, &v2z))
  {
      return NULL;
  }
  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res; E_Float norm;

  res =
    K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, cn, eltType);

  if (res == 1 || res == 2)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      if (res == 2) delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "symetrize: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    E_Int npts = f->getSize();

    // Normalisation des vecteurs du plan
    norm = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
    v1x = v1x / norm;
    v1y = v1y / norm;
    v1z = v1z / norm;
    norm = sqrt(v2x*v2x + v2y*v2y + v2z*v2z);
    v2x = v2x / norm;
    v2y = v2y / norm;
    v2z = v2z / norm;

    // 3 points du plan
    E_Float p1[3]; E_Float p2[3]; E_Float p3[3];
    p1[0] = x0; p1[1] = y0; p1[2] = z0;
    p2[0] = x0+v1x; p2[1] = y0+v1y; p2[2] = z0+v1z;
    p3[0] = x0+v2x; p3[1] = y0+v2y; p3[2] = z0+v2z;

    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);

    // Symetry
#pragma omp parallel
    {
        E_Float sigma1, sigma0, dist2;
        E_Bool in;
        E_Float xint, yint, zint;
        E_Float p[3];
        #pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
            p[0] = xp[i]; p[1] = yp[i]; p[2] = zp[i];
            K_COMPGEOM::distanceToTriangle(p1, p2, p3, p, 0,
				     dist2, in, xint, yint, zint,
				     sigma0, sigma1);

            xp[i] = 2*xint - xp[i];
            yp[i] = 2*yint - yp[i];
            zp[i] = 2*zint - zp[i];
        }
    }
    RELEASESHAREDB(res, array, f, cn);
    Py_INCREF(Py_None);
    return Py_None;

  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "symetrize: not a valid array.");
    return NULL;
  }
}

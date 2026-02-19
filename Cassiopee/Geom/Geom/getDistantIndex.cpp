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
/* Get the distant index */
// ============================================================================
PyObject* K_GEOM::getDistantIndex(PyObject* self, PyObject* args)
{
  E_Float eps = 1.e-10;
  PyObject* array;
  E_Int ind;
  E_Float l;

  if (!PYPARSETUPLE_(args, O_ I_ R_,
                    &array, &ind, &l))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  E_Int imjm, kind, jind;
  E_Float t, dx, dy, dz;
  E_Int is, bsup, binf;
  E_Float length;

  if (res == 1)
  {
    if (jm != 1 || km != 1)
      printf("Warning: getDistantIndex: only line j=1, k=1 is taken into account.\n");
    // Data check
    if (ind > im*jm*km || ind < 1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "getDistantIndex: index is out of mesh bounds.");
      return NULL;
    }

    E_Int posx = K_ARRAY::isCoordinateXPresent( varString );
    E_Int posy = K_ARRAY::isCoordinateYPresent( varString );
    E_Int posz = K_ARRAY::isCoordinateZPresent( varString );
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "getDistantIndex: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    is = ind-1;
    imjm = im*jm;
    kind = is / imjm;
    jind = (is - kind*imjm)/im;

    bsup = im-1 + jind*im + kind*imjm;
    binf = jind*im + kind*imjm;

    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);

    length = 0.;
    if (l >= 0.)
    {
      while (length < l && is < bsup)
      {
        dx = xt[is+1] - xt[is];
        dy = yt[is+1] - yt[is];
        dz = zt[is+1] - zt[is];
        t = sqrt( dx*dx + dy*dy + dz*dz );
        length = length + t;
        is++;
      }
      if (is == bsup && length < l-eps)
      {
        RELEASESHAREDS(array, f);
        PyErr_SetString(PyExc_ValueError,
                        "getDistantIndex: max bound reached without matching length.");
        return NULL;
      }
    }
    else
    {
      while (length > l && is > binf)
      {
        dx = xt[is] - xt[is-1];
        dy = yt[is] - yt[is-1];
        dz = zt[is] - zt[is-1];
        t = sqrt( dx*dx + dy*dy + dz*dz);
        length = length - t;
        is--;
      }
      if (is == binf && length > l+eps)
      {
         RELEASESHAREDS(array, f);
        PyErr_SetString(PyExc_ValueError,
                        "getDistantIndex: min bound reached with matching length.");
        return NULL;
      }
    }

    RELEASESHAREDS(array, f);
    return Py_BuildValue(I_, is);
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getDistantIndex: can not be used on an unstructured array.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getDistantIndex: invalid array.");
    return NULL;
  }
}

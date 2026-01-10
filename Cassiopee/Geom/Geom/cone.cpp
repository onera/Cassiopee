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

// Analytical geometries creation

#include "geom.h"
#include <math.h>
using namespace K_FLD;
using namespace K_CONST;

// ============================================================================
/* Create a cone of center C, basis radius Rb, vertex radius Rv,
   and height H */
// ============================================================================
PyObject* K_GEOM::cone(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float xc, yc, zc;
  E_Float Rb, Rv, H;
  if (!PYPARSETUPLE_(args, TRRR_ RRR_ I_,
                    &xc, &yc, &zc, &Rb, &Rv, &H, &N))
  {
    return NULL;
  }
  E_Float pi = 4*atan(1.);

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError,
                    "cone: insufficient number of points.");
    return NULL;
  }
  if (K_FUNC::fEqualZero(H))
  {
    PyErr_SetString(PyExc_ValueError,
                    "cone: H must be non null.");
    return NULL;
  }

  // Create a cone
  E_Int api = 1; // TODO
  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", N, N, 1, api);
  FldArrayF* f;
  K_ARRAY::getFromArray3(tpl, f);
  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  E_Float alpha;
  E_Float x1, y1, z1;
  E_Float delta = 2.*pi/(N-1.);
  E_Float rapport = (Rb-Rv)/H;
  E_Float hk = H/(N-1.);
  E_Float Rk;
  z1 = 0.;

  for (E_Int k = 0; k < N; k++)
  {
    Rk = Rb - z1*rapport;
    for (E_Int i = 0; i < N; i++)
    {
      alpha = i*delta;
      x1 = Rk*cos(alpha);
      y1 = Rk*sin(alpha);
      xt[i+k*N] = xc + x1;
      yt[i+k*N] = yc + y1;
      zt[i+k*N] = zc + z1;
    }
    z1 = z1 + hk;
  }

  RELEASESHAREDS(tpl, f);
  return tpl;
}

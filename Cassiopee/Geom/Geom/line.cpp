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
/* Create a line of N points passing by P1 and P2 */
// ============================================================================
PyObject* K_GEOM::line(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float x1, y1, z1;
  E_Float x2, y2, z2;

  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ I_,
                    &x1, &y1, &z1, &x2, &y2, &z2, &N))
  {
    return NULL;
  }

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError, "line: insufficient number of point.");
    return NULL;
  }

  // Create a line
  E_Int api = 1; // TODO
  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", N, 1, 1, api);
  FldArrayF* f;
  K_ARRAY::getFromArray3(tpl, f);
  E_Float* coordx = f->begin(1);
  E_Float* coordy = f->begin(2);
  E_Float* coordz = f->begin(3);

  E_Float delta = 1./(N-1.);
  E_Float dx12 = delta * (x2 - x1);
  E_Float dy12 = delta * (y2 - y1);
  E_Float dz12 = delta * (z2 - z1);

  for (E_Int i = 0; i < N; i++)
  {
    coordx[i] = x1 + i*dx12;
    coordy[i] = y1 + i*dy12;
    coordz[i] = z1 + i*dz12;
  }

  RELEASESHAREDS(tpl, f);
  return tpl;
}

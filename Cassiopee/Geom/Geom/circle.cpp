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
/* Create a circle of center C and radius R,
   between tetas and tetae angles */
// ============================================================================
PyObject* K_GEOM::circle(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float xc, yc, zc;
  E_Float R, tetas, tetae;

  if (!PYPARSETUPLE_(args, TRRR_ RRR_ I_,
                    &xc, &yc, &zc, &R, &tetas, &tetae, &N))
  {
    return NULL;
  }
  E_Float pi = 4*atan(1.);
  E_Float t1 = tetas*pi/180.;
  E_Float t2 = tetae*pi/180.;

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError,
                    "circle: insufficient number of points.");
    return NULL;
  }

  // Create a portion of circle
  E_Int api = 1; // TODO
  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", N, 1, 1, api);
  FldArrayF* f;
  K_ARRAY::getFromArray3(tpl, f);
  E_Float* coordx = f->begin(1);
  E_Float* coordy = f->begin(2);
  E_Float* coordz = f->begin(3);

  E_Float delta = 1./(N-1.);
  E_Float alpha;
  for (E_Int i = 0; i < N; i++)
  {
    alpha = t1 + i*delta*(t2 - t1);
    coordx[i] = xc + R*cos(alpha);
    coordy[i] = yc + R*sin(alpha);
    coordz[i] = zc;
  }

  RELEASESHAREDS(tpl, f);
  return tpl;
}

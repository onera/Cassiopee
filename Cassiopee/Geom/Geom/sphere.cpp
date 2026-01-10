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
/* Create a sphere of center C and radius R */
// ============================================================================
PyObject* K_GEOM::sphere(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float xc, yc, zc;
  E_Float R;
  if (!PYPARSETUPLE_(args, TRRR_ R_ I_,
                    &xc, &yc, &zc, &R, &N))
  {
    return NULL;
  }

  E_Float pi = 4*atan(1.);

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError,
                    "sphere: insufficient number of points.");
    return NULL;
  }

  // Create a sphere
  E_Int api = 1; // TODO
  E_Int P = 2*N;
  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", N, P, 1, api);
  FldArrayF* f;
  K_ARRAY::getFromArray3(tpl, f);
  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  E_Float alpha, beta, cbeta, sbeta, x1, y1, z1;
  E_Int ind;

  E_Float delta = pi/(N-1.);
  E_Float deltap = (2*pi)/(P-1.);

  for (E_Int j = 0; j < P; j++)
    for (E_Int i = 0; i < N; i++)
    {
      alpha = i*delta;
      beta = j*deltap;
      x1 = R*cos(alpha);
      y1 = R*sin(alpha);
      z1 = 0.;
      ind = i + j*N;
      cbeta = cos(beta);
      sbeta = sin(beta);
      xt[ind] = xc + x1;
      yt[ind] = yc + cbeta*y1 - sbeta*z1;
      zt[ind] = zc + sbeta*y1 + cbeta*z1;
    }

  RELEASESHAREDS(tpl, f);
  return tpl;
}

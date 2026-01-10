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
/* Create a triangle */
// ============================================================================
PyObject* K_GEOM::triangle(PyObject* self, PyObject* args)
{
  E_Float x1, y1, z1;
  E_Float x2, y2, z2;
  E_Float x3, y3, z3;
  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ TRRR_,
                    &x1, &y1, &z1, &x2, &y2, &z2, &x3, &y3, &z3))
  {
    return NULL;
  }

  E_Int api = 1; // TODO
  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", 3, 1, "TRI", false, api);
  FldArrayF* f; FldArrayI* cn;
  K_ARRAY::getFromArray3(tpl, f, cn);
  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);

  xt[0] = x1; yt[0] = y1; zt[0] = z1;
  xt[1] = x2; yt[1] = y2; zt[1] = z2;
  xt[2] = x3; yt[2] = y3; zt[2] = z3;
  cn1[0] = 1; cn2[0] = 2; cn3[0] = 3;

  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}

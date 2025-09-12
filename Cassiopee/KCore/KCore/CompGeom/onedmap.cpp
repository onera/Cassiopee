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
# include "CompGeom/compGeom.h"

//=============================================================================
// map a 1D distribution over a profile
//=============================================================================
void K_COMPGEOM::onedmap(
  const E_Int ni,
  const E_Float* x, const E_Float* y, const E_Float* z,
  const E_Int no, const E_Float* d,
  E_Float* xo, E_Float* yo, E_Float* zo,
  E_Float* s, E_Float* dx, E_Float* dy, E_Float* dz
)
{
  E_Float stota;

  // Parametrisation de la ligne entrante
  slope(ni, x, y, z, dx, dy, dz);

  paramFunc(ni, x, y, z, dx, dy, dz, stota, s);

  // Projection
  interp(ni, no, stota, s, d, x, y, z, dx, dy, dz, xo, yo, zo);
}

//=============================================================================
void K_COMPGEOM::onedmapbar(
  const E_Int npts,
  const E_Float* x, const E_Float* y, const E_Float* z,
  const E_Int no, const E_Float* d,
  const E_Int net, const E_Int* cn1, const E_Int* cn2,
  const E_Int neto, E_Int* cn1o, E_Int* cn2o,
  E_Float* xo, E_Float* yo, E_Float* zo,
  E_Float* s, E_Float* dx, E_Float* dy, E_Float* dz
)
{
  E_Float stota;

  // Parametrisation de la ligne entrante
  slopebar(npts, net, cn1, cn2, x, y, z, dx, dy, dz);

  paramFuncBar(npts, net, cn1, cn2, x, y, z, dx, dy, dz, stota, s);

  // Projection
  interpbar(npts, no, stota, s, d, x, y, z, dx, dy, dz,
            xo, yo, zo, neto, cn1o, cn2o);
}

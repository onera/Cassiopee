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
# include "metric.h"


//=============================================================================
// Compute barycenter of tetra cells.
// IN: ind1, ind2, ind3, ind4: Sommets du tetraedre
// IN: xt, yt, zt: Vertex coordinates
// OUT: vol: Cell volume
//=============================================================================
void K_METRIC::compVolOfTetraCell(
  const E_Int ind1, const E_Int ind2, const E_Int ind3, const E_Int ind4,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& vol
)
{
  E_Float x1, x2, x3, x4;
  E_Float y1, y2, y3, y4;
  E_Float z1, z2, z3, z4;

  E_Float xint[4], yint[4], zint[4];
  E_Float snx[4], sny[4], snz[4], surf[4];
  E_Float xloc[4], yloc[4], zloc[4];
  E_Int cnloc[4];

  // Calcul des centres des facettes 123, 124, 234, 134, respectively
  x1 = xt[ind1]; y1 = yt[ind1]; z1 = zt[ind1];
  x2 = xt[ind2]; y2 = yt[ind2]; z2 = zt[ind2];
  x3 = xt[ind3]; y3 = yt[ind3]; z3 = zt[ind3];
  x4 = xt[ind4]; y4 = yt[ind4]; z4 = zt[ind4];

  // Face 123
  xint[0] = K_CONST::ONE_THIRD * (x1 + x2 + x3);
  yint[0] = K_CONST::ONE_THIRD * (y1 + y2 + y3);
  zint[0] = K_CONST::ONE_THIRD * (z1 + z2 + z3);

  // Face 124
  xint[1] = K_CONST::ONE_THIRD * (x1 + x2 + x4);
  yint[1] = K_CONST::ONE_THIRD * (y1 + y2 + y4);
  zint[1] = K_CONST::ONE_THIRD * (z1 + z2 + z4);

  // Face 234
  xint[2] = K_CONST::ONE_THIRD * (x2 + x3 + x4);
  yint[2] = K_CONST::ONE_THIRD * (y2 + y3 + y4);
  zint[2] = K_CONST::ONE_THIRD * (z2 + z3 + z4);

  // Face 134
  xint[3] = K_CONST::ONE_THIRD * (x1 + x3 + x4);
  yint[3] = K_CONST::ONE_THIRD * (y1 + y3 + y4);
  zint[3] = K_CONST::ONE_THIRD * (z1 + z3 + z4);

  xloc[0] = x1; xloc[1] = x2; xloc[2] = x3; xloc[3] = x4;
  yloc[0] = y1; yloc[1] = y2; yloc[2] = y3; yloc[3] = y4;
  zloc[0] = z1; zloc[1] = z2; zloc[2] = z3; zloc[3] = z4;

  cnloc[0] = 1; cnloc[1] = 2; cnloc[2] = 3; cnloc[3] = 4;

  // compTetraSurf(cnloc, xloc, yloc, zloc, TODO
  //               snx, sny, snz, surf);

  vol = K_CONST::E_ZERO_FLOAT;
  for (E_Int edge = 0; edge < 4; edge++)
  {
    vol += xint[edge] * snx[edge]
         + yint[edge] * sny[edge]
         + zint[edge] * snz[edge];
  }
  vol *= K_CONST::ONE_THIRD;
}

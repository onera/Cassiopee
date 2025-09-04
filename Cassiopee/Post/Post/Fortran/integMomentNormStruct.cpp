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
# include "post.h"

// ============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
//       and F have the same size
// ============================================================================
void K_POST::integMomentNormStruct(
  const E_Int ni, const E_Int nj, const E_Float cx,
  const E_Float cy, const E_Float cz, const E_Float* ratio, const E_Float* xt,
  const E_Float* yt, const E_Float* zt, const E_Float* sx, const E_Float* sy,
  const E_Float* sz, const E_Float* field, E_Float* result)
{
  E_Int ni1;
  E_Int ind, ind1, ind2, ind3, ind4;
  E_Float f, f1, f2, f3, f4;
  E_Float mx, my, mz;
  E_Float centerx, centery, centerz;
  E_Float res1, res2, res3;

  res1 = 0.0;
  res2 = 0.0;
  res3 = 0.0;

  ni1 = ni - 1;
  for (E_Int j = 0; j < nj - 1; j++)
  {
    for (E_Int i = 0; i < ni - 1; i++)
    {
      ind1 = i + j * ni;
      ind2 = ind1 + ni;
      ind3 = ind1 + 1;
      ind4 = ind3 + ni;
      ind = i + j * ni1;

      f1 = ratio[ind1] * field[ind1];
      f2 = ratio[ind2] * field[ind2];
      f3 = ratio[ind3] * field[ind3];
      f4 = ratio[ind4] * field[ind4];

      f = K_CONST::ONE_FOURTH * (f1 + f2 + f3 + f4);

      centerx = K_CONST::ONE_FOURTH * (xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4]);
      centerx = centerx - cx;

      centery = K_CONST::ONE_FOURTH * (yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4]);
      centery = centery - cy;

      centerz = K_CONST::ONE_FOURTH * (zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4]);
      centerz = centerz - cz;

      mx = centery * sz[ind] - centerz * sy[ind];
      my = centerz * sx[ind] - centerx * sz[ind];
      mz = centerx * sy[ind] - centery * sx[ind];

      res1 += f * mx;
      res2 += f * my;
      res3 += f * mz;
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

//=============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
//       are defined in nodes and F is defined in center
//=============================================================================
void K_POST::integMomentNormStructNodeCenter(
  const E_Int ni, const E_Int nj,
  const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* sx,
  const E_Float* sy, const E_Float* sz, const E_Float* field, E_Float* result)
{
  E_Int ni1;
  E_Int ind, ind1, ind2, ind3, ind4;
  E_Float f, sx0, sy0, sz0;
  E_Float mx, my, mz, resx, resy, resz;
  E_Float centerx, centery, centerz;

  resx = 0.0;
  resy = 0.0;
  resz = 0.0;

  ni1 = ni - 1;
  for (E_Int j = 0; j < nj - 1; j++)
  {
    for (E_Int i = 0; i < ni - 1; i++)
    {
      ind1 = i + j * ni;
      ind2 = ind1 + ni;
      ind3 = ind1 + 1;
      ind4 = ind3 + ni;
      ind = i + j * ni1;

      f = ratio[ind] * field[ind];

      sx0 = sx[ind];
      sy0 = sy[ind];
      sz0 = sz[ind];

      centerx = xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4];
      centerx = K_CONST::ONE_FOURTH * centerx - cx;

      centery = yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4];
      centery = K_CONST::ONE_FOURTH * centery - cy;

      centerz = zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4];
      centerz = K_CONST::ONE_FOURTH * centerz - cz;

      mx = centery * sz0 - centerz * sy0;
      my = centerz * sx0 - centerx * sz0;
      mz = centerx * sy0 - centery * sx0;

      resx += f * mx;
      resy += f * my;
      resz += f * mz;
    }
  }

  result[0] = resx;
  result[1] = resy;
  result[2] = resz;
}
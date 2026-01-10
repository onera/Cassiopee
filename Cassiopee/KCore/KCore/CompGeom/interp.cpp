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

// Cubic Hermite interpolation
# include <math.h>
# include <stdio.h>
# include "CompGeom/compGeom.h"

//=============================================================================
// Hermite cubic interpolation helper
//=============================================================================
inline E_Float K_COMPGEOM::valcub(
  E_Float a, E_Float b, E_Float c, E_Float d, E_Float t
)
{
  E_Float onemt = 1.0 - t;
  E_Float t2 = t * t;
  E_Float onemt2 = onemt * onemt;
  return (a * (1.0 + 2.0 * t) + c * t) * onemt2 +
         (b * (3.0 - 2.0 * t) + d * (t - 1.0)) * t2;
}

//=============================================================================
// Cubic Hermite interpolation
// By F.MONTIGNY-RANNOU 1994
//=============================================================================
void K_COMPGEOM::interp(
  const E_Int im0, const E_Int im, const E_Float stota,
  const E_Float* s0, const E_Float* s,
  const E_Float* tabx0, const E_Float* taby0, const E_Float* tabz0,
  const E_Float* dx0, const E_Float* dy0, const E_Float* dz0,
  E_Float* tabx, E_Float* taby, E_Float* tabz
)
{
  tabx[0] = tabx0[0];
  taby[0] = taby0[0];
  tabz[0] = tabz0[0];
  tabx[im-1] = tabx0[im0-1];
  taby[im-1] = taby0[im0-1];
  tabz[im-1] = tabz0[im0-1];

  if (im == 2) return;
  else
  {
    E_Int k;
    E_Float t, coef, d1, d2;

    for (E_Int i = 1; i < im-1; i++)
    {
      for (k = 0; k < im0-1; k++)
      {
        if (s[i] <= s0[k+1]) break;
      }
      if (k >= im0-1) k = im0-2;
      t = (s[i] - s0[k]) / (s0[k+1] - s0[k]);
      coef = stota * (s0[k+1] - s0[k]);
      d1 = dx0[k] * coef;
      d2 = dx0[k+1] * coef;
      tabx[i] = valcub(tabx0[k], tabx0[k+1], d1, d2, t);
      d1 = dy0[k] * coef;
      d2 = dy0[k+1] * coef;
      taby[i] = valcub(taby0[k], taby0[k+1], d1, d2, t);
      d1 = dz0[k] * coef;
      d2 = dz0[k+1] * coef;
      tabz[i] = valcub(tabz0[k], tabz0[k+1], d1, d2, t);
    }
  }
}

//=============================================================================
// BAR version: here the BAR is isomorphic to an i-array.
//=============================================================================
void K_COMPGEOM::interpbar(
  const E_Int im0, const E_Int im, const E_Float stota,
  const E_Float* s0, const E_Float* s,
  const E_Float* tabx0, const E_Float* taby0, const E_Float* tabz0,
  const E_Float* dx0, const E_Float* dy0, const E_Float* dz0,
  E_Float* tabx, E_Float* taby, E_Float* tabz, 
  const E_Int net, E_Int* cn1, E_Int* cn2
)
{
  tabx[0] = tabx0[0];
  taby[0] = taby0[0];
  tabz[0] = tabz0[0];
  tabx[im-1] = tabx0[im0-1];
  taby[im-1] = taby0[im0-1];
  tabz[im-1] = tabz0[im0-1];

  if (im == 2) return;
  else
  {
    E_Int k;
    E_Float t, coef, d1, d2;

    for (E_Int i = 1; i < im-1; i++)
    {
      for (k = 0; k < im0-1; k++)
      {
        if (s[i] <= s0[k+1]) break;
      }

      if (k >= im0-1) k = im0-2;
      t = (s[i] - s0[k]) / (s0[k+1] - s0[k]);
      coef = stota * (s0[k+1] - s0[k]);
      d1 = dx0[k] * coef;
      d2 = dx0[k+1] * coef;
      tabx[i] = valcub(tabx0[k], tabx0[k+1], d1, d2, t);
      d1 = dy0[k] * coef;
      d2 = dy0[k+1] * coef;
      taby[i] = valcub(taby0[k], taby0[k+1], d1, d2, t);
      d1 = dz0[k] * coef;
      d2 = dz0[k+1] * coef;
      tabz[i] = valcub(tabz0[k], tabz0[k+1], d1, d2, t);
    }
  }

  // Build BAR connectivity
  E_Int i = 0;
  for (E_Int et = 0; et < net; et++)
  {
    cn1[et] = i + 1;
    cn2[et] = i + 2;
    i += 1;
  }
}

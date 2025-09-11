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

// Cubic Hermite interpolation

# include <math.h>
# include <stdio.h>
# include "CompGeom/compGeom.h"

//****************************************************************************
// Cubic Hermite interpolation
// By F.MONTIGNY-RANNOU 1994
//****************************************************************************
void K_COMPGEOM::interp(E_Int im0, E_Int im, E_Float stota,
			const E_Float* tabx0, const E_Float* taby0, const E_Float* tabz0, const E_Float* s0,
			const E_Float* dx0, const E_Float* dy0, const E_Float* dz0,
			E_Float* tabx, E_Float* taby, E_Float* tabz, const E_Float* s)
{
  // Hermite cubic interpolation helper
  #define VALCUB(a, b, c, d, t) ((a*(1.0+2.0*t)+c*t)*pow(1.0-t,2) + (b*(3.0-2.0*t)+d*(t-1.0))*pow(t,2))

  tabx[0] = tabx0[0];
  taby[0] = taby0[0];
  tabz[0] = tabz0[0];
  tabx[im-1] = tabx0[im0-1];
  taby[im-1] = taby0[im0-1];
  tabz[im-1] = tabz0[im0-1];

  if (im == 2)
  {
    return;
  }
  else
  {
    for (E_Int i = 1; i < im-1; i++)
    {
      E_Int k;
      for (k = 0; k < im0-1; k++)
      {
        if (s[i] <= s0[k+1]) break;
      }
      if (k >= im0-1) k = im0-2;
      E_Float t = (s[i] - s0[k]) / (s0[k+1] - s0[k]);
      E_Float coef = stota * (s0[k+1] - s0[k]);
      E_Float d1 = dx0[k] * coef;
      E_Float d2 = dx0[k+1] * coef;
      tabx[i] = VALCUB(tabx0[k], tabx0[k+1], d1, d2, t);
      d1 = dy0[k] * coef;
      d2 = dy0[k+1] * coef;
      taby[i] = VALCUB(taby0[k], taby0[k+1], d1, d2, t);
      d1 = dz0[k] * coef;
      d2 = dz0[k+1] * coef;
      tabz[i] = VALCUB(tabz0[k], tabz0[k+1], d1, d2, t);
    }
  }
  #undef VALCUB
}

//=============================================================================
// BAR version: here the BAR is isomorphic to an i-array
void K_COMPGEOM::interpbar(E_Int im0, E_Int im, E_Float stota,
			   E_Int net0, const E_Int* cn10, const E_Int* cn20,
			   const E_Float* tabx0, const E_Float* taby0, const E_Float* tabz0, const E_Float* s0,
			   const E_Float* dx0, const E_Float* dy0, const E_Float* dz0,
			   E_Int net, E_Int* cn1, E_Int* cn2,
			   E_Float* tabx, E_Float* taby, E_Float* tabz, const E_Float* s)
{
  #define VALCUB(a, b, c, d, t) ((a*(1.0+2.0*t)+c*t)*pow(1.0-t,2) + (b*(3.0-2.0*t)+d*(t-1.0))*pow(t,2))

  tabx[0] = tabx0[0];
  taby[0] = taby0[0];
  tabz[0] = tabz0[0];
  tabx[im-1] = tabx0[im0-1];
  taby[im-1] = taby0[im0-1];
  tabz[im-1] = tabz0[im0-1];

  if (im == 2)
  {
    return;
  }
  else
  {
    for (E_Int i = 1; i < im-1; i++)
    {
      E_Int k;
      for (k = 0; k < im0-1; k++)
      {
        if (s[i] <= s0[k+1]) break;
      }
      if (k >= im0-1) k = im0-2;
      E_Float t = (s[i] - s0[k]) / (s0[k+1] - s0[k]);
      E_Float coef = stota * (s0[k+1] - s0[k]);
      E_Float d1 = dx0[k] * coef;
      E_Float d2 = dx0[k+1] * coef;
      tabx[i] = VALCUB(tabx0[k], tabx0[k+1], d1, d2, t);
      d1 = dy0[k] * coef;
      d2 = dy0[k+1] * coef;
      taby[i] = VALCUB(taby0[k], taby0[k+1], d1, d2, t);
      d1 = dz0[k] * coef;
      d2 = dz0[k+1] * coef;
      tabz[i] = VALCUB(tabz0[k], tabz0[k+1], d1, d2, t);
    }
  }
  // Build BAR connectivity
  E_Int i = 0;
  for (E_Int et = 0; et < net; et++)
  {
    cn1[et] = i + 1;
    cn2[et] = i + 2;
    i = i + 1;
  }
  #undef VALCUB
}

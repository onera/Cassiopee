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
// Compute surface integral of the moment M (OM^F), coordinates 
//     and field have the same size
//     I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//     Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
// ============================================================================
void K_POST::integMomentStruct(const E_Int ni, const E_Int nj,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result)
{
  E_Int ni1;
  E_Int ind1, ind2, ind3, ind4, ind;
  E_Float f1x, f2x, f3x, f4x;
  E_Float f1y, f2y, f3y, f4y;
  E_Float f1z, f2z, f3z, f4z;
  E_Float m1x, m2x, m3x, m4x;
  E_Float m1y, m2y, m3y, m4y;
  E_Float m1z, m2z, m3z, m4z;
  E_Float dx, dy, dz;
  E_Float rind, sind;
  E_Float res1, res2, res3;

  res1 = 0.0;
  res2 = 0.0;
  res3 = 0.0;

  ni1 = ni - 1;

  for (E_Int j = 0; j < nj - 1; j++)
  {
    for (E_Int i = 0; i < ni - 1; i++)
    {
      ind1 = i + j*ni;
      ind2 = ind1 + ni;
      ind3 = ind1 + 1;
      ind4 = ind3 + ni;
      ind = i + j*ni1;

      rind = ratio[ind1];
      dx = xt[ind1] - cx;
      dy = yt[ind1] - cy;
      dz = zt[ind1] - cz;

      f1x = rind * vx[ind1];
      f1y = rind * vy[ind1];
      f1z = rind * vz[ind1];
      m1x = dy * f1z - dz * f1y;
      m1y = dz * f1x - dx * f1z;
      m1z = dx * f1y - dy * f1x;

      rind = ratio[ind2];
      dx = xt[ind2] - cx;
      dy = yt[ind2] - cy;
      dz = zt[ind2] - cz;

      f2x = rind * vx[ind2];
      f2y = rind * vy[ind2];
      f2z = rind * vz[ind2];
      m2x = dy * f2z - dz * f2y;
      m2y = dz * f2x - dx * f2z;
      m2z = dx * f2y - dy * f2x;

      rind = ratio[ind3];
      dx = xt[ind3] - cx;
      dy = yt[ind3] - cy;
      dz = zt[ind3] - cz;

      f3x = rind * vx[ind3];
      f3y = rind * vy[ind3];
      f3z = rind * vz[ind3];
      m3x = dy * f3z - dz * f3y;
      m3y = dz * f3x - dx * f3z;
      m3z = dx * f3y - dy * f3x;

      rind = ratio[ind4];
      dx = xt[ind4] - cx;
      dy = yt[ind4] - cy;
      dz = zt[ind4] - cz;

      f4x = rind * vx[ind4];
      f4y = rind * vy[ind4];
      f4z = rind * vz[ind4];
      m4x = dy * f4z - dz * f4y;
      m4y = dz * f4x - dx * f4z;
      m4z = dx * f4y - dy * f4x;

      sind = surf[ind];
      res1 += sind * (m1x + m2x + m3x + m4x);
      res2 += sind * (m1y + m2y + m3y + m4y);
      res3 += sind * (m1z + m2z + m3z + m4z);
    }
  }

  result[0] = K_CONST::ONE_FOURTH * res1;
  result[1] = K_CONST::ONE_FOURTH * res2;
  result[2] = K_CONST::ONE_FOURTH * res3;
}

//=============================================================================
// Compute linear integral of the moment M (OM^F), coordinates 
//     and field have the same size
//     I(AB) = LENGTH(ABCD)*(F(A)+F(B))/2
// ============================================================================
void K_POST::integMomentStruct1D(
  const E_Int ni, const E_Float cx, const E_Float cy,
  const E_Float cz, const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* length, const E_Float* vx, const E_Float* vy,
  const E_Float* vz, E_Float* result
)
{
  E_Int ind1, ind2, ind;
  E_Float f1x, f2x;
  E_Float f1y, f2y;
  E_Float f1z, f2z;
  E_Float m1x, m2x;
  E_Float m1y, m2y;
  E_Float m1z, m2z;
  E_Float ri, li, dx, dy, dz;
  E_Float res1, res2, res3;

  res1 = 0.0;
  res2 = 0.0;
  res3 = 0.0;

  for (E_Int i = 0; i < ni - 1; i++)
  {
    ind1 = i;
    ind2 = i + 1;

    ind = i;
    li = length[ind];

    ri = ratio[ind1];
    f1x = ri * vx[ind1];
    f1y = ri * vy[ind1];
    f1z = ri * vz[ind1];

    dx = xt[ind1] - cx;
    dy = yt[ind1] - cy;
    dz = zt[ind1] - cz;

    m1x = dy * f1z - dz * f1y;
    m1y = dz * f1x - dx * f1z;
    m1z = dx * f1y - dy * f1x;

    ri = ratio[ind2];
    f2x = ri * vx[ind2];
    f2y = ri * vy[ind2];
    f2z = ri * vz[ind2];

    dx = xt[ind2] - cx;
    dy = yt[ind2] - cy;
    dz = zt[ind2] - cz;

    m2x = dy * f2z - dz * f2y;
    m2y = dz * f2x - dx * f2z;
    m2z = dx * f2y - dy * f2x;

    res1 += li * (m1x + m2x);
    res2 += li * (m1y + m2y);
    res3 += li * (m1z + m2z);
  }

  result[0] = K_CONST::ONE_HALF * res1;
  result[1] = K_CONST::ONE_HALF * res2;
  result[2] = K_CONST::ONE_HALF * res3;
}


//=============================================================================
// Compute surface integral of the moment M (OM^F), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
void K_POST::integMomentStructNodeCenter(const E_Int ni, const E_Int nj,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* surf, const E_Float* vx,
  const E_Float* vy, const E_Float* vz, E_Float* result)
{
  E_Int ni1;
  E_Int ind, ind1, ind2, ind3, ind4;
  E_Float mx, my, mz;
  E_Float res1, res2, res3;
  E_Float centerx, centery, centerz;
  E_Float srind;

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
      ind4 = ind2 + 1;
      ind = i + j*ni1;
      srind = surf[ind]*ratio[ind];
      
      centerx = xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4];
      centerx = K_CONST::ONE_FOURTH*centerx;
      centerx = centerx - cx;

      centery = yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4];
      centery = K_CONST::ONE_FOURTH*centery;
      centery = centery - cy;

      centerz = zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4];
      centerz = K_CONST::ONE_FOURTH*centerz;
      centerz = centerz - cz;

      mx = centery * vz[ind] - centerz * vy[ind];
      my = centerz * vx[ind] - centerx * vz[ind];
      mz = centerx * vy[ind] - centery * vx[ind];

      res1 += srind * mx;
      res2 += srind * my;
      res3 += srind * mz;
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

//=============================================================================
// Compute linear integral of the moment M (OM^F), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
void K_POST::integMomentStructNodeCenter1D(
  const E_Int ni, const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* length, const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result)
{
  E_Int ind, ind1, ind2;
  E_Float resx, resy, resz;
  E_Float mx, my, mz, rlind;
  E_Float centerx, centery, centerz;

  resx = 0.0;
  resy = 0.0;
  resz = 0.0;

  for (E_Int i = 0; i < ni - 1; i++)
  {
    ind1 = i;
    ind2 = i + 1;
    ind = i;

    rlind = length[ind] * ratio[ind];

    centerx = K_CONST::ONE_HALF * (xt[ind1] + xt[ind2]) - cx;
    centery = K_CONST::ONE_HALF * (yt[ind1] + yt[ind2]) - cy;
    centerz = K_CONST::ONE_HALF * (zt[ind1] + zt[ind2]) - cz;

    mx = centery * vz[ind] - centerz * vy[ind];
    my = centerz * vx[ind] - centerx * vz[ind];
    mz = centerx * vy[ind] - centery * vx[ind];

    resx += rlind * mx;
    resy += rlind * my;
    resz += rlind * mz;
  }

  result[0] = resx;
  result[1] = resy;
  result[2] = resz;
}
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
# include "Array/Array.h"

//=============================================================================
// Integre "surfaciquement" les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integMomentStruct2D(E_Int ni, E_Int nj, E_Int nk, 
                     E_Int center2node, 
                     E_Int posx, E_Int posy, E_Int posz,
                     E_Float cx, E_Float cy, E_Float cz, 
                     FldArrayF& coord, 
                     FldArrayF& F, FldArrayF& ratio, 
                     FldArrayF& resultat)
{
  E_Int NI, NJ;
  FldArrayF res(3);
  
  if      (nk == 1) {NI = ni; NJ = nj;}
  else if (nj == 1) {NI = ni; NJ = nk;}
  else if (ni == 1) {NI = nj; NJ = nk;}
  else return 0;

  
  E_Int ncells = (NI-1)*(NJ-1);
  FldArrayF surf(ncells);
  
  // Compute surface of each "block" i cell, with coordinates coord
  K_METRIC::compSurfStruct2D(
    NI, NJ, 1,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    surf.begin());

  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node 
    // and field F in center 
    integMomentStructCellCenter2D(
      NI, NJ, cx, cy, cz, ratio.begin(), 
      coord.begin(posx), coord.begin(posy), coord.begin(posz),
      surf.begin(), F.begin(1), F.begin(2), F.begin(3), 
      res.begin()
    );
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    integMomentStructNodeCenter2D(
      NI, NJ, cx, cy, cz, ratio.begin(), 
      coord.begin(posx), coord.begin(posy), coord.begin(posz),
      surf.begin(), F.begin(1), F.begin(2), F.begin(3), 
      res.begin()
    );
  }
  
  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];   
  return 1;
}

//=============================================================================
// Integre "lineairement" les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integMomentStruct1D(E_Int ni, E_Int nj, E_Int nk, 
                       E_Int center2node, 
                       E_Int posx, E_Int posy, E_Int posz,
                       E_Float cx, E_Float cy, E_Float cz, 
                       FldArrayF& coord, 
                       FldArrayF& F, FldArrayF& ratio, 
                       FldArrayF& resultat)
{
  E_Int NI;
  FldArrayF res(3);
  resultat.setAllValuesAtNull();

  if      (ni > 1) NI = ni;
  else if (nj > 1) NI = nj;
  else if (nk > 1) NI = nk;
  else return 0;
  
  FldArrayF length(NI-1);
  // Compute surface of each "block" i cell, with coordinates coord
  K_METRIC::compSurfStruct1D(
    NI, 1, 1,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    length.begin());
  
  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node and field F in center 
    integMomentStructCellCenter1D(
      NI, cx, cy, cz, ratio.begin(), 
      coord.begin(posx), coord.begin(posy), coord.begin(posz), 
      length.begin(), F.begin(1), F.begin(2), F.begin(3), 
      res.begin()
    );
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    integMomentStructNodeCenter1D(
      NI, cx, cy, cz, ratio.begin(),
      coord.begin(posx), coord.begin(posy), coord.begin(posz),
      length.begin(), F.begin(1), F.begin(2), F.begin(3),
      res.begin()
    );
  }
  
  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];
  
  return 1;
}
  
// ============================================================================
// Compute surface integral of the moment M (OM^F), coordinates 
//     and field have the same size
//     I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//     Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
// ============================================================================
void K_POST::integMomentStructNodeCenter2D(const E_Int ni, const E_Int nj,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  E_Int ni1 = ni - 1;

  #pragma omp parallel for collapse(2) reduction(+:res1,res2,res3)
  for (E_Int j = 0; j < nj - 1; j++)
  {
    for (E_Int i = 0; i < ni - 1; i++)
    {
      E_Int ind1 = i + j*ni;
      E_Int ind2 = ind1 + ni;
      E_Int ind3 = ind1 + 1;
      E_Int ind4 = ind3 + ni;
      E_Int ind = i + j*ni1;

      E_Float rind, dx, dy, dz;

      rind = ratio[ind1];
      dx = xt[ind1] - cx;
      dy = yt[ind1] - cy;
      dz = zt[ind1] - cz;

      E_Float f1x = rind * vx[ind1];
      E_Float f1y = rind * vy[ind1];
      E_Float f1z = rind * vz[ind1];
      E_Float m1x = dy * f1z - dz * f1y;
      E_Float m1y = dz * f1x - dx * f1z;
      E_Float m1z = dx * f1y - dy * f1x;

      rind = ratio[ind2];
      dx = xt[ind2] - cx;
      dy = yt[ind2] - cy;
      dz = zt[ind2] - cz;

      E_Float f2x = rind * vx[ind2];
      E_Float f2y = rind * vy[ind2];
      E_Float f2z = rind * vz[ind2];
      E_Float m2x = dy * f2z - dz * f2y;
      E_Float m2y = dz * f2x - dx * f2z;
      E_Float m2z = dx * f2y - dy * f2x;

      rind = ratio[ind3];
      dx = xt[ind3] - cx;
      dy = yt[ind3] - cy;
      dz = zt[ind3] - cz;

      E_Float f3x = rind * vx[ind3];
      E_Float f3y = rind * vy[ind3];
      E_Float f3z = rind * vz[ind3];
      E_Float m3x = dy * f3z - dz * f3y;
      E_Float m3y = dz * f3x - dx * f3z;
      E_Float m3z = dx * f3y - dy * f3x;

      rind = ratio[ind4];
      dx = xt[ind4] - cx;
      dy = yt[ind4] - cy;
      dz = zt[ind4] - cz;

      E_Float f4x = rind * vx[ind4];
      E_Float f4y = rind * vy[ind4];
      E_Float f4z = rind * vz[ind4];
      E_Float m4x = dy * f4z - dz * f4y;
      E_Float m4y = dz * f4x - dx * f4z;
      E_Float m4z = dx * f4y - dy * f4x;

      E_Float sind = surf[ind];
      res1 += sind * (m1x + m2x + m3x + m4x);
      res2 += sind * (m1y + m2y + m3y + m4y);
      res3 += sind * (m1z + m2z + m3z + m4z);
    }
  }

  result[0] = 0.25 * res1;
  result[1] = 0.25 * res2;
  result[2] = 0.25 * res3;
}

//=============================================================================
// Compute linear integral of the moment M (OM^F), coordinates 
//     and field have the same size
//     I(AB) = LENGTH(ABCD)*(F(A)+F(B))/2
// ============================================================================
void K_POST::integMomentStructNodeCenter1D(
  const E_Int ni, const E_Float cx, const E_Float cy,
  const E_Float cz, const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* length, const E_Float* vx, const E_Float* vy,
  const E_Float* vz, E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  for (E_Int i = 0; i < ni - 1; i++)
  {
    E_Int ind1 = i;
    E_Int ind2 = i + 1;
    E_Int ind = i;

    E_Float rind, dx, dy, dz;

    rind = ratio[ind1];
    E_Float f1x = rind * vx[ind1];
    E_Float f1y = rind * vy[ind1];
    E_Float f1z = rind * vz[ind1];

    dx = xt[ind1] - cx;
    dy = yt[ind1] - cy;
    dz = zt[ind1] - cz;

    E_Float m1x = dy * f1z - dz * f1y;
    E_Float m1y = dz * f1x - dx * f1z;
    E_Float m1z = dx * f1y - dy * f1x;

    rind = ratio[ind2];
    E_Float f2x = rind * vx[ind2];
    E_Float f2y = rind * vy[ind2];
    E_Float f2z = rind * vz[ind2];

    dx = xt[ind2] - cx;
    dy = yt[ind2] - cy;
    dz = zt[ind2] - cz;

    E_Float m2x = dy * f2z - dz * f2y;
    E_Float m2y = dz * f2x - dx * f2z;
    E_Float m2z = dx * f2y - dy * f2x;

    E_Float lind = length[ind];
    res1 += lind * (m1x + m2x);
    res2 += lind * (m1y + m2y);
    res3 += lind * (m1z + m2z);
  }

  result[0] = 0.5 * res1;
  result[1] = 0.5 * res2;
  result[2] = 0.5 * res3;
}

//=============================================================================
// Compute surface integral of the moment M (OM^F), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
void K_POST::integMomentStructCellCenter2D(const E_Int ni, const E_Int nj,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* surf, const E_Float* vx,
  const E_Float* vy, const E_Float* vz, E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  E_Int ni1 = ni - 1;

  #pragma omp parallel for collapse(2) reduction(+:res1,res2,res3)
  for (E_Int j = 0; j < nj - 1; j++)
  {
    for (E_Int i = 0; i < ni - 1; i++)
    {
      E_Int ind1 = i + j * ni;
      E_Int ind2 = ind1 + ni;
      E_Int ind3 = ind1 + 1;
      E_Int ind4 = ind2 + 1;
      E_Int ind = i + j*ni1;

      E_Float srind = surf[ind]*ratio[ind];
      
      E_Float centerx = xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4];
      centerx = 0.25 * centerx - cx;

      E_Float centery = yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4];
      centery = 0.25 * centery - cy;

      E_Float centerz = zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4];
      centerz = 0.25 * centerz - cz;

      E_Float mx = centery * vz[ind] - centerz * vy[ind];
      E_Float my = centerz * vx[ind] - centerx * vz[ind];
      E_Float mz = centerx * vy[ind] - centery * vx[ind];

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
void K_POST::integMomentStructCellCenter1D(
  const E_Int ni, const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* length, const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  for (E_Int i = 0; i < ni - 1; i++)
  {
    E_Int ind1 = i;
    E_Int ind2 = i + 1;
    E_Int ind = i;

    E_Float rlind = length[ind] * ratio[ind];

    E_Float centerx = 0.5 * (xt[ind1] + xt[ind2]) - cx;
    E_Float centery = 0.5 * (yt[ind1] + yt[ind2]) - cy;
    E_Float centerz = 0.5 * (zt[ind1] + zt[ind2]) - cz;

    E_Float mx = centery * vz[ind] - centerz * vy[ind];
    E_Float my = centerz * vx[ind] - centerx * vz[ind];
    E_Float mz = centerx * vy[ind] - centery * vx[ind];

    res1 += rlind * mx;
    res2 += rlind * my;
    res3 += rlind * mz;
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}
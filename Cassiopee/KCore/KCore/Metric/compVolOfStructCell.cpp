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
# include "metric.h"
# include "Math/math.h"


//=============================================================================
// Calcule l aire d'une cellule d un maillage surfacique nk=1. N'est pas necessairement dans le plan 
// On rentre soit l indice de la cellule indcell, soit indnode l indice du premier point
// d indices i et j min de la cellule. Si indnode est different de -1, c'est lui qui prime
// IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
//=============================================================================
void K_METRIC::compVolOfStructCell2D(
  const E_Int ni, const E_Int nj,
  const E_Int indcell, E_Int indnode,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& area
)
{
  E_Int nic = K_FUNC::E_max(E_Int(1.), ni - 1);

  if (indcell > -1)
  {
    E_Int i, j;
    if (indnode < 0)
    {
      j = indcell / nic;
      i = indcell - j * nic;
      indnode = i + j * ni;
    }
  }
  else
  {
    if (indnode < 0)
    {
      fprintf(stderr, "ERROR: compVolOfStructCell2D:\n");
      fprintf(stderr, "One of indcell or indnode must be > -1.\n");
      exit(0);
    }
  }

  E_Int ind1 = indnode;
  E_Int ind2 = ind1 + ni;
  E_Int ind3 = ind1 + 1;
  E_Int ind4 = ind2 + 1;

  E_Float l1x = xt[ind1] - xt[ind2];
  E_Float l1y = yt[ind1] - yt[ind2];
  E_Float l1z = zt[ind1] - zt[ind2];
  E_Float l2x = xt[ind1] - xt[ind3];
  E_Float l2y = yt[ind1] - yt[ind3];
  E_Float l2z = zt[ind1] - zt[ind3];

  E_Float surf1x = (l1y * l2z - l1z * l2y);
  E_Float surf1y = (l1z * l2x - l1x * l2z);
  E_Float surf1z = (l1x * l2y - l1y * l2x);
  E_Float surface1 = sqrt(surf1x * surf1x + surf1y * surf1y + surf1z * surf1z);
  
  E_Float l3x = xt[ind4] - xt[ind3];
  E_Float l3y = yt[ind4] - yt[ind3];
  E_Float l3z = zt[ind4] - zt[ind3];
  E_Float l4x = xt[ind4] - xt[ind2];
  E_Float l4y = yt[ind4] - yt[ind2];
  E_Float l4z = zt[ind4] - zt[ind2]; 
  
  E_Float surf2x = (l3y * l4z - l3z * l4y);
  E_Float surf2y = (l3z * l4x - l3x * l4z);
  E_Float surf2z = (l3x * l4y - l3y * l4x);
  E_Float surface2 = sqrt(surf2x * surf2x + surf2y * surf2y + surf2z * surf2z);
  
  E_Float ps = surf1x * surf2x + surf1y * surf2y + surf1z * surf2z;
  area = K_CONST::ONE_HALF * K_FUNC::E_sign(ps) * (surface1 + surface2);
}

void K_METRIC::compVolOfStructCell3D(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Int indcell, E_Int indnode,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& vol
)
{
  E_Int ninj, i, j, k, nicnjc, nic, njc;
  E_Int indA, indB, indC, indD, indE, indF, indG, indH;
  E_Float xA, xB, xC, xD, xE, xF, xG, xH;
  E_Float yA, yB, yC, yD, yE, yF, yG, yH;
  E_Float zA, zB, zC, zD, zE, zF, zG, zH;
  E_Float xv1, yv1, zv1, xv2, yv2, zv2;
  E_Float xcc, ycc, zcc;
  E_Float nx, ny, nz;
  E_Float v1, v2, v3, v4, v5, v6;
  const E_Float ONE_24TH = K_CONST::ONE_EIGHTH * K_CONST::ONE_THIRD;

  ninj = ni * nj;
  nic = (ni > 1) ? ni - 1 : 1;
  njc = (nj > 1) ? nj - 1 : 1;
  nicnjc = nic * njc;

  // Warning: here i, j, k start at 0
  if (indcell > -1)
  {
    if (indnode < 0)
    {
      k = indcell / nicnjc;
      j = (indcell - k * nicnjc) / nic;
      i = indcell - j * nic - k * nicnjc;
      indnode = i + j * ni + k * ninj;
    }
  }
  else
  {
    if (indnode < 0)
    {
      fprintf(stderr, "ERROR: compVolOfStructCell3D:\n");
      fprintf(stderr, "One of indcell or indnode must be > -1.\n");
      exit(0);
    }
  }

  k = indnode / ninj;
  j = indnode / ni - k * nj;
  i = indnode - k * ninj - j * ni;

  // Decalages a gauche
  if (i == ni-1) indnode = indnode - 1;
  if (j == nj-1) indnode = indnode - ni;
  if (k == nk-1) indnode = indnode - ninj;
  
  indA = indnode;
  indB = indA + 1;
  indC = indB + ni;
  indD = indA + ni;
  indE = indA + ninj;
  indF = indB + ninj;
  indG = indC + ninj;
  indH = indD + ninj;

  xA = xt[indA]; yA = yt[indA]; zA = zt[indA];
  xB = xt[indB]; yB = yt[indB]; zB = zt[indB];
  xC = xt[indC]; yC = yt[indC]; zC = zt[indC];
  xD = xt[indD]; yD = yt[indD]; zD = zt[indD];
  xE = xt[indE]; yE = yt[indE]; zE = zt[indE];
  xF = xt[indF]; yF = yt[indF]; zF = zt[indF];
  xG = xt[indG]; yG = yt[indG]; zG = zt[indG];
  xH = xt[indH]; yH = yt[indH]; zH = zt[indH];

  // Face ABCD
  xcc = xA + xB + xC + xD;
  ycc = yA + yB + yC + yD;
  zcc = zA + zB + zC + zD;
  xv1 = xC - xA; yv1 = yC - yA; zv1 = zC - zA;
  xv2 = xD - xB; yv2 = yD - yB; zv2 = zD - zB;
  K_MATH::cross(xv1, yv1, zv1, xv2, yv2, zv2, nx, ny, nz);
  v1 = xcc * nx + ycc * ny + zcc * nz;

  // Face EFGH
  xcc = xE + xF + xG + xH;
  ycc = yE + yF + yG + yH;
  zcc = zE + zF + zG + zH;
  xv1 = xG - xE; yv1 = yG - yE; zv1 = zG - zE;
  xv2 = xH - xF; yv2 = yH - yF; zv2 = zH - zF;
  K_MATH::cross(xv1, yv1, zv1, xv2, yv2, zv2, nx, ny, nz);
  v2 = xcc * nx + ycc * ny + zcc * nz;

  // Face EADH
  xcc = xA + xD + xH + xE;
  ycc = yA + yD + yH + yE;
  zcc = zA + zD + zH + zE;
  xv1 = xD - xE; yv1 = yD - yE; zv1 = zD - zE;
  xv2 = xH - xA; yv2 = yH - yA; zv2 = zH - zA;
  K_MATH::cross(xv1, yv1, zv1, xv2, yv2, zv2, nx, ny, nz);
  v3 = xcc * nx + ycc * ny + zcc * nz;

  // Face FBCG
  xcc = xB + xF + xG + xC;
  ycc = yB + yF + yG + yC;
  zcc = zB + zF + zG + zC;
  xv1 = xC - xF; yv1 = yC - yF; zv1 = zC - zF;
  xv2 = xG - xB; yv2 = yG - yB; zv2 = zG - zB;
  K_MATH::cross(xv1, yv1, zv1, xv2, yv2, zv2, nx, ny, nz);
  v4 = xcc * nx + ycc * ny + zcc * nz;

  // Face EFBA
  xcc = xA + xB + xF + xE;
  ycc = yA + yB + yF + yE;
  zcc = zA + zB + zF + zE;
  xv1 = xB - xE; yv1 = yB - yE; zv1 = zB - zE;
  xv2 = xA - xF; yv2 = yA - yF; zv2 = zA - zF;
  K_MATH::cross(xv1, yv1, zv1, xv2, yv2, zv2, nx, ny, nz);
  v5 = xcc * nx + ycc * ny + zcc * nz;

  // Face HGCD
  xcc = xD + xC + xG + xH;
  ycc = yD + yC + yG + yH;
  zcc = zD + zC + zG + zH;
  xv1 = xC - xH; yv1 = yC - yH; zv1 = zC - zH;
  xv2 = xD - xG; yv2 = yD - yG; zv2 = zD - zG;
  K_MATH::cross(xv1, yv1, zv1, xv2, yv2, zv2, nx, ny, nz);
  v6 = xcc * nx + ycc * ny + zcc * nz;

  vol = K_FUNC::E_abs(v2 - v1 + v4 - v3 + v6 - v5) * ONE_24TH;
}
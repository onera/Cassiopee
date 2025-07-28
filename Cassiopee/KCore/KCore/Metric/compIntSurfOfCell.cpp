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
// Calcul des vecteurs surfaces des 6 interfaces d'une cellule 3D structuree.
// IN: ind: indice du premier sommet de la cellule
// IN: ni, nj, nk: dimensions du maillage en noeuds
// IN: xt, yt, zt: Vertex coordinates
// OUT: surf: vecteur surface de la cellule
//=============================================================================
void K_METRIC::compIntSurfOfCell(
  const E_Int ind, const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surf
)
{
  E_Int i, j, k;
  E_Int ninj;
  E_Int l1, l2, l3, l4;
  E_Float d13x, d13y, d13z, d24x, d24y, d24z;

  ninj = ni * nj;

  k = ind / ninj;
  j = ind / ni - k * nj;
  i = ind - j * ni - k * ninj;

  i = i + 1;
  j = j + 1;
  k = k + 1;

  // interface in i
  l1 = ind;                  // i,j,k
  l2 = l1 + ni;              // i,j+1,k
  l3 = l2 + ninj;            // i,j+1,k+1
  l4 = l1 + ninj;            // i,j,k+1

  d13x = xt[l3] - xt[l1];
  d13y = yt[l3] - yt[l1];
  d13z = zt[l3] - zt[l1];

  d24x = xt[l4] - xt[l2];
  d24y = yt[l4] - yt[l2];
  d24z = zt[l4] - zt[l2];

  surf[0*3+0] = -K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
  surf[0*3+1] = -K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
  surf[0*3+2] = -K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

  // interface in i+1
  l1 = ind + 1;              // i+1,j,k
  l2 = l1 + ni;              // i+1,j+1,k
  l3 = l2 + ninj;            // i+1,j+1,k+1
  l4 = l1 + ninj;            // i+1,j,k+1

  d13x = xt[l3] - xt[l1];
  d13y = yt[l3] - yt[l1];
  d13z = zt[l3] - zt[l1];

  d24x = xt[l4] - xt[l2];
  d24y = yt[l4] - yt[l2];
  d24z = zt[l4] - zt[l2];

  surf[1*3+0] = K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
  surf[1*3+1] = K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
  surf[1*3+2] = K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

  // interface in j
  l1 = ind + 1;              // i+1,j,k
  l2 = l1 - 1;               // i,j,k
  l3 = l2 + ninj;            // i,j,k+1
  l4 = l1 + ninj;            // i+1,j,k+1

  d13x = xt[l3] - xt[l1];
  d13y = yt[l3] - yt[l1];
  d13z = zt[l3] - zt[l1];

  d24x = xt[l4] - xt[l2];
  d24y = yt[l4] - yt[l2];
  d24z = zt[l4] - zt[l2];

  surf[2*3+0] = -K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
  surf[2*3+1] = -K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
  surf[2*3+2] = -K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

  // interface in j+1
  l1 = ind + ni + 1;         // i+1,j+1,k
  l2 = l1 - 1;               // i,j+1,k
  l3 = l2 + ninj;            // i,j+1,k+1
  l4 = l1 + ninj;            // i,j+1,k+1

  d13x = xt[l3] - xt[l1];
  d13y = yt[l3] - yt[l1];
  d13z = zt[l3] - zt[l1];

  d24x = xt[l4] - xt[l2];
  d24y = yt[l4] - yt[l2];
  d24z = zt[l4] - zt[l2];

  surf[3*3+0] = K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
  surf[3*3+1] = K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
  surf[3*3+2] = K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

  // interface in k
  l1 = ind;                  // i,j,k
  l2 = l1 + 1;               // i+1,j,k
  l3 = l2 + ni;              // i+1,j+1,k
  l4 = l1 + ni;              // i,j+1,k

  d13x = xt[l3] - xt[l1];
  d13y = yt[l3] - yt[l1];
  d13z = zt[l3] - zt[l1];

  d24x = xt[l4] - xt[l2];
  d24y = yt[l4] - yt[l2];
  d24z = zt[l4] - zt[l2];

  surf[4*3+0] = -K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
  surf[4*3+1] = -K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
  surf[4*3+2] = -K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

  // interface in k+1
  l1 = ind + ninj;           // i,j,k+1
  l2 = l1 + 1;               // i+1,j,k+1
  l3 = l2 + ni;              // i+1,j+1,k+1
  l4 = l1 + ni;              // i,j+1,k+1

  d13x = xt[l3] - xt[l1];
  d13y = yt[l3] - yt[l1];
  d13z = zt[l3] - zt[l1];

  d24x = xt[l4] - xt[l2];
  d24y = yt[l4] - yt[l2];
  d24z = zt[l4] - zt[l2];

  surf[5*3+0] = K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
  surf[5*3+1] = K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
  surf[5*3+2] = K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);
}
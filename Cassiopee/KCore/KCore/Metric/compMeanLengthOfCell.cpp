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
#include <math.h>


//=============================================================================
// Calcul d'une longueur caracteristique d'une cellule structuree: 
// moyenne des longueurs des cotes de la cellule
// IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
// IN: indA: First vertex of the cell
// IN: xt, yt, zt: Vertex coordinates
// OUT: meanl: Cell mean length
//=============================================================================
void K_METRIC::compMeanLengthOfStructCell(
  const E_Int ni, const E_Int nj, const E_Int nk, E_Int indA,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& meanl
)
{
  const E_Float ONE_12TH = K_CONST::ONE_FOURTH * K_CONST::ONE_THIRD;
  E_Int ninj = ni * nj;

  // Warning: here i, j, k start at 0
  E_Int k = indA / ninj;
  E_Int j = indA / ni - k * nj;
  E_Int i = indA - k * ninj - j * ni;

  // Increments
  E_Int alpha = 1;
  E_Int beta = 1;
  E_Int gamma = 1;

  if (i == ni - 1) i -= 1;
  if (j == nj - 1) j -= 1;
  if (k == nk - 1) k -= 1;

  // 2D case
  if (nk == 1) gamma = 0;
  if (ni == 1) alpha = 0;
  if (nj == 1) beta = 0;

  indA = i + j * ni + k * ninj;
  E_Int indB = (i + alpha) + j * ni + k * ninj;
  E_Int indC = (i + alpha) + (j + beta) * ni + k * ninj;
  E_Int indD = i + (j + beta) * ni + k * ninj;

  E_Int indE = i + j * ni + (k + gamma) * ninj;
  E_Int indF = (i + alpha) + j * ni + (k + gamma) * ninj;
  E_Int indG = (i + alpha) + (j + beta) * ni + (k + gamma) * ninj;
  E_Int indH = i + (j + beta) * ni + (k + gamma) * ninj;

  // premier centre: facette ABCD
  E_Float x21 = xt[indE] + xt[indF] + xt[indG] + xt[indH]
              - (xt[indA] + xt[indB] + xt[indC] + xt[indD]);
  E_Float y21 = yt[indE] + yt[indF] + yt[indG] + yt[indH]
              - (yt[indA] + yt[indB] + yt[indC] + yt[indD]);
  E_Float z21 = zt[indE] + zt[indF] + zt[indG] + zt[indH]
              - (zt[indA] + zt[indB] + zt[indC] + zt[indD]);

  E_Float x43 = xt[indA] + xt[indD] + xt[indE] + xt[indH]
              - (xt[indB] + xt[indC] + xt[indF] + xt[indG]);
  E_Float y43 = yt[indA] + yt[indD] + yt[indE] + yt[indH]
              - (yt[indB] + yt[indC] + yt[indF] + yt[indG]);
  E_Float z43 = zt[indA] + zt[indD] + zt[indE] + zt[indH]
              - (zt[indB] + zt[indC] + zt[indF] + zt[indG]);

  E_Float x65 = xt[indG] + xt[indH] + xt[indC] + xt[indD]
              - (xt[indA] + xt[indB] + xt[indF] + xt[indE]);
  E_Float y65 = yt[indG] + yt[indH] + yt[indC] + yt[indD]
              - (yt[indA] + yt[indB] + yt[indF] + yt[indE]);
  E_Float z65 = zt[indG] + zt[indH] + zt[indC] + zt[indD]
              - (zt[indA] + zt[indB] + zt[indF] + zt[indE]);

  E_Float d1 = sqrt(x21 * x21 + y21 * y21 + z21 * z21);
  E_Float d2 = sqrt(x43 * x43 + y43 * y43 + z43 * z43);
  E_Float d3 = sqrt(x65 * x65 + y65 * y65 + z65 * z65);

  meanl = (d1 + d2 + d3) * ONE_12TH;
}

//=============================================================================
// Calcul d'une longueur caracteristique d'une cellule tetra: 
// moyenne des longueurs des cotes de la cellule
// IN: indA, indB, indC, indD: Sommets du tetraedre
// IN: xt, yt, zt: Vertex coordinates
// OUT: meanl: Cell mean length
//=============================================================================
void K_METRIC::compMeanLengthOfTetraCell(
  const E_Int indA, const E_Int indB, const E_Int indC, const E_Int indD,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& meanl
)
{
  const E_Float ONE_SIXTH = K_CONST::ONE_HALF * K_CONST::ONE_THIRD;
  E_Float xA = xt[indA], yA = yt[indA], zA = zt[indA];
  E_Float xB = xt[indB], yB = yt[indB], zB = zt[indB];
  E_Float xC = xt[indC], yC = yt[indC], zC = zt[indC];
  E_Float xD = xt[indD], yD = yt[indD], zD = zt[indD];

  E_Float dAB = sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA) + (zB - zA) * (zB - zA));
  E_Float dAC = sqrt((xC - xA) * (xC - xA) + (yC - yA) * (yC - yA) + (zC - zA) * (zC - zA));
  E_Float dDA = sqrt((xA - xD) * (xA - xD) + (yA - yD) * (yA - yD) + (zA - zD) * (zA - zD));
  E_Float dBC = sqrt((xC - xB) * (xC - xB) + (yC - yB) * (yC - yB) + (zC - zB) * (zC - zB));
  E_Float dBD = sqrt((xD - xB) * (xD - xB) + (yD - yB) * (yD - yB) + (zD - zB) * (zD - zB));
  E_Float dCD = sqrt((xD - xC) * (xD - xC) + (yD - yC) * (yD - yC) + (zD - zC) * (zD - zC));

  meanl = (dAB + dAC + dBC + dCD + dDA + dBD) * ONE_SIXTH;
}
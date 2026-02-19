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
# include <math.h>

//=============================================================================
// Calcul d'une longueur caracteristique de la cellule : minimum  des
// longueurs des cotes de la cellule
// IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
// IN: indA: Index of the first cell vertex
// IN: xt, yt, zt: Vertex coordinates
// OUT: minl: Minimum length of the cell
//=============================================================================
void K_METRIC::compMinLengthOfCell(
  const E_Int ni, const E_Int nj, const E_Int nk, E_Int indA,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& minl
)
{
  E_Int ninj, i, j, k, alpha, beta, gamma;
  E_Int indB, indC, indD, indE, indF, indG, indH;
  E_Float xA, xB, xC, xD, xE, xF, xG, xH;
  E_Float yA, yB, yC, yD, yE, yF, yG, yH;
  E_Float zA, zB, zC, zD, zE, zF, zG, zH;
  E_Float dAB, dCD, dEF, dGH, dAE, dBF, dDH, dCG, dAD, dBC, dEH, dFG;

  ninj = ni * nj;

  // Warning: here i, j, k start at 0
  k = indA / ninj;
  j = indA / ni - k * nj;
  i = indA - k * ninj - j * ni;

  // Increments
  alpha = 1;
  beta = 1;
  gamma = 1;

  if (i == ni - 1) i -= 1;
  if (j == nj - 1) j -= 1;
  if (k == nk - 1) k -= 1;

  // cas 2D
  if (ni == 1) alpha = 0;
  if (nj == 1) beta = 0;
  if (nk == 1) gamma = 0;

  indA = i + j * ni + k * ninj;
  indB = (i + alpha) + j * ni + k * ninj;
  indC = (i + alpha) + (j + beta) * ni + k * ninj;
  indD = i + (j + beta) * ni + k * ninj;
  indE = i + j * ni + (k + gamma) * ninj;
  indF = (i + alpha) + j * ni + (k + gamma) * ninj;
  indG = (i + alpha) + (j + beta) * ni + (k + gamma) * ninj;
  indH = i + (j + beta) * ni + (k + gamma) * ninj;

  xA = xt[indA]; yA = yt[indA]; zA = zt[indA];
  xB = xt[indB]; yB = yt[indB]; zB = zt[indB];
  xC = xt[indC]; yC = yt[indC]; zC = zt[indC];
  xD = xt[indD]; yD = yt[indD]; zD = zt[indD];
  xE = xt[indE]; yE = yt[indE]; zE = zt[indE];
  xF = xt[indF]; yF = yt[indF]; zF = zt[indF];
  xG = xt[indG]; yG = yt[indG]; zG = zt[indG];
  xH = xt[indH]; yH = yt[indH]; zH = zt[indH];

  if (alpha == 0)
  {
    dAB = K_CONST::E_MAX_FLOAT;
    dCD = K_CONST::E_MAX_FLOAT;
    dEF = K_CONST::E_MAX_FLOAT;
    dGH = K_CONST::E_MAX_FLOAT;
  }
  else
  {
    dAB = sqrt((xB - xA)*(xB - xA) + (yB - yA)*(yB - yA) + (zB - zA)*(zB - zA));
    dCD = sqrt((xD - xC)*(xD - xC) + (yD - yC)*(yD - yC) + (zD - zC)*(zD - zC));
    dEF = sqrt((xF - xE)*(xF - xE) + (yF - yE)*(yF - yE) + (zF - zE)*(zF - zE));
    dGH = sqrt((xH - xG)*(xH - xG) + (yH - yG)*(yH - yG) + (zH - zG)*(zH - zG));
  }

  if (beta == 0)
  {
    dAD = K_CONST::E_MAX_FLOAT;
    dBC = K_CONST::E_MAX_FLOAT;
    dEH = K_CONST::E_MAX_FLOAT;
    dFG = K_CONST::E_MAX_FLOAT;
  }
  else
  {
    dAD = sqrt((xD - xA)*(xD - xA) + (yD - yA)*(yD - yA) + (zD - zA)*(zD - zA));
    dBC = sqrt((xC - xB)*(xC - xB) + (yC - yB)*(yC - yB) + (zC - zB)*(zC - zB));
    dEH = sqrt((xH - xE)*(xH - xE) + (yH - yE)*(yH - yE) + (zH - zE)*(zH - zE));
    dFG = sqrt((xG - xF)*(xG - xF) + (yG - yF)*(yG - yF) + (zG - zF)*(zG - zF));
  }

  if (gamma == 0)
  {
    dAE = K_CONST::E_MAX_FLOAT;
    dBF = K_CONST::E_MAX_FLOAT;
    dDH = K_CONST::E_MAX_FLOAT;
    dCG = K_CONST::E_MAX_FLOAT;
  }
  else
  {
    dAE = sqrt((xE - xA)*(xE - xA) + (yE - yA)*(yE - yA) + (zE - zA)*(zE - zA));
    dBF = sqrt((xF - xB)*(xF - xB) + (yF - yB)*(yF - yB) + (zF - zB)*(zF - zB));
    dDH = sqrt((xH - xD)*(xH - xD) + (yH - yD)*(yH - yD) + (zH - zD)*(zH - zD));
    dCG = sqrt((xG - xC)*(xG - xC) + (yG - yC)*(yG - yC) + (zG - zC)*(zG - zC));
  }

  minl = dAB;
  minl = K_FUNC::E_min(minl, dEF);
  minl = K_FUNC::E_min(minl, dGH);
  minl = K_FUNC::E_min(minl, dAD);
  minl = K_FUNC::E_min(minl, dBC);
  minl = K_FUNC::E_min(minl, dEH);
  minl = K_FUNC::E_min(minl, dFG);
  minl = K_FUNC::E_min(minl, dAE);
  minl = K_FUNC::E_min(minl, dBF);
  minl = K_FUNC::E_min(minl, dDH);
  minl = K_FUNC::E_min(minl, dCG);
  minl = K_FUNC::E_min(minl, dCD);
}

//=============================================================================
// Calcul d'une longueur caracteristique de la cellule : minimum  des
// longueurs des cotes de la cellule
// IN: indA, indB, indC, indD: Sommets du tetraedre
// IN: xt, yt, zt: Vertex coordinates
// OUT: minl: Minimum length of the cell
//=============================================================================
void K_METRIC::compMinLengthOfTetraCell(
  const E_Int indA, const E_Int indB, const E_Int indC, const E_Int indD,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& minl
)
{
  E_Float xA, xB, xC, xD;
  E_Float yA, yB, yC, yD;
  E_Float zA, zB, zC, zD;
  E_Float dAB, dAC, dDA, dBC, dBD, dCD;

  xA = xt[indA]; yA = yt[indA]; zA = zt[indA];
  xB = xt[indB]; yB = yt[indB]; zB = zt[indB];
  xC = xt[indC]; yC = yt[indC]; zC = zt[indC];
  xD = xt[indD]; yD = yt[indD]; zD = zt[indD];

  dAB = sqrt((xB - xA)*(xB - xA) + (yB - yA)*(yB - yA) + (zB - zA)*(zB - zA));
  dAC = sqrt((xC - xA)*(xC - xA) + (yC - yA)*(yC - yA) + (zC - zA)*(zC - zA));
  dDA = sqrt((xA - xD)*(xA - xD) + (yA - yD)*(yA - yD) + (zA - zD)*(zA - zD));
  dBC = sqrt((xC - xB)*(xC - xB) + (yC - yB)*(yC - yB) + (zC - zB)*(zC - zB));
  dBD = sqrt((xD - xB)*(xD - xB) + (yD - yB)*(yD - yB) + (zD - zB)*(zD - zB));
  dCD = sqrt((xD - xC)*(xD - xC) + (yD - yC)*(yD - yC) + (zD - zC)*(zD - zC));

  minl = dAB;
  minl = K_FUNC::E_min(minl, dDA);
  minl = K_FUNC::E_min(minl, dBC);
  minl = K_FUNC::E_min(minl, dBD);
  minl = K_FUNC::E_min(minl, dCD);
  minl = K_FUNC::E_min(minl, dAC);
}

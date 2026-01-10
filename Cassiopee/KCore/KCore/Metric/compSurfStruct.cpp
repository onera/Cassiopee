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
// Calcul de la surface pour une grille surfacique structuree
// IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
// IN: xt, yt, zt: Vertex coordinates
// OUT: surface: Mesh area
//=============================================================================
void K_METRIC::compSurfStruct2D(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surface
)
{
  E_Int ni1, nj1, nk1;
  E_Int ind1, ind2, ind3, ind4, indcell;
  E_Float l1x, l1y, l1z;
  E_Float l2x, l2y, l2z;
  E_Float l3x, l3y, l3z;
  E_Float l4x, l4y, l4z;
  E_Float surf1x, surf1y, surf1z;
  E_Float surf2x, surf2y, surf2z;
  E_Float surface1, surface2, ps;
  
  ni1 = ni - 1;
  nj1 = nj - 1;
  nk1 = nk - 1;

  if (ni == 1)
  {
    for (E_Int k = 0; k < nk1; k++)
    {
      for (E_Int j = 0; j < nj1; j++)
      {
        indcell = j + k * nj1;
        ind1 = j + k * nj;
        ind2 = ind1 + nj;
        ind3 = ind1 + 1;
        ind4 = ind2 + 1;
        l1x = xt[ind1] - xt[ind2];
        l1y = yt[ind1] - yt[ind2];
        l1z = zt[ind1] - zt[ind2];
        l2x = xt[ind1] - xt[ind3];
        l2y = yt[ind1] - yt[ind3];
        l2z = zt[ind1] - zt[ind3];
        surf1x = l1y * l2z - l1z * l2y;
        surf1y = l1z * l2x - l1x * l2z;
        surf1z = l1x * l2y - l1y * l2x;
        surface1 = sqrt(surf1x * surf1x + surf1y * surf1y + surf1z * surf1z);
        l3x = xt[ind4] - xt[ind3];
        l3y = yt[ind4] - yt[ind3];
        l3z = zt[ind4] - zt[ind3];
        l4x = xt[ind4] - xt[ind2];
        l4y = yt[ind4] - yt[ind2];
        l4z = zt[ind4] - zt[ind2];
        surf2x = l3y * l4z - l3z * l4y;
        surf2y = l3z * l4x - l3x * l4z;
        surf2z = l3x * l4y - l3y * l4x;
        surface2 = sqrt(surf2x * surf2x + surf2y * surf2y + surf2z * surf2z);
        ps = surf1x * surf2x + surf1y * surf2y + surf1z * surf2z;
        surface[indcell] = K_CONST::ONE_HALF * K_FUNC::E_sign(ps) * (surface1 + surface2);
      }
    }
  }
  else if (nj == 1)
  {
    for (E_Int i = 0; i < ni1; i++)
    {
      for (E_Int k = 0; k < nk1; k++)
      {
        indcell = i + k * ni1;
        ind1 = i + k * ni;
        ind2 = ind1 + ni;
        ind3 = ind1 + 1;
        ind4 = ind2 + 1;
        l1x = xt[ind1] - xt[ind2];
        l1y = yt[ind1] - yt[ind2];
        l1z = zt[ind1] - zt[ind2];
        l2x = xt[ind1] - xt[ind3];
        l2y = yt[ind1] - yt[ind3];
        l2z = zt[ind1] - zt[ind3];
        surf1x = l1y * l2z - l1z * l2y;
        surf1y = l1z * l2x - l1x * l2z;
        surf1z = l1x * l2y - l1y * l2x;
        surface1 = sqrt(surf1x * surf1x + surf1y * surf1y + surf1z * surf1z);
        l3x = xt[ind4] - xt[ind3];
        l3y = yt[ind4] - yt[ind3];
        l3z = zt[ind4] - zt[ind3];
        l4x = xt[ind4] - xt[ind2];
        l4y = yt[ind4] - yt[ind2];
        l4z = zt[ind4] - zt[ind2];
        surf2x = l3y * l4z - l3z * l4y;
        surf2y = l3z * l4x - l3x * l4z;
        surf2z = l3x * l4y - l3y * l4x;
        surface2 = sqrt(surf2x * surf2x + surf2y * surf2y + surf2z * surf2z);
        ps = surf1x * surf2x + surf1y * surf2y + surf1z * surf2z;
        surface[indcell] = K_CONST::ONE_HALF * K_FUNC::E_sign(ps) * (surface1 + surface2);
      }
    }
  }
  else if (nk == 1)
  {
    for (E_Int j = 0; j < nj1; j++)
    {
      for (E_Int i = 0; i < ni1; i++)
      {
        indcell = i + j * ni1;
        ind1 = i + j * ni;
        ind2 = ind1 + ni;
        ind3 = ind1 + 1;
        ind4 = ind2 + 1;
        l1x = xt[ind1] - xt[ind2];
        l1y = yt[ind1] - yt[ind2];
        l1z = zt[ind1] - zt[ind2];
        l2x = xt[ind1] - xt[ind3];
        l2y = yt[ind1] - yt[ind3];
        l2z = zt[ind1] - zt[ind3];
        surf1x = l1y * l2z - l1z * l2y;
        surf1y = l1z * l2x - l1x * l2z;
        surf1z = l1x * l2y - l1y * l2x;
        surface1 = sqrt(surf1x * surf1x + surf1y * surf1y + surf1z*surf1z);
        l3x = xt[ind4] - xt[ind3];
        l3y = yt[ind4] - yt[ind3];
        l3z = zt[ind4] - zt[ind3];
        l4x = xt[ind4] - xt[ind2];
        l4y = yt[ind4] - yt[ind2];
        l4z = zt[ind4] - zt[ind2];
        surf2x = l3y * l4z - l3z * l4y;
        surf2y = l3z * l4x - l3x * l4z;
        surf2z = l3x * l4y - l3y * l4x;
        surface2 = sqrt(surf2x * surf2x + surf2y * surf2y + surf2z * surf2z);
        ps = surf1x * surf2x + surf1y * surf2y + surf1z * surf2z;
        surface[indcell] = K_CONST::ONE_HALF * K_FUNC::E_sign(ps) * (surface1 + surface2);
      }
    }
  }
  else
  {
    fprintf(stderr, "In compSurfStruct: 2D field is required.\n");
    exit(0);
  }
}


//=============================================================================
// Calcul de la longueur entre chaque sommet pour une ligne structuree.
// IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
// IN: xt, yt, zt: Vertex coordinates
// OUT: length: Longueur entre chaque sommet
//=============================================================================
void K_METRIC::compSurfStruct1D(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* length
)
{
  E_Int i, j, k;
  E_Int ind1, ind2, indcell;
  E_Float l1x, l1y, l1z;

  if ((nj == 1) && (nk == 1))
  {
    for (i = 0; i < ni - 1; i++)
    {
      indcell = i;
      ind1 = i;
      ind2 = i + 1;
      l1x = xt[ind1] - xt[ind2];
      l1y = yt[ind1] - yt[ind2];
      l1z = zt[ind1] - zt[ind2];
      length[indcell] = sqrt(l1x * l1x + l1y * l1y + l1z * l1z);
    }
  }
  else if ((ni == 1) && (nk == 1))
  {
    for (j = 0; j < nj - 1; j++)
    {
      indcell = j;
      ind1 = j;
      ind2 = j + 1;
      l1x = xt[ind1] - xt[ind2];
      l1y = yt[ind1] - yt[ind2];
      l1z = zt[ind1] - zt[ind2];
      length[indcell] = sqrt(l1x * l1x + l1y * l1y + l1z * l1z);
    }
  }
  else if ((ni == 1) && (nj == 1))
  {
    for (k = 0; k < nk - 1; k++)
    {
      indcell = k;
      ind1 = k;
      ind2 = k + 1;
      l1x = xt[ind1] - xt[ind2];
      l1y = yt[ind1] - yt[ind2];
      l1z = zt[ind1] - zt[ind2];
      length[indcell] = sqrt(l1x * l1x + l1y * l1y + l1z * l1z);
    }
  }
  else
  {
    fprintf(stderr, "In compSurfStruct1D: 1D field is required.\n");
    exit(0);
  }
}


//=============================================================================
// Calcul de la surface pour une grille surfacique structuree.
// IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
// IN: xt, yt, zt: Vertex coordinates
// OUT: surface: Aire
//=============================================================================
void K_METRIC::compSurfOfStructCell(
  const E_Int ni, const E_Int nj, const E_Int nk, const E_Int indcell,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& surface
)
{
  E_Int i, j, k, ni1, nj1, nk1;
  E_Int ind1, ind2, ind3, ind4;
  E_Float l1x, l1y, l1z;
  E_Float l2x, l2y, l2z;
  E_Float l3x, l3y, l3z;
  E_Float l4x, l4y, l4z;
  E_Float surf1x, surf1y, surf1z;
  E_Float surf2x, surf2y, surf2z;
  E_Float surface1, surface2, ps;

  ni1 = ni - 1;
  nj1 = nj - 1;
  nk1 = nk - 1;

  i = indcell % ni;
  j = (indcell % (ni*nj)) / ni;
  k = indcell / (ni*nj);

  if (i == ni1) i = i - 1;
  if (j == nj1) j = j - 1;
  if (k == nk1) k = k - 1;

  if (ni == 1)
  {
    ind1 = j + k*nj;
    ind2 = ind1 + nj;
    ind3 = ind1 + 1;
    ind4 = ind2 + 1;

    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];
    l2x = xt[ind1] - xt[ind3];
    l2y = yt[ind1] - yt[ind3];
    l2z = zt[ind1] - zt[ind3];

    surf1x = l1y * l2z - l1z * l2y;
    surf1y = l1z * l2x - l1x * l2z;
    surf1z = l1x * l2y - l1y * l2x;

    surface1 = sqrt(surf1x * surf1x + surf1y * surf1y + surf1z * surf1z);

    l3x = xt[ind4] - xt[ind3];
    l3y = yt[ind4] - yt[ind3];
    l3z = zt[ind4] - zt[ind3];
    l4x = xt[ind4] - xt[ind2];
    l4y = yt[ind4] - yt[ind2];
    l4z = zt[ind4] - zt[ind2];

    surf2x = l3y * l4z - l3z * l4y;
    surf2y = l3z * l4x - l3x * l4z;
    surf2z = l3x * l4y - l3y * l4x;

    surface2 = sqrt(surf2x * surf2x + surf2y * surf2y + surf2z * surf2z);
    ps = surf1x * surf2x + surf1y * surf2y + surf1z * surf2z;

    surface = K_CONST::ONE_HALF * K_FUNC::E_sign(ps) * (surface1 + surface2);
  }
  else if (nj == 1)
  {
    ind1 = i + k*ni;
    ind2 = ind1 + ni;
    ind3 = ind1 + 1;
    ind4 = ind2 + 1;

    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];
    l2x = xt[ind1] - xt[ind3];
    l2y = yt[ind1] - yt[ind3];
    l2z = zt[ind1] - zt[ind3];

    surf1x = l1y * l2z - l1z * l2y;
    surf1y = l1z * l2x - l1x * l2z;
    surf1z = l1x * l2y - l1y * l2x;

    surface1 = sqrt(surf1x * surf1x + surf1y * surf1y + surf1z * surf1z);

    l3x = xt[ind4] - xt[ind3];
    l3y = yt[ind4] - yt[ind3];
    l3z = zt[ind4] - zt[ind3];
    l4x = xt[ind4] - xt[ind2];
    l4y = yt[ind4] - yt[ind2];
    l4z = zt[ind4] - zt[ind2];

    surf2x = l3y * l4z - l3z * l4y;
    surf2y = l3z * l4x - l3x * l4z;
    surf2z = l3x * l4y - l3y * l4x;

    surface2 = sqrt(surf2x * surf2x + surf2y * surf2y + surf2z * surf2z);
    ps = surf1x * surf2x + surf1y * surf2y + surf1z * surf2z;

    surface = K_CONST::ONE_HALF * K_FUNC::E_sign(ps) * (surface1 + surface2);
  }
  else if (nk == 1)
  {
    ind1 = i + j*ni;
    ind2 = ind1 + ni;
    ind3 = ind1 + 1;
    ind4 = ind2 + 1;

    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];
    l2x = xt[ind1] - xt[ind3];
    l2y = yt[ind1] - yt[ind3];
    l2z = zt[ind1] - zt[ind3];

    surf1x = l1y * l2z - l1z * l2y;
    surf1y = l1z * l2x - l1x * l2z;
    surf1z = l1x * l2y - l1y * l2x;

    surface1 = sqrt(surf1x * surf1x + surf1y * surf1y + surf1z * surf1z);

    l3x = xt[ind4] - xt[ind3];
    l3y = yt[ind4] - yt[ind3];
    l3z = zt[ind4] - zt[ind3];
    l4x = xt[ind4] - xt[ind2];
    l4y = yt[ind4] - yt[ind2];
    l4z = zt[ind4] - zt[ind2];

    surf2x = l3y * l4z - l3z * l4y;
    surf2y = l3z * l4x - l3x * l4z;
    surf2z = l3x * l4y - l3y * l4x;

    surface2 = sqrt(surf2x * surf2x + surf2y * surf2y + surf2z * surf2z);
    ps = surf1x * surf2x + surf1y * surf2y + surf1z * surf2z;

    surface = K_CONST::ONE_HALF * K_FUNC::E_sign(ps) * (surface1 + surface2);
  }
  else
  {
    fprintf(stderr, "In compSurfOfStructCell: 2D field is required.\n");
    exit(0);
  }
}


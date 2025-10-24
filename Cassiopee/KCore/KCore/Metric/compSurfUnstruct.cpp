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
# include "Array/Array.h"
#include "Connect/connect.h"

//=============================================================================
// Calcul des normales pour un maillage Multi-Elements.
// Les normales aux surfaces sont orientees vers l'exterieur de l'element.
// IN: xt, yt, zt: pointeurs sur les coordonnees du maillage
//=============================================================================
void K_METRIC::compSurfUnstruct(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface
)
{
  E_Int fctOffset = 0;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Number of facets per element
  std::vector<E_Int> nfpe;
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
  if (ierr != 0) return;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int nfpc;  // number of facets per connectivity

    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      compTriSurf(cm, fctOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
    }
    else if (strcmp(eltTypes[ic], "QUAD") == 0)
    {
      compQuadSurf(cm, fctOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
    }
    else if (strcmp(eltTypes[ic], "TETRA") == 0)
    {
      compTetraSurf(cm, fctOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
    }
    else if (strcmp(eltTypes[ic], "PYRA") == 0)
    {
      compPyraSurf(cm, fctOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
    }
    else if (strcmp(eltTypes[ic], "PENTA") == 0)
    {
      compPentaSurf(cm, fctOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
    }
    else if (strcmp(eltTypes[ic], "HEXA") == 0)
    {
      compHexaSurf(cm, fctOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
    }
    else
    {
      fprintf(stderr, "Error: in K_METRIC::compUnstrSurf.\n");
      fprintf(stderr, "Unknown type of element.\n");
      exit(0);
    }

    // Increment the face offset
    nfpc = nfpe[ic]*nelts;
    fctOffset += nfpc;
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

void K_METRIC::compTriSurf(
  K_FLD::FldArrayI& cm, const E_Int fctOffset,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface
)
{
  E_Int nelts = cm.getSize();
  E_Int ind1, ind2, ind3, pos;
  E_Float surf, surfx, surfy, surfz;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z;

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cm(i, 1) - 1;
    ind2 = cm(i, 2) - 1;
    ind3 = cm(i, 3) - 1;

    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];

    l2x = xt[ind1] - xt[ind3];
    l2y = yt[ind1] - yt[ind3];
    l2z = zt[ind1] - zt[ind3];

    surfx = (l1y * l2z - l1z * l2y);
    surfy = (l1z * l2x - l1x * l2z);
    surfz = (l1x * l2y - l1y * l2x);
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos = fctOffset + i;  // fctOffset + i = fctOffset + i * nfpe
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;
  }
}

void K_METRIC::compQuadSurf(
  K_FLD::FldArrayI& cm, const E_Int fctOffset,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface
)
{
  E_Int nelts = cm.getSize();
  E_Int ind1, ind2, ind3, ind4, pos;
  E_Float surf, surfx, surfy, surfz;
  E_Float surf1x, surf1y, surf1z, surf2x, surf2y, surf2z;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z;

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cm(i, 1) - 1;
    ind2 = cm(i, 2) - 1;
    ind3 = cm(i, 3) - 1;
    ind4 = cm(i, 4) - 1;

    // AB x AC
    l1x = xt[ind2] - xt[ind1];
    l1y = yt[ind2] - yt[ind1];
    l1z = zt[ind2] - zt[ind1];

    l2x = xt[ind3] - xt[ind1];
    l2y = yt[ind3] - yt[ind1];
    l2z = zt[ind3] - zt[ind1];

    surf1x = (l1y * l2z - l1z * l2y);
    surf1y = (l1z * l2x - l1x * l2z);
    surf1z = (l1x * l2y - l1y * l2x);

    // AC x AD
    l1x = xt[ind3] - xt[ind1];
    l1y = yt[ind3] - yt[ind1];
    l1z = zt[ind3] - zt[ind1];

    l2x = xt[ind4] - xt[ind1];
    l2y = yt[ind4] - yt[ind1];
    l2z = zt[ind4] - zt[ind1];

    surf2x = (l1y * l2z - l1z * l2y);
    surf2y = (l1z * l2x - l1x * l2z);
    surf2z = (l1x * l2y - l1y * l2x);

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos = fctOffset + i;  // fctOffset + i = fctOffset + i * nfpe
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;
  }
}

void K_METRIC::compTetraSurf(
  K_FLD::FldArrayI& cm, const E_Int fctOffset,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface
)
{
  E_Int nelts = cm.getSize();
  E_Int nfpe = 4;
  E_Int ind1, ind2, ind3, ind4, pos;
  E_Float surf, surfx, surfy, surfz;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z;

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cm(i, 1) - 1;  // A1
    ind2 = cm(i, 2) - 1;  // A2
    ind3 = cm(i, 3) - 1;  // A3
    ind4 = cm(i, 4) - 1;  // A4

    // Face A1A2A3
    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];
    l2x = xt[ind3] - xt[ind2];
    l2y = yt[ind3] - yt[ind2];
    l2z = zt[ind3] - zt[ind2];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos = fctOffset + i * nfpe;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Face A1A2A4
    l1x = xt[ind2] - xt[ind1];
    l1y = yt[ind2] - yt[ind1];
    l1z = zt[ind2] - zt[ind1];
    l2x = xt[ind4] - xt[ind1];
    l2y = yt[ind4] - yt[ind1];
    l2z = zt[ind4] - zt[ind1];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Face A2A3A4
    l1x = xt[ind3] - xt[ind2];
    l1y = yt[ind3] - yt[ind2];
    l1z = zt[ind3] - zt[ind2];
    l2x = xt[ind4] - xt[ind2];
    l2y = yt[ind4] - yt[ind2];
    l2z = zt[ind4] - zt[ind2];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Face A1A3A4
    l1x = xt[ind4] - xt[ind1];
    l1y = yt[ind4] - yt[ind1];
    l1z = zt[ind4] - zt[ind1];
    l2x = xt[ind3] - xt[ind1];
    l2y = yt[ind3] - yt[ind1];
    l2z = zt[ind3] - zt[ind1];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;  
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;
  }
}

void K_METRIC::compPyraSurf(
  K_FLD::FldArrayI& cm, const E_Int fctOffset,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface
)
{
  E_Int nelts = cm.getSize();
  E_Int nfpe = 5;
  E_Int indt1, indt2, indt3;
  E_Int ind1, ind2, ind3, ind4, ind5, pos;
  E_Float surf, surfx, surfy, surfz;
  E_Float surf1x, surf1y, surf1z, surf2x, surf2y, surf2z;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z;

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cm(i, 1) - 1;  // A1
    ind2 = cm(i, 2) - 1;  // A2
    ind3 = cm(i, 3) - 1;  // A3
    ind4 = cm(i, 4) - 1;  // A4
    ind5 = cm(i, 5) - 1;  // A5

    // First face: A1A2A3A4
    // A2A1 x A2A3
    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];
    l2x = xt[ind3] - xt[ind2];
    l2y = yt[ind3] - yt[ind2];
    l2z = zt[ind3] - zt[ind2];
    surf1x = l1y * l2z - l1z * l2y;
    surf1y = l1z * l2x - l1x * l2z;
    surf1z = l1x * l2y - l1y * l2x;

    // A4A3 x A4A1
    l1x = xt[ind3] - xt[ind4];
    l1y = yt[ind3] - yt[ind4];
    l1z = zt[ind3] - zt[ind4];
    l2x = xt[ind1] - xt[ind4];
    l2y = yt[ind1] - yt[ind4];
    l2z = zt[ind1] - zt[ind4];
    surf2x = l1y * l2z - l1z * l2y;
    surf2y = l1z * l2x - l1x * l2z;
    surf2z = l1x * l2y - l1y * l2x;

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos = fctOffset + i * nfpe;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Second face: triangle A1A2A5
    indt1 = ind1;
    indt2 = ind2;
    indt3 = ind5;
    l1x = xt[indt1] - xt[indt2];
    l1y = yt[indt1] - yt[indt2];
    l1z = zt[indt1] - zt[indt2];
    l2x = xt[indt1] - xt[indt3];
    l2y = yt[indt1] - yt[indt3];
    l2z = zt[indt1] - zt[indt3];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Third face: triangle A2A3A5
    indt1 = ind2;
    indt2 = ind3;
    indt3 = ind5;
    l1x = xt[indt1] - xt[indt2];
    l1y = yt[indt1] - yt[indt2];
    l1z = zt[indt1] - zt[indt2];
    l2x = xt[indt1] - xt[indt3];
    l2y = yt[indt1] - yt[indt3];
    l2z = zt[indt1] - zt[indt3];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Fourth face: triangle A3A4A5
    indt1 = ind3;
    indt2 = ind4;
    indt3 = ind5;
    l1x = xt[indt1] - xt[indt2];
    l1y = yt[indt1] - yt[indt2];
    l1z = zt[indt1] - zt[indt2];
    l2x = xt[indt1] - xt[indt3];
    l2y = yt[indt1] - yt[indt3];
    l2z = zt[indt1] - zt[indt3];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Fifth face: triangle A4A1A5
    indt1 = ind4;
    indt2 = ind1;
    indt3 = ind5;
    l1x = xt[indt1] - xt[indt2];
    l1y = yt[indt1] - yt[indt2];
    l1z = zt[indt1] - zt[indt2];
    l2x = xt[indt1] - xt[indt3];
    l2y = yt[indt1] - yt[indt3];
    l2z = zt[indt1] - zt[indt3];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;
  }
}

void K_METRIC::compPentaSurf(
  K_FLD::FldArrayI& cm, const E_Int fctOffset,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface
)
{
  E_Int nelts = cm.getSize();
  E_Int nfpe = 5;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, pos;
  E_Float surf, surfx, surfy, surfz;
  E_Float surf1x, surf1y, surf1z, surf2x, surf2y, surf2z;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z;

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cm(i, 1) - 1;  // A1
    ind2 = cm(i, 2) - 1;  // A2
    ind3 = cm(i, 3) - 1;  // A3
    ind4 = cm(i, 4) - 1;  // A4
    ind5 = cm(i, 5) - 1;  // A5
    ind6 = cm(i, 6) - 1;  // A6

    // First face: triangle A1A2A3
    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];
    l2x = xt[ind3] - xt[ind2];
    l2y = yt[ind3] - yt[ind2];
    l2z = zt[ind3] - zt[ind2];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos = fctOffset + i * nfpe;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Second face: triangle A4A5A6
    l1x = xt[ind5] - xt[ind4];
    l1y = yt[ind5] - yt[ind4];
    l1z = zt[ind5] - zt[ind4];
    l2x = xt[ind6] - xt[ind4];
    l2y = yt[ind6] - yt[ind4];
    l2z = zt[ind6] - zt[ind4];
    surfx = l1y * l2z - l1z * l2y;
    surfy = l1z * l2x - l1x * l2z;
    surfz = l1x * l2y - l1y * l2x;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Third face: quad 1254
    l1x = xt[ind2] - xt[ind1];
    l1y = yt[ind2] - yt[ind1];
    l1z = zt[ind2] - zt[ind1];
    l2x = xt[ind5] - xt[ind1];
    l2y = yt[ind5] - yt[ind1];
    l2z = zt[ind5] - zt[ind1];
    surf1x = l1y * l2z - l1z * l2y;
    surf1y = l1z * l2x - l1x * l2z;
    surf1z = l1x * l2y - l1y * l2x;

    l1x = xt[ind5] - xt[ind1];
    l1y = yt[ind5] - yt[ind1];
    l1z = zt[ind5] - zt[ind1];
    l2x = xt[ind4] - xt[ind1];
    l2y = yt[ind4] - yt[ind1];
    l2z = zt[ind4] - zt[ind1];
    surf2x = l1y * l2z - l1z * l2y;
    surf2y = l1z * l2x - l1x * l2z;
    surf2z = l1x * l2y - l1y * l2x;

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Fourth face: quad 2365
    l1x = xt[ind3] - xt[ind2];
    l1y = yt[ind3] - yt[ind2];
    l1z = zt[ind3] - zt[ind2];
    l2x = xt[ind6] - xt[ind2];
    l2y = yt[ind6] - yt[ind2];
    l2z = zt[ind6] - zt[ind2];
    surf1x = l1y * l2z - l1z * l2y;
    surf1y = l1z * l2x - l1x * l2z;
    surf1z = l1x * l2y - l1y * l2x;

    l1x = xt[ind6] - xt[ind2];
    l1y = yt[ind6] - yt[ind2];
    l1z = zt[ind6] - zt[ind2];
    l2x = xt[ind5] - xt[ind2];
    l2y = yt[ind5] - yt[ind2];
    l2z = zt[ind5] - zt[ind2];
    surf2x = l1y * l2z - l1z * l2y;
    surf2y = l1z * l2x - l1x * l2z;
    surf2z = l1x * l2y - l1y * l2x;

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // Fifth face: quad 3146
    l1x = xt[ind1] - xt[ind3];
    l1y = yt[ind1] - yt[ind3];
    l1z = zt[ind1] - zt[ind3];
    l2x = xt[ind4] - xt[ind3];
    l2y = yt[ind4] - yt[ind3];
    l2z = zt[ind4] - zt[ind3];
    surf1x = l1y * l2z - l1z * l2y;
    surf1y = l1z * l2x - l1x * l2z;
    surf1z = l1x * l2y - l1y * l2x;

    l1x = xt[ind4] - xt[ind3];
    l1y = yt[ind4] - yt[ind3];
    l1z = zt[ind4] - zt[ind3];
    l2x = xt[ind6] - xt[ind3];
    l2y = yt[ind6] - yt[ind3];
    l2z = zt[ind6] - zt[ind3];
    surf2x = l1y * l2z - l1z * l2y;
    surf2y = l1z * l2x - l1x * l2z;
    surf2z = l1x * l2y - l1y * l2x;

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;
  }
}

void K_METRIC::compHexaSurf(
  K_FLD::FldArrayI& cm, const E_Int fctOffset,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface
)
{
  E_Int nelts = cm.getSize();
  E_Int nfpe = 6;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, pos;
  E_Float surf, surfx, surfy, surfz;
  E_Float surf1x, surf1y, surf1z, surf2x, surf2y, surf2z;
  E_Float l1x, l1y, l1z, l2x, l2y, l2z;

  for (E_Int i = 0; i < nelts; i++)
  {
    ind1 = cm(i, 1) - 1;  // A1
    ind2 = cm(i, 2) - 1;  // A2
    ind3 = cm(i, 3) - 1;  // A3
    ind4 = cm(i, 4) - 1;  // A4
    ind5 = cm(i, 5) - 1;  // A5
    ind6 = cm(i, 6) - 1;  // A6
    ind7 = cm(i, 7) - 1;  // A7
    ind8 = cm(i, 8) - 1;  // A8

    // premiere facette A1A2A3A4
    // A2A1 x A2A3
    l1x = xt[ind1] - xt[ind2];
    l1y = yt[ind1] - yt[ind2];
    l1z = zt[ind1] - zt[ind2];

    l2x = xt[ind3] - xt[ind2];
    l2y = yt[ind3] - yt[ind2];
    l2z = zt[ind3] - zt[ind2];

    surf1x = (l1y * l2z - l1z * l2y);
    surf1y = (l1z * l2x - l1x * l2z);
    surf1z = (l1x * l2y - l1y * l2x);

    // A4A3 x A4A1
    l1x = xt[ind3] - xt[ind4];
    l1y = yt[ind3] - yt[ind4];
    l1z = zt[ind3] - zt[ind4];

    l2x = xt[ind1] - xt[ind4];
    l2y = yt[ind1] - yt[ind4];
    l2z = zt[ind1] - zt[ind4];

    surf2x = (l1y * l2z - l1z * l2y);
    surf2y = (l1z * l2x - l1x * l2z);
    surf2z = (l1x * l2y - l1y * l2x);

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos = fctOffset + i * nfpe;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // deuxieme facette A5A6A7A8
    // A5A6 x A5A7
    l1x = xt[ind6] - xt[ind5];
    l1y = yt[ind6] - yt[ind5];
    l1z = zt[ind6] - zt[ind5];

    l2x = xt[ind7] - xt[ind5];
    l2y = yt[ind7] - yt[ind5];
    l2z = zt[ind7] - zt[ind5];

    surf1x = (l1y * l2z - l1z * l2y);
    surf1y = (l1z * l2x - l1x * l2z);
    surf1z = (l1x * l2y - l1y * l2x);

    // A5A7 x A5A8
    l1x = xt[ind7] - xt[ind5];
    l1y = yt[ind7] - yt[ind5];
    l1z = zt[ind7] - zt[ind5];

    l2x = xt[ind8] - xt[ind5];
    l2y = yt[ind8] - yt[ind5];
    l2z = zt[ind8] - zt[ind5];

    surf2x = (l1y * l2z - l1z * l2y);
    surf2y = (l1z * l2x - l1x * l2z);
    surf2z = (l1x * l2y - l1y * l2x);

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // troisieme facette 4158
    // A4A1 x A4A5
    l1x = xt[ind1] - xt[ind4];
    l1y = yt[ind1] - yt[ind4];
    l1z = zt[ind1] - zt[ind4];

    l2x = xt[ind5] - xt[ind4];
    l2y = yt[ind5] - yt[ind4];
    l2z = zt[ind5] - zt[ind4];

    surf1x = (l1y * l2z - l1z * l2y);
    surf1y = (l1z * l2x - l1x * l2z);
    surf1z = (l1x * l2y - l1y * l2x);

    // A4A5 x A4A8
    l1x = xt[ind5] - xt[ind4];
    l1y = yt[ind5] - yt[ind4];
    l1z = zt[ind5] - zt[ind4];

    l2x = xt[ind8] - xt[ind4];
    l2y = yt[ind8] - yt[ind4];
    l2z = zt[ind8] - zt[ind4];

    surf2x = (l1y * l2z - l1z * l2y);
    surf2y = (l1z * l2x - l1x * l2z);
    surf2z = (l1x * l2y - l1y * l2x);

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // quatrieme facette A2A3A7A6
    // A2A3x A2A7
    l1x = xt[ind3] - xt[ind2];
    l1y = yt[ind3] - yt[ind2];
    l1z = zt[ind3] - zt[ind2];

    l2x = xt[ind7] - xt[ind2];
    l2y = yt[ind7] - yt[ind2];
    l2z = zt[ind7] - zt[ind2];

    surf1x = (l1y * l2z - l1z * l2y);
    surf1y = (l1z * l2x - l1x * l2z);
    surf1z = (l1x * l2y - l1y * l2x);

    // A2A7 x A2A6
    l1x = xt[ind7] - xt[ind2];
    l1y = yt[ind7] - yt[ind2];
    l1z = zt[ind7] - zt[ind2];

    l2x = xt[ind6] - xt[ind2];
    l2y = yt[ind6] - yt[ind2];
    l2z = zt[ind6] - zt[ind2];

    surf2x = (l1y * l2z - l1z * l2y);
    surf2y = (l1z * l2x - l1x * l2z);
    surf2z = (l1x * l2y - l1y * l2x);

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // cinquieme facette A1A2A6A5
    // A1A2 x A1A6
    l1x = xt[ind2] - xt[ind1];
    l1y = yt[ind2] - yt[ind1];
    l1z = zt[ind2] - zt[ind1];

    l2x = xt[ind6] - xt[ind1];
    l2y = yt[ind6] - yt[ind1];
    l2z = zt[ind6] - zt[ind1];

    surf1x = (l1y * l2z - l1z * l2y);
    surf1y = (l1z * l2x - l1x * l2z);
    surf1z = (l1x * l2y - l1y * l2x);

    // A1A6 x A1A5
    l1x = xt[ind6] - xt[ind1];
    l1y = yt[ind6] - yt[ind1];
    l1z = zt[ind6] - zt[ind1];

    l2x = xt[ind5] - xt[ind1];
    l2y = yt[ind5] - yt[ind1];
    l2z = zt[ind5] - zt[ind1];

    surf2x = (l1y * l2z - l1z * l2y);
    surf2y = (l1z * l2x - l1x * l2z);
    surf2z = (l1x * l2y - l1y * l2x);

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;

    // sixieme facette A3A4A8A7
    // A3A4 x A3A8
    l1x = xt[ind4] - xt[ind3];
    l1y = yt[ind4] - yt[ind3];
    l1z = zt[ind4] - zt[ind3];

    l2x = xt[ind8] - xt[ind3];
    l2y = yt[ind8] - yt[ind3];
    l2z = zt[ind8] - zt[ind3];

    surf1x = (l1y * l2z - l1z * l2y);
    surf1y = (l1z * l2x - l1x * l2z);
    surf1z = (l1x * l2y - l1y * l2x);

    // A3A8 x A3A7
    l1x = xt[ind8] - xt[ind3];
    l1y = yt[ind8] - yt[ind3];
    l1z = zt[ind8] - zt[ind3];

    l2x = xt[ind7] - xt[ind3];
    l2y = yt[ind7] - yt[ind3];
    l2z = zt[ind7] - zt[ind3];

    surf2x = (l1y * l2z - l1z * l2y);
    surf2y = (l1z * l2x - l1x * l2z);
    surf2z = (l1x * l2y - l1y * l2x);

    surfx = surf1x + surf2x;
    surfy = surf1y + surf2y;
    surfz = surf1z + surf2z;
    surf = sqrt(surfx * surfx + surfy * surfy + surfz * surfz);

    pos += 1;
    surfnx[pos] = K_CONST::ONE_HALF * surfx;
    surfny[pos] = K_CONST::ONE_HALF * surfy;
    surfnz[pos] = K_CONST::ONE_HALF * surfz;
    surface[pos] = K_CONST::ONE_HALF * surf;
  }
}

//=============================================================================
// Calcul de la longueur entre chaque sommet pour une ligne non structuree
//=============================================================================
void K_METRIC::compUnstructSurf1d(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* length
)
{
  E_Int offset = 0;
  E_Int ind1, ind2;
  E_Float lx, ly, lz;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    
    if (strcmp(eltTypes[ic], "BAR") == 0)
    {
      E_Int nelts = cm.getSize();
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        lx = xt[ind1] - xt[ind2];
        ly = yt[ind1] - yt[ind2];
        lz = zt[ind1] - zt[ind2];
        length[offset+i] = sqrt(lx * lx + ly * ly + lz * lz);
      }
      offset += nelts;
    }
    else
    {
      fprintf(stderr, "Error in K_METRIC::compUnstructSurf1d.\n");
      fprintf(stderr, "Element type can be BAR only, not %s.\n", eltTypes[ic]);
      exit(0);
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}
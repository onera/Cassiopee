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

//=============================================================================
// Calcul des centres des interfaces pour des mailles non structurees.
// IN: cn: Element-Node connectivity
// IN: xt, yt, zt: Vertex coordinates
// OUT: xint, yint, zint: Coordonnees du centre des facettes
//=============================================================================
void K_METRIC::compUnstructCenterInt(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float* xint, E_Float* yint, E_Float* zint
)
{
  E_Int fctOffset = 0;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int nfpe;  // number of facets per element
    E_Int nfpc;  // number of facets per connectivity

    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      nfpe = 1;

      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, pos;
        E_Float x1, x2, x3;
        E_Float y1, y2, y3;
        E_Float z1, z2, z3;

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;

          x1 = xt[ind1]; x2 = xt[ind2]; x3 = xt[ind3];
          y1 = yt[ind1]; y2 = yt[ind2]; y3 = yt[ind3];
          z1 = zt[ind1]; z2 = zt[ind2]; z3 = zt[ind3];

          pos = fctOffset + i;  // fctOffset + i = fctOffset + i * nfpe
          xint[pos] = K_CONST::ONE_THIRD * (x1 + x2 + x3);
          yint[pos] = K_CONST::ONE_THIRD * (y1 + y2 + y3);
          zint[pos] = K_CONST::ONE_THIRD * (z1 + z2 + z3);
        }
      }
    }
    else if (strcmp(eltTypes[ic], "QUAD") == 0)
    {
      nfpe = 1;

      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, ind4, pos;
        E_Float x1, x2, x3, x4;
        E_Float y1, y2, y3, y4;
        E_Float z1, z2, z3, z4;
        
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;
          ind4 = cm(i, 4) - 1;

          x1 = xt[ind1]; x2 = xt[ind2]; x3 = xt[ind3]; x4 = xt[ind4];
          y1 = yt[ind1]; y2 = yt[ind2]; y3 = yt[ind3]; y4 = yt[ind4];
          z1 = zt[ind1]; z2 = zt[ind2]; z3 = zt[ind3]; z4 = zt[ind4];

          pos = fctOffset + i;  // fctOffset + i = fctOffset + i * nfpe
          xint[pos] = K_CONST::ONE_FOURTH * (x1 + x2 + x3 + x4);
          yint[pos] = K_CONST::ONE_FOURTH * (y1 + y2 + y3 + y4);
          zint[pos] = K_CONST::ONE_FOURTH * (z1 + z2 + z3 + z4);
        }
      }
    }
    else if (strcmp(eltTypes[ic], "TETRA") == 0)
    {
      nfpe = 4;

      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, ind4, pos;
        E_Float x1, x2, x3, x4;
        E_Float y1, y2, y3, y4;
        E_Float z1, z2, z3, z4;
      
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;
          ind4 = cm(i, 4) - 1;

          x1 = xt[ind1]; x2 = xt[ind2]; x3 = xt[ind3]; x4 = xt[ind4];
          y1 = yt[ind1]; y2 = yt[ind2]; y3 = yt[ind3]; y4 = yt[ind4];
          z1 = zt[ind1]; z2 = zt[ind2]; z3 = zt[ind3]; z4 = zt[ind4];
          
          // facette 123
          pos = fctOffset + i * nfpe;
          xint[pos] = K_CONST::ONE_THIRD * (x1 + x2 + x3);
          yint[pos] = K_CONST::ONE_THIRD * (y1 + y2 + y3);
          zint[pos] = K_CONST::ONE_THIRD * (z1 + z2 + z3);

          // facette 124
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x1 + x2 + x4);
          yint[pos] = K_CONST::ONE_THIRD * (y1 + y2 + y4);
          zint[pos] = K_CONST::ONE_THIRD * (z1 + z2 + z4);

          // facette 234
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x2 + x3 + x4);
          yint[pos] = K_CONST::ONE_THIRD * (y2 + y3 + y4);
          zint[pos] = K_CONST::ONE_THIRD * (z2 + z3 + z4);

          // facette 134
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x1 + x3 + x4);
          yint[pos] = K_CONST::ONE_THIRD * (y1 + y3 + y4);
          zint[pos] = K_CONST::ONE_THIRD * (z1 + z3 + z4);
        }
      }
    }
    else if (strcmp(eltTypes[ic], "PYRA") == 0)
    {
      nfpe = 5;

      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, ind4, ind5, pos;
        E_Float x1, x2, x3, x4, x5;
        E_Float y1, y2, y3, y4, y5;
        E_Float z1, z2, z3, z4, z5;

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;
          ind4 = cm(i, 4) - 1;
          ind5 = cm(i, 5) - 1;

          x1 = xt[ind1]; x2 = xt[ind2]; x3 = xt[ind3]; x4 = xt[ind4]; x5 = xt[ind5];
          y1 = yt[ind1]; y2 = yt[ind2]; y3 = yt[ind3]; y4 = yt[ind4]; y5 = yt[ind5];
          z1 = zt[ind1]; z2 = zt[ind2]; z3 = zt[ind3]; z4 = zt[ind4]; z5 = zt[ind5];

          // facette 1234 : quad
          pos = fctOffset + i * nfpe;
          xint[pos] = K_CONST::ONE_FOURTH * (x1 + x2 + x3 + x4);
          yint[pos] = K_CONST::ONE_FOURTH * (y1 + y2 + y3 + y4);
          zint[pos] = K_CONST::ONE_FOURTH * (z1 + z2 + z3 + z4);

          // facette 125
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x1 + x2 + x5);
          yint[pos] = K_CONST::ONE_THIRD * (y1 + y2 + y5);
          zint[pos] = K_CONST::ONE_THIRD * (z1 + z2 + z5);

          // facette 235
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x2 + x3 + x5);
          yint[pos] = K_CONST::ONE_THIRD * (y2 + y3 + y5);
          zint[pos] = K_CONST::ONE_THIRD * (z2 + z3 + z5);

          // facette 345
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x5 + x3 + x4);
          yint[pos] = K_CONST::ONE_THIRD * (y5 + y3 + y4);
          zint[pos] = K_CONST::ONE_THIRD * (z5 + z3 + z4);

          // facette 415
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x1 + x5 + x4);
          yint[pos] = K_CONST::ONE_THIRD * (y1 + y5 + y4);
          zint[pos] = K_CONST::ONE_THIRD * (z1 + z5 + z4);
        }
      }
    }
    else if (strcmp(eltTypes[ic], "HEXA") == 0)
    {
      nfpe = 6;

      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, pos;
        E_Float x1, x2, x3, x4, x5, x6, x7, x8;
        E_Float y1, y2, y3, y4, y5, y6, y7, y8;
        E_Float z1, z2, z3, z4, z5, z6, z7, z8;

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;
          ind4 = cm(i, 4) - 1;
          ind5 = cm(i, 5) - 1;
          ind6 = cm(i, 6) - 1;
          ind7 = cm(i, 7) - 1;
          ind8 = cm(i, 8) - 1;

          x1 = xt[ind1]; x2 = xt[ind2]; x3 = xt[ind3]; x4 = xt[ind4];
          x5 = xt[ind5]; x6 = xt[ind6]; x7 = xt[ind7]; x8 = xt[ind8];

          y1 = yt[ind1]; y2 = yt[ind2]; y3 = yt[ind3]; y4 = yt[ind4];
          y5 = yt[ind5]; y6 = yt[ind6]; y7 = yt[ind7]; y8 = yt[ind8];

          z1 = zt[ind1]; z2 = zt[ind2]; z3 = zt[ind3]; z4 = zt[ind4];
          z5 = zt[ind5]; z6 = zt[ind6]; z7 = zt[ind7]; z8 = zt[ind8];

          // facette 1234
          pos = fctOffset + i * nfpe;
          xint[pos] = K_CONST::ONE_FOURTH * (x1 + x2 + x3 + x4);
          yint[pos] = K_CONST::ONE_FOURTH * (y1 + y2 + y3 + y4);
          zint[pos] = K_CONST::ONE_FOURTH * (z1 + z2 + z3 + z4);

          // facette 5678
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x5 + x6 + x7 + x8);
          yint[pos] = K_CONST::ONE_FOURTH * (y5 + y6 + y7 + y8);
          zint[pos] = K_CONST::ONE_FOURTH * (z5 + z6 + z7 + z8);

          // facette 1485
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x1 + x4 + x8 + x5);
          yint[pos] = K_CONST::ONE_FOURTH * (y1 + y4 + y8 + y5);
          zint[pos] = K_CONST::ONE_FOURTH * (z1 + z4 + z8 + z5);

          // facette 2376
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x2 + x3 + x7 + x6);
          yint[pos] = K_CONST::ONE_FOURTH * (y2 + y3 + y7 + y6);
          zint[pos] = K_CONST::ONE_FOURTH * (z2 + z3 + z7 + z6);

          // facette 1265
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x1 + x2 + x6 + x5);
          yint[pos] = K_CONST::ONE_FOURTH * (y1 + y2 + y6 + y5);
          zint[pos] = K_CONST::ONE_FOURTH * (z1 + z2 + z6 + z5);

          // facette 4378
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x4 + x3 + x7 + x8);
          yint[pos] = K_CONST::ONE_FOURTH * (y4 + y3 + y7 + y8);
          zint[pos] = K_CONST::ONE_FOURTH * (z4 + z3 + z7 + z8);
        }
      }
    }
    else if (strcmp(eltTypes[ic], "PENTA") == 0)
    {
      nfpe = 5;

      #pragma omp parallel
      {
        E_Int ind1, ind2, ind3, ind4, ind5, ind6, pos;
        E_Float x1, x2, x3, x4, x5, x6;
        E_Float y1, y2, y3, y4, y5, y6;
        E_Float z1, z2, z3, z4, z5, z6;

        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm(i, 1) - 1;
          ind2 = cm(i, 2) - 1;
          ind3 = cm(i, 3) - 1;
          ind4 = cm(i, 4) - 1;
          ind5 = cm(i, 5) - 1;
          ind6 = cm(i, 6) - 1;

          x1 = xt[ind1]; x2 = xt[ind2]; x3 = xt[ind3];
          x4 = xt[ind4]; x5 = xt[ind5]; x6 = xt[ind6];

          y1 = yt[ind1]; y2 = yt[ind2]; y3 = yt[ind3];
          y4 = yt[ind4]; y5 = yt[ind5]; y6 = yt[ind6];

          z1 = zt[ind1]; z2 = zt[ind2]; z3 = zt[ind3];
          z4 = zt[ind4]; z5 = zt[ind5]; z6 = zt[ind6];

          // facette triangle 1 : 123
          pos = fctOffset + i * nfpe;
          xint[pos] = K_CONST::ONE_THIRD * (x1 + x2 + x3);
          yint[pos] = K_CONST::ONE_THIRD * (y1 + y2 + y3);
          zint[pos] = K_CONST::ONE_THIRD * (z1 + z2 + z3);

          // facette triangle 2 : 456
          pos += 1;
          xint[pos] = K_CONST::ONE_THIRD * (x4 + x5 + x6);
          yint[pos] = K_CONST::ONE_THIRD * (y4 + y5 + y6);
          zint[pos] = K_CONST::ONE_THIRD * (z4 + z5 + z6);

          // troisieme facette quad 1254
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x1 + x2 + x5 + x4);
          yint[pos] = K_CONST::ONE_FOURTH * (y1 + y2 + y5 + y4);
          zint[pos] = K_CONST::ONE_FOURTH * (z1 + z2 + z5 + z4);

          // quatrieme facette quad 2365
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x2 + x3 + x6 + x5);
          yint[pos] = K_CONST::ONE_FOURTH * (y2 + y3 + y6 + y5);
          zint[pos] = K_CONST::ONE_FOURTH * (z2 + z3 + z6 + z5);

          // cinquieme facette quad 3146
          pos += 1;
          xint[pos] = K_CONST::ONE_FOURTH * (x3 + x1 + x4 + x6);
          yint[pos] = K_CONST::ONE_FOURTH * (y3 + y1 + y4 + y6);
          zint[pos] = K_CONST::ONE_FOURTH * (z3 + z1 + z4 + z6);
        }
      }
    }
    else
    {
      fprintf(stderr, "Error in K_METRIC::compUnstructCenterInt.\n");
      fprintf(stderr, "Unknown type of element, %s.\n", eltTypes[ic]);
      exit(0);
    }

    // Increment the face offset
    nfpc = nfpe*nelts;
    fctOffset += nfpc;
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}
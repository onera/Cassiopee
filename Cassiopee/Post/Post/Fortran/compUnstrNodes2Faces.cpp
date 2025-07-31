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

void compUnstrNodes2Faces(
  K_FLD::FldArrayI& cn, const E_Float* fieldn,
  E_Float* fieldf
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  E_Int elOffset = 0;

  # pragma omp parallel
  {
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();

      if (strcmp(eltTypes[ic], "TRI") == 0)
      {
        E_Int ind1, ind2, ind3, pos;
        E_Float f1, f2, f3;

        # pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm[i*3] - 1;
          ind2 = cm[i*3 + 1] - 1;
          ind3 = cm[i*3 + 2] - 1;

          f1 = fieldn[ind1];
          f2 = fieldn[ind2];
          f3 = fieldn[ind3];

          pos = elOffset + i * 3;
          fieldf[pos] = K_CONST::ONE_HALF * (f1 + f2);
          fieldf[pos + 1] = K_CONST::ONE_HALF * (f2 + f3);
          fieldf[pos + 2] = K_CONST::ONE_HALF * (f3 + f1);
        }
      }
      else if (strcmp(eltTypes[ic], "QUAD") == 0)
      {
        E_Int ind1, ind2, ind3, ind4, pos;
        E_Float f1, f2, f3, f4;

        # pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm[i*4] - 1;
          ind2 = cm[i*4 + 1] - 1;
          ind3 = cm[i*4 + 2] - 1;
          ind4 = cm[i*4 + 3] - 1;

          f1 = fieldn[ind1];
          f2 = fieldn[ind2];
          f3 = fieldn[ind3];
          f4 = fieldn[ind4];

          pos = elOffset + i * 4;
          fieldf[pos] = K_CONST::ONE_HALF * (f1 + f2);
          fieldf[pos + 1] = K_CONST::ONE_HALF * (f2 + f3);
          fieldf[pos + 2] = K_CONST::ONE_HALF * (f3 + f4);
          fieldf[pos + 3] = K_CONST::ONE_HALF * (f4 + f1);
        }
      }
      else if (strcmp(eltTypes[ic], "TETRA") == 0)
      {
        E_Int ind1, ind2, ind3, ind4, pos;
        E_Float f1, f2, f3, f4;

        # pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm[i*4] - 1;
          ind2 = cm[i*4 + 1] - 1;
          ind3 = cm[i*4 + 2] - 1;
          ind4 = cm[i*4 + 3] - 1;

          f1 = fieldn[ind1];
          f2 = fieldn[ind2];
          f3 = fieldn[ind3];
          f4 = fieldn[ind4];

          pos = elOffset + i * 3;
          fieldf[pos] = K_CONST::ONE_THIRD * (f1 + f2 + f3);      // A1A2A3
          fieldf[pos + 1] = K_CONST::ONE_THIRD * (f1 + f2 + f4);  // A1A2A4
          fieldf[pos + 2] = K_CONST::ONE_THIRD * (f2 + f3 + f4);  // A2A3A4
          fieldf[pos + 3] = K_CONST::ONE_THIRD * (f1 + f3 + f4);  // A1A3A4
        }
      }
      else if (strcmp(eltTypes[ic], "HEXA") == 0)
      {
        nfpe = 6;
        compHexaSurf(cm, fcOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
        k6comphexafield(npts, nelts, cn, fieldn, fieldf);
      }
      else if (strcmp(eltTypes[ic], "PENTA") == 0)
      {
        nfpe = 5;
        compPentaSurf(cm, fcOffset, xt, yt, zt, surfnx, surfny, surfnz, surface);
        k6comppentafield(npts, nelts, cn, fieldn, fieldf);
      }
      else
      {
        fprintf(stderr, "Error: in K_POST::compUnstrNodes2Faces.\n");
        fprintf(stderr, "Unknown type of element.\n");
        exit(0);
      }
      // Increment element offset
      elOffset += nelts;
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

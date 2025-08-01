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

/* Calcul des valeurs d un champ aux faces des elements a partir des 
   des valeurs aux noeuds.
   IN: cn: connectivite elts-noeuds
   IN: fieldn: champs aux noeuds
   OUT: fieldf: champs aux centres des faces
*/
void K_POST::compUnstrNodes2Faces(
  K_FLD::FldArrayI& cn, const char* eltType, const E_Float* fieldn,
  E_Float* fieldf
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Pre-compute element offset
  vector<E_Int> nepc(nc+1); nepc[0] = 0;;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
  }

  # pragma omp parallel
  {
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int elOffset = nepc[ic];

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
      else if (strcmp(eltTypes[ic], "PYRA") == 0)
      {
        E_Int ind1, ind2, ind3, ind4, ind5, pos;
        E_Float f1, f2, f3, f4, f5;

        # pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm[i*5] - 1;
          ind2 = cm[i*5 + 1] - 1;
          ind3 = cm[i*5 + 2] - 1;
          ind4 = cm[i*5 + 3] - 1;
          ind5 = cm[i*5 + 4] - 1;

          f1 = fieldn[ind1];
          f2 = fieldn[ind2];
          f3 = fieldn[ind3];
          f4 = fieldn[ind4];
          f5 = fieldn[ind5];

          pos = elOffset + i * 5;
          fieldf[pos] = K_CONST::ONE_FOURTH * (f1 + f2 + f3 + f4);  // A1A2A3A4
          fieldf[pos + 1] = K_CONST::ONE_THIRD * (f1 + f2 + f5);    // A1A2A5
          fieldf[pos + 2] = K_CONST::ONE_THIRD * (f2 + f3 + f5);    // A2A3A5
          fieldf[pos + 3] = K_CONST::ONE_THIRD * (f3 + f4 + f5);    // A3A4A5
          fieldf[pos + 4] = K_CONST::ONE_THIRD * (f4 + f1 + f5);    // A4A1A5
        }
      }
      
      else if (strcmp(eltTypes[ic], "PENTA") == 0)
      {
        E_Int ind1, ind2, ind3, ind4, ind5, ind6, pos;
        E_Float f1, f2, f3, f4, f5, f6;

        # pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm[i*6] - 1;
          ind2 = cm[i*6 + 1] - 1;
          ind3 = cm[i*6 + 2] - 1;
          ind4 = cm[i*6 + 3] - 1;
          ind5 = cm[i*6 + 4] - 1;
          ind6 = cm[i*6 + 5] - 1;

          f1 = fieldn[ind1];
          f2 = fieldn[ind2];
          f3 = fieldn[ind3];
          f4 = fieldn[ind4];
          f5 = fieldn[ind5];
          f6 = fieldn[ind6];

          pos = elOffset + i * 5;
          fieldf[pos] = K_CONST::ONE_THIRD * (f1 + f2 + f3);            // A1A2A3
          fieldf[pos + 1] = K_CONST::ONE_THIRD * (f4 + f5 + f6);        // A4A5A6
          fieldf[pos + 2] = K_CONST::ONE_FOURTH * (f1 + f2 + f5 + f4);  // 1254
          fieldf[pos + 3] = K_CONST::ONE_FOURTH * (f2 + f3 + f6 + f5);  // 2365
          fieldf[pos + 4] = K_CONST::ONE_FOURTH * (f3 + f1 + f4 + f6);  // 3146
        }
      }
      else if (strcmp(eltTypes[ic], "HEXA") == 0)
      {
        E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, pos;
        E_Float f1, f2, f3, f4, f5, f6, f7, f8;

        # pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cm[i*8] - 1;
          ind2 = cm[i*8 + 1] - 1;
          ind3 = cm[i*8 + 2] - 1;
          ind4 = cm[i*8 + 3] - 1;
          ind5 = cm[i*8 + 4] - 1;
          ind6 = cm[i*8 + 5] - 1;
          ind7 = cm[i*8 + 6] - 1;
          ind8 = cm[i*8 + 7] - 1;

          f1 = fieldn[ind1];
          f2 = fieldn[ind2];
          f3 = fieldn[ind3];
          f4 = fieldn[ind4];
          f5 = fieldn[ind5];
          f6 = fieldn[ind6];
          f7 = fieldn[ind7];
          f8 = fieldn[ind8];

          pos = elOffset + i * 6;
          fieldf[pos] = K_CONST::ONE_FOURTH * (f1 + f2 + f3 + f4);      // A1A2A3A4
          fieldf[pos + 1] = K_CONST::ONE_FOURTH * (f5 + f6 + f7 + f8);  // A5A6A7A8
          fieldf[pos + 2] = K_CONST::ONE_FOURTH * (f4 + f1 + f5 + f8);  // 4158
          fieldf[pos + 3] = K_CONST::ONE_FOURTH * (f2 + f3 + f7 + f6);  // A2A3A7A6
          fieldf[pos + 4] = K_CONST::ONE_FOURTH * (f1 + f2 + f6 + f5);  // A1A2A6A5
          fieldf[pos + 5] = K_CONST::ONE_FOURTH * (f3 + f4 + f8 + f7);  // A3A4A8A7
        }
      }
      else
      {
        fprintf(stderr, "Error: in K_POST::compUnstrNodes2Faces.\n");
        fprintf(stderr, "Unknown type of element.\n");
        exit(0);
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

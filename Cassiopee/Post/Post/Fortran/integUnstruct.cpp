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

// ============================================================================
// Compute surface integral of the field F, unstructured case
// Coordinates and field have the same size
//   I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3        TRI
//   I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C)+F(D))/4   QUAD
//   Aire(ABC) = ||AB^AC||/2
// ============================================================================
void K_POST::integUnstruct(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* ratio, const E_Float* surf, const E_Float* field,
  E_Float& result
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  
  result = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();

    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        E_Int ind1 = cm(i, 1) - 1;
        E_Int ind2 = cm(i, 2) - 1;
        E_Int ind3 = cm(i, 3) - 1;

        E_Float f1 = ratio[ind1] * field[ind1];
        E_Float f2 = ratio[ind2] * field[ind2];
        E_Float f3 = ratio[ind3] * field[ind3];

        result += K_CONST::ONE_THIRD * surf[i] * (f1 + f2 + f3);
      }
    }
    else if (strcmp(eltTypes[ic], "QUAD") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        E_Int ind1 = cm(i, 1) - 1;
        E_Int ind2 = cm(i, 2) - 1;
        E_Int ind3 = cm(i, 3) - 1;
        E_Int ind4 = cm(i, 4) - 1;

        E_Float f1 = ratio[ind1] * field[ind1];
        E_Float f2 = ratio[ind2] * field[ind2];
        E_Float f3 = ratio[ind3] * field[ind3];
        E_Float f4 = ratio[ind4] * field[ind4];

        result += K_CONST::ONE_FOURTH * surf[i] * (f1 + f2 + f3 + f4);
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute linear integral of field F, unstructured 1D case
//   I(AB) = Length(AB)*(F(A)+F(B))/2
// ============================================================================
void K_POST::integUnstruct1d(
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* ratio, const E_Float* length, const E_Float* field,
  E_Float& result
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  result = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    if (strcmp(eltTypes[ic], "BAR") != 0) continue;
  
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();

    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int ind1 = cm(i, 1) - 1;
      E_Int ind2 = cm(i, 2) - 1;

      E_Float f1 = ratio[ind1] * field[ind1];
      E_Float f2 = ratio[ind2] * field[ind2];

      result += length[i] * (f1 + f2);
    }
  }

  result *= K_CONST::ONE_HALF;

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute surface integral of field F, coordinates in nodes,
// field defined in centers, unstructured case
// ============================================================================
void K_POST::integUnstructNodeCenter(
  K_FLD::FldArrayI& cn,
  const E_Float* ratio, const E_Float* surf, const E_Float* field,
  E_Float& result
)
{
  E_Int nc = cn.getNConnect();
  result = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();

    for (E_Int i = 0; i < nelts; i++)
    {
      E_Float f = ratio[i] * field[i];
      result += surf[i] * f;
    }
  }
}
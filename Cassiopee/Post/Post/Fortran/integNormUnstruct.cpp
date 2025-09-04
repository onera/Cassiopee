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

// ============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
// and field have the same size
// I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3        - TRI
// I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)F(D))/4  - QUAD
// Aire(ABCD) = ||AB^AC||/2
// ============================================================================
void K_POST::integNormUnstruct(
  FldArrayI& cn, const char* eltType,
  const E_Float *ratio,
  const E_Float *sx, const E_Float *sy, const E_Float *sz, 
  const E_Float *field, E_Float *result
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int ind1, ind2, ind3;
    E_Float sum, f1, f2, f3;
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();

    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;

        f1 = ratio[ind1] * field[ind1];
        f2 = ratio[ind2] * field[ind2];
        f3 = ratio[ind3] * field[ind3];

        sum = K_CONST::ONE_THIRD * (f1 + f2 + f3);
        res1 += sx[i] * sum;
        res2 += sy[i] * sum;
        res3 += sz[i] * sum;
      }
    }
    else if (strcmp(eltTypes[ic], "QUAD") == 0)
    {
      E_Int ind4;
      E_Float f4;

      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;
        ind4 = cm(i, 4) - 1;

        f1 = ratio[ind1] * field[ind1];
        f2 = ratio[ind2] * field[ind2];
        f3 = ratio[ind3] * field[ind3];
        f4 = ratio[ind4] * field[ind4];

        sum = K_CONST::ONE_FOURTH * (f1 + f2 + f3 + f4);
        res1 += sx[i] * sum;
        res2 += sy[i] * sum;
        res3 += sz[i] * sum;
      }
    }
    else
    {
      fprintf(stderr, "Error: in K_POST::integNormUnstruct.\n");
      fprintf(stderr, "Unsupported type of element, %s.\n", eltTypes[ic]);
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute surface integral of the field F, coordinates are defined
// in nodes and F is defined in center, unstructured case
// ============================================================================
void K_POST::integNormUnstructNodeCenter(
  const E_Int nelts, const E_Float *ratio,
  const E_Float *nsurfx, const E_Float *nsurfy, const E_Float *nsurfz,
  const E_Float *field, E_Float *result
)
{
  E_Float f;
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  for (E_Int i = 0; i < nelts; i++)
  {
    f = ratio[i] * field[i];
    res1 += nsurfx[i] * f;
    res2 += nsurfy[i] * f;
    res3 += nsurfz[i] * f;
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

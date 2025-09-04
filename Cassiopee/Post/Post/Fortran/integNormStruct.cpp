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
//  Compute surface integral of the field F.vect(n), coordinates 
//      and field have the same size
//      I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//      Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
//=============================================================================
void K_POST::integNormStruct(
  const E_Int ni, const E_Int nj, const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* field, E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  E_Int ni1 = ni - 1;

  for (E_Int j = 0; j < nj - 1; j++)
  {
    for (E_Int i = 0; i < ni - 1; i++)
    {
      E_Int ind1 = i + j * ni;
      E_Int ind2 = ind1 + ni;
      E_Int ind3 = ind1 + 1;
      E_Int ind4 = ind3 + ni;
      E_Int ind = i + j * ni1;

      E_Float f1 = ratio[ind1] * field[ind1];
      E_Float f2 = ratio[ind2] * field[ind2];
      E_Float f3 = ratio[ind3] * field[ind3];
      E_Float f4 = ratio[ind4] * field[ind4];

      res1 += sx[ind] * (f1 + f2 + f3 + f4);
      res2 += sy[ind] * (f1 + f2 + f3 + f4);
      res3 += sz[ind] * (f1 + f2 + f3 + f4);
    }
  }

  result[0] = K_CONST::ONE_FOURTH * res1;
  result[1] = K_CONST::ONE_FOURTH * res2;
  result[2] = K_CONST::ONE_FOURTH * res3;
}

//=============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
void K_POST::integNormStructNodeCenter(
  const E_Int ni, const E_Int nj,
  const E_Float* ratio, const E_Float* sx, const E_Float* sy,
  const E_Float* sz, const E_Float* field, E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  for (E_Int j = 0; j < nj; j++)
  {
    for (E_Int i = 0; i < ni; i++)
    {
      E_Int ind = i + j * ni;
      E_Float ri = ratio[ind] * field[ind];

      res1 += ri * sx[ind];
      res2 += ri * sy[ind];
      res3 += ri * sz[ind];
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}
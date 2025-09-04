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
//     and field have the same size
//     I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//     Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
// ============================================================================
void K_POST::integNormProdStruct(
  const E_Int ni, const E_Int nj, const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result)
{
  E_Int ni1;
  E_Int ind1, ind2, ind3, ind4, ind;
  E_Float fx, fy, fz;
  E_Float r1, r2, r3, r4;

  ni1 = ni - 1;
  result = 0.0;

  for (E_Int j = 0; j < nj - 1; j++)
  {
    for (E_Int i = 0; i < ni - 1; i++)
    {
      ind1 = i + j * ni;
      ind2 = ind1 + ni;
      ind3 = ind1 + 1;
      ind4 = ind3 + ni;

      ind = i + j * ni1;

      r1 = ratio[ind1];
      r2 = ratio[ind2];
      r3 = ratio[ind3];
      r4 = ratio[ind4];

      fx = r1 * vx[ind1] + r2 * vx[ind2] + r3 * vx[ind3] + r4 * vx[ind4];
      fy = r1 * vy[ind1] + r2 * vy[ind2] + r3 * vy[ind3] + r4 * vy[ind4];
      fz = r1 * vz[ind1] + r2 * vz[ind2] + r3 * vz[ind3] + r4 * vz[ind4];

      result += sx[ind] * fx + sy[ind] * fy + sz[ind] * fz;
    }
  }

  result = K_CONST::ONE_FOURTH * result;
}

//==============================================================================
// Compute surface integral of the product vect(F).vect(n), coordinates 
// are defined in nodes and F is defined in nodes (center-based formulation)
//==============================================================================
void K_POST::integNormProdStructNodeCenter(
  const E_Int ni, const E_Int nj, const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result)
{
  E_Int ind;
  E_Float sum;

  result = 0.0;

  for (E_Int j = 0; j < nj; j++)
  {
    for (E_Int i = 0; i < ni; i++)
    {
      ind = i + j * ni;
      sum = sx[ind] * vx[ind] + sy[ind] * vy[ind] + sz[ind] * vz[ind];
      result += ratio[ind] * sum;
    }
  }
}
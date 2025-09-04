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
//     and field have the same size
//     I(ABC) = Aire(ABC) * (F(A) + F(B) + F(C)) / 3
//     Aire(ABC) = ||AB ^ AC|| / 2
// ============================================================================
void K_POST::integNormProdUnstruct(
  FldArrayI& cn, const char* eltType,
  const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  
  E_Int ind1, ind2, ind3;
  E_Float fx, fy, fz;
  E_Float r1, r2, r3;

  result = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    
    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;

        r1 = ratio[ind1];
        r2 = ratio[ind2];
        r3 = ratio[ind3];

        fx = r1 * vx[ind1] + r2 * vx[ind2] + r3 * vx[ind3];
        fy = r1 * vy[ind1] + r2 * vy[ind2] + r3 * vy[ind3];
        fz = r1 * vz[ind1] + r2 * vz[ind2] + r3 * vz[ind3];

        result += K_CONST::ONE_THIRD * (sx[i] * fx + sy[i] * fy + sz[i] * fz);
      }
    }
    else
    {
      fprintf(stderr, "Error: in K_POST::integNormProdUnstruct.\n");
      fprintf(stderr, "Unsupported type of element, %s.\n", eltTypes[ic]);
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute surface integral of the field F, coordinates are defined
// in nodes and F is defined in center, unstructured case
// ============================================================================
void K_POST::integNormProdUnstructNodeCenter(
  const E_Int nbt, const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result
)
{
  E_Float ri, sum;

  result = 0.0;
  for (E_Int i = 0; i < nbt; i++)
  {
    ri = ratio[i];
    sum = sx[i] * vx[i] + sy[i] * vy[i] + sz[i] * vz[i];
    result += ri * sum;
  }
}

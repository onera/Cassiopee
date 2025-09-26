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

//=============================================================================
// Integre les grandeurs de vect(F).vect(n)
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integNormProdUnstruct2D(E_Int center2node,
                                      E_Int posx, E_Int posy, E_Int posz,
                                      FldArrayI& cn, const char* eltType, FldArrayF& coord, 
                                      FldArrayF& F, FldArrayF& ratio, 
                                      E_Float& resultat)
{
  E_Float result = 0.;
  /*E_Int size = coord.getSize();*/

  E_Int ntotElts = 0;
  E_Int nc = cn.getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }
  
  // Compute surface of each "block" i cell, with coordinates coord
  FldArrayF nsurf(ntotElts, 3);
  K_METRIC::compNormUnstructSurf(
    cn, eltType,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    nsurf.begin(1), nsurf.begin(2), nsurf.begin(3)
  );

  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node and field F in center 
    integNormProdUnstructNodeCenter(
      ntotElts,
      ratio.begin(), 
      nsurf.begin(1), nsurf.begin(2), nsurf.begin(3),
      F.begin(1), F.begin(2), F.begin(3),
      result
    );
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    integNormProdUnstructCellCenter(
      cn, eltType,
      ratio.begin(),
      nsurf.begin(1), nsurf.begin(2), nsurf.begin(3),
      F.begin(1), F.begin(2), F.begin(3),
      result
    );
  }
  resultat += result;
  return 1;
}

// ============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
//     and field have the same size
//     I(ABC) = Aire(ABC) * (F(A) + F(B) + F(C)) / 3
//     Aire(ABC) = ||AB ^ AC|| / 2
// ============================================================================
void K_POST::integNormProdUnstructCellCenter(
  FldArrayI& cn, const char* eltType,
  const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  std::vector<E_Int> nepc(nc+1);
  nepc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
  }
  
  result = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int ind1, ind2, ind3, ind4;
    E_Float fx, fy, fz;
    E_Float r1, r2, r3, r4;
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int elOffset = nepc[ic];
    
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

        result += K_CONST::ONE_THIRD * (sx[i+elOffset] * fx + sy[i+elOffset] * fy + sz[i+elOffset] * fz);
      }
    }
    else if (strcmp(eltTypes[ic], "QUAD") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;
        ind4 = cm(i, 4) - 1;

        r1 = ratio[ind1];
        r2 = ratio[ind2];
        r3 = ratio[ind3];
        r4 = ratio[ind4];

        fx = r1 * vx[ind1] + r2 * vx[ind2] + r3 * vx[ind3] + r4 * vx[ind4];
        fy = r1 * vy[ind1] + r2 * vy[ind2] + r3 * vy[ind3] + r4 * vy[ind4];
        fz = r1 * vz[ind1] + r2 * vz[ind2] + r3 * vz[ind3] + r4 * vz[ind4];

        result += K_CONST::ONE_FOURTH * (sx[i+elOffset] * fx + sy[i+elOffset] * fy + sz[i+elOffset] * fz);
      }
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
  E_Float& result)
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

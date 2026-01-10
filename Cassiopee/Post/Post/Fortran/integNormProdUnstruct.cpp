/*    
    Copyright 2013-2026 ONERA.

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
    integNormProdUnstructCellCenter(
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
    integNormProdUnstructNodeCenter(
      cn,
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
void K_POST::integNormProdUnstructNodeCenter(
  FldArrayI& cn,
  const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result)
{
  E_Int nc = cn.getNConnect();

  std::vector<E_Int> nepc(nc+1);
  nepc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
  }
  
  result = 0.0;

  #pragma omp parallel
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int elOffset = nepc[ic];
    E_Int nvpe = cm.getNfld();
    E_Float nvpeinv = 1./nvpe;

    E_Int ind;
    E_Float fx, fy, fz;

    #pragma omp for reduction(+:result)
    for (E_Int i = 0; i < nelts; i++)
    {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
      for (E_Int j = 1; j <= nvpe; j++)
      {
        ind = cm(i, j) - 1;
        fx += ratio[ind] * vx[ind];
        fy += ratio[ind] * vy[ind];
        fz += ratio[ind] * vz[ind];
      }
      result += nvpeinv * (sx[i+elOffset] * fx + sy[i+elOffset] * fy + sz[i+elOffset] * fz);
    }
  }
}

// ============================================================================
// Compute surface integral of the field F, coordinates are defined
// in nodes and F is defined in center, unstructured case
// ============================================================================
void K_POST::integNormProdUnstructCellCenter(
  const E_Int nbt, const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result)
{
  result = 0.0;

  #pragma omp parallel for reduction(+:result)
  for (E_Int i = 0; i < nbt; i++)
  {
    E_Float ri = ratio[i];
    E_Float sum = sx[i] * vx[i] + sy[i] * vy[i] + sz[i] * vz[i];
    result += ri * sum;
  }
}

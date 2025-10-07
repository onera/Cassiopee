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
// Integre les grandeurs de F.vect(n)
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integNormUnstruct2D(E_Int center2node,
                             E_Int posx, E_Int posy, E_Int posz,
                             FldArrayI& cn, const char* eltType, FldArrayF& coord,
                             FldArrayF& F, FldArrayF& ratio,
                             FldArrayF& resultat)
{
  FldArrayF result(3);
  result.setAllValuesAtNull();

  E_Int nvars = F.getNfld();

  E_Float* res1 = resultat.begin(1);
  E_Float* res2 = resultat.begin(2);
  E_Float* res3 = resultat.begin(3);

  E_Int ntotElts = 0;
  E_Int nc = cn.getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }

  FldArrayF nsurf(ntotElts, 3);
  E_Float* nsurf1 = nsurf.begin(1);
  E_Float* nsurf2 = nsurf.begin(2);
  E_Float* nsurf3 = nsurf.begin(3);

  // Compute surface of each "block" i cell, with coordinates coord
  K_METRIC::compNormUnstructSurf(
    cn, eltType,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    nsurf1, nsurf2, nsurf3
  );

  if (center2node == 1)
  {
    for (E_Int n = 1; n <= nvars; n++)
    {
      // Compute integral, coordinates defined in node and field F in center
      integNormUnstructCellCenter(
        ntotElts, ratio.begin(),
        nsurf1, nsurf2, nsurf3, F.begin(n),
        result.begin());

      res1[n-1] += result[0];
      res2[n-1] += result[1];
      res3[n-1] += result[2];
    }
  }
  else
  {
    for (E_Int n = 1; n <= nvars; n++)
    {
      // Compute integral, coordinates and field have the same size
      integNormUnstructNodeCenter(
        cn,
        ratio.begin(),
        nsurf1, nsurf2, nsurf3, F.begin(n),
        result.begin());

      res1[n-1] += result[0];
      res2[n-1] += result[1];
      res3[n-1] += result[2];
    }
  }
  return 1;
}

// ============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
// and field have the same size
// I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3        - TRI
// I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)F(D))/4  - QUAD
// Aire(ABCD) = ||AB^AC||/2
// ============================================================================
void K_POST::integNormUnstructNodeCenter(
  FldArrayI& cn,
  const E_Float *ratio,
  const E_Float *sx, const E_Float *sy, const E_Float *sz, 
  const E_Float *field, E_Float *result)
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
  
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  #pragma omp parallel
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int elOffset = nepc[ic];
    E_Int nvpe = cm.getNfld();
    E_Float nvpeinv = 1./nvpe;

    E_Int ind;
    E_Float sum, f;

    #pragma omp for reduction(+:res1,res2,res3)
    for (E_Int i = 0; i < nelts; i++)
    {
      f = 0.0;
      for (E_Int j = 1; j <= nvpe; j++)
      {
        ind = cm(i, j) - 1;
        f += ratio[ind] * field[ind];
      }
      sum = nvpeinv * f;
      res1 += sx[i+elOffset] * sum;
      res2 += sy[i+elOffset] * sum;
      res3 += sz[i+elOffset] * sum;
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

// ============================================================================
// Compute surface integral of the field F, coordinates are defined
// in nodes and F is defined in center, unstructured case
// ============================================================================
void K_POST::integNormUnstructCellCenter(
  const E_Int nelts, const E_Float *ratio,
  const E_Float *nsurfx, const E_Float *nsurfy, const E_Float *nsurfz,
  const E_Float *field, E_Float *result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  #pragma omp parallel for reduction(+:res1,res2,res3)
  for (E_Int i = 0; i < nelts; i++)
  {
    E_Float f = ratio[i] * field[i];
    res1 += nsurfx[i] * f;
    res2 += nsurfy[i] * f;
    res3 += nsurfz[i] * f;
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

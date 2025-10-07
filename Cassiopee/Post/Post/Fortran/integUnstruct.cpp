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

//=============================================================================
// Integre les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
// Attention : cette routine n'integre que sur des elements triangulaires
//=============================================================================
E_Int K_POST::integUnstruct2D(E_Int center2node,
                              E_Int posx, E_Int posy, E_Int posz,
                              FldArrayI& cn, const char* eltType, FldArrayF& coord, 
                              FldArrayF& F, FldArrayF& ratio, 
                              FldArrayF& resultat)
{
  E_Float result = 0.;
  E_Int numberOfVariables = F.getNfld();
  E_Int ntotElts = 0;
  E_Int nc = cn.getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }
  FldArrayF surf(ntotElts);
  FldArrayF snx(ntotElts), sny(ntotElts), snz(ntotElts); // normale a la surface

  K_METRIC::compSurfUnstruct(
    cn, eltType,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    snx.begin(), sny.begin(), snz.begin(), surf.begin());

  if (center2node == 1) 
  {
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      // Compute integral, coordinates defined in node 
      // and field F in center 
      K_POST::integUnstructCellCenter(
        ntotElts,
        ratio.begin(), surf.begin(), F.begin(n),
        result
      );
      resultat[n-1] += result;
    }
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      K_POST::integUnstructNodeCenter(
        cn,
        ratio.begin(), surf.begin(), F.begin(n),
        result
      );
      resultat[n-1] += result;
    }
  }
  return 1;
}

//=============================================================================
// Integre les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
// Attention : cette routine n'integre que sur des elements "bar"
//=============================================================================
E_Int K_POST::integUnstruct1D(E_Int center2node,
                              E_Int posx, E_Int posy, E_Int posz,
                              FldArrayI& cn, const char* eltType, FldArrayF& coord, 
                              FldArrayF& F, FldArrayF& ratio, 
                              FldArrayF& resultat)
{
  E_Float result = 0.;
  E_Int numberOfVariables = F.getNfld();
  E_Int ntotElts = 0;
  E_Int nc = cn.getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }
  FldArrayF length(ntotElts);
  
  K_METRIC::compUnstructSurf1d(
    cn, eltType,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    length.begin());

  if (center2node == 1) 
  {
    for (E_Int n = 1 ; n <= numberOfVariables ; n++)
    {
      // Compute integral, coordinates defined in node 
      // and field F in center 
      K_POST::integUnstructCellCenter(
        ntotElts,
        ratio.begin(), length.begin(), F.begin(n),
        result
      );
      resultat[n-1] += result;
    }
  }
  else
  {
    for (E_Int n = 1 ; n <= numberOfVariables ; n++)    
    {
      // Compute integral, coordinates and field have the same size
      K_POST::integUnstructNodeCenter(
        cn,
        ratio.begin(), length.begin(), F.begin(n),
        result
      );
      resultat[n-1] += result;
    }
  }
  return 1;
}

// ============================================================================
// Compute surface integral of the field F, unstructured case
// Coordinates and field have the same size
//   I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3        TRI
//   I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C)+F(D))/4   QUAD
//   Aire(ABC) = ||AB^AC||/2
// ============================================================================
void K_POST::integUnstructNodeCenter(
  K_FLD::FldArrayI& cn,
  const E_Float* ratio, const E_Float* surf, const E_Float* field,
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
    E_Float f;

    #pragma omp for reduction(+:result)
    for (E_Int i = 0; i < nelts; i++)
    {
      f = 0.0;
      for (E_Int j = 1; j <= nvpe; j++)
      {
        ind = cm(i, j) - 1;
        f += ratio[ind] * field[ind];
      }
      result += nvpeinv * surf[i+elOffset] * f;
    }
  }
}

// ============================================================================
// Compute surface integral of field F, coordinates in nodes,
// field defined in centers, unstructured case
// ============================================================================
void K_POST::integUnstructCellCenter(
  const E_Int nelts, const E_Float* ratio, 
  const E_Float* surf, const E_Float* field,
  E_Float& result)
{
  result = 0.0;

  #pragma omp parallel for reduction(+:result)
  for (E_Int i = 0; i < nelts; i++)
  {
    E_Float f = ratio[i] * field[i];
    result += surf[i] * f;
  }
}
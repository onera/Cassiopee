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
# include "Array/Array.h"

//=============================================================================
// Integre "surfaciquement" les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
//=============================================================================
E_Int K_POST::integStruct2D(E_Int ni, E_Int nj, E_Int nk, 
                            E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                            FldArrayF& coord, FldArrayF& F, 
                            FldArrayF& ratio, FldArrayF& resultat)
{
  E_Int NI, NJ;
  E_Float result = 0.;
  E_Int numberOfVariables = F.getNfld();

  if      (nk == 1) {NI = ni; NJ = nj;}
  else if (nj == 1) {NI = ni; NJ = nk;}
  else if (ni == 1) {NI = nj; NJ = nk;}
  else return 0;

  E_Int ncells = (NI-1)*(NJ-1);
  
  // Compute surface of each "block" i cell, with coordinates coord
  FldArrayF surf(ncells);

  K_METRIC::compSurfStruct2D(
    NI, NJ, 1,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    surf.begin());
  
  if (center2node == 1) 
  {
    // Compute integral, coordinates defined in node 
    // and field F in center 
    for (E_Int n = 1; n <= numberOfVariables; n++)
    { 
      K_POST::integStructCellCenter2D(
        NI-1, NJ-1,
        ratio.begin(), surf.begin(), F.begin(n),
        result);
      
      resultat[n-1] += result;
    }
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      K_POST::integStructNodeCenter2D(
        NI, NJ,
        ratio.begin(), surf.begin(), F.begin(n),
        result);
      
      resultat[n-1] += result;
    }
  } 
  return 1;
}

//=============================================================================
// Integre "surfaciquement" les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
//=============================================================================
E_Int K_POST::integStruct1D(E_Int ni, E_Int nj, E_Int nk, 
                            E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                            FldArrayF& coord, FldArrayF& F, 
                            FldArrayF& ratio, FldArrayF& resultat)
{
  E_Int NI;
  E_Float result = 0.;
  E_Int numberOfVariables = F.getNfld();

  if      (ni > 1) NI = ni;
  else if (nj > 1) NI = nj;
  else if (nk > 1) NI = nk;
  else return 0;

  E_Int ncells = (NI-1);

  // Compute surface of each "block" i cell, with coordinates coord
  FldArrayF length(ncells);
  K_METRIC::compSurfStruct1D(
    NI, 1, 1,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    length.begin());
  
  if (center2node == 1) 
  {
    // Compute integral, coordinates defined in node 
    // and field F in center 
    for (E_Int n = 1; n <= numberOfVariables; n++)
    { 
      K_POST::integStructCellCenter1D(
        NI-1,
        ratio.begin(), length.begin(), F.begin(n),
        result);

      resultat[n-1] += result;
    }
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      K_POST::integStructNodeCenter1D(
        NI,
        ratio.begin(), length.begin(), F.begin(n),
        result);
      
      resultat[n-1] += result;
    }
  } 
  return 1;
}

// ============================================================================
// Compute surface integral of the field F, coordinates 
// and field have the same size
//   I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//   Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
// ============================================================================
void K_POST::integStructNodeCenter2D(
  const E_Int ni, const E_Int nj,
  const E_Float* ratio, const E_Float* surf, const E_Float* field,
  E_Float& result)
{
  E_Int ni1 = ni - 1;
  result = 0.0;

  #pragma omp parallel
  {
    E_Int ind, ind1, ind2, ind3, ind4;
    E_Float f1, f2, f3, f4;

    #pragma omp for collapse(2) reduction(+:result)
    for (E_Int j = 0; j <= nj - 2; j++)
    {
      for (E_Int i = 0; i <= ni - 2; i++)
      {
        ind1 = i + j * ni;
        ind2 = ind1 + ni;
        ind3 = ind1 + 1;
        ind4 = ind3 + ni;
        ind = i + j * ni1;

        f1 = ratio[ind1] * field[ind1];
        f2 = ratio[ind2] * field[ind2];
        f3 = ratio[ind3] * field[ind3];
        f4 = ratio[ind4] * field[ind4];

        result += surf[ind] * (f1 + f2 + f3 + f4);
      }
    }
  }
  result *= 0.25;
}

// ============================================================================
// Compute linear integral of field F, structured 1D case
//   I(AB) = Length(AB)*(F(A)+F(B))/2
// ============================================================================
void K_POST::integStructNodeCenter1D(
  const E_Int ni,
  const E_Float* ratio, const E_Float* length, const E_Float* field,
  E_Float& result)
{
  result = 0.0;

  E_Int ind, ind1, ind2;
  E_Float f1, f2;

  for (E_Int i = 0; i <= ni - 2; i++)
  {
    ind1 = i;
    ind2 = i + 1;
    ind = i;

    f1 = ratio[ind1] * field[ind1];
    f2 = ratio[ind2] * field[ind2];

    result += length[ind] * (f1 + f2);
  }
  result *= 0.5;
}

// ============================================================================
// Compute surface integral of field F, coordinates in nodes,
// field defined in centers, structured case
// IN: ni1, nj1 : dim en centres
// ============================================================================
void K_POST::integStructCellCenter2D(
  const E_Int ni1, const E_Int nj1,
  const E_Float* ratio, const E_Float* surf, const E_Float* field,
  E_Float& result)
{
  result = 0.0;

  #pragma omp parallel
  {
    E_Int ind;

    #pragma omp for collapse(2) reduction(+:result)
    for (E_Int j = 0; j <= nj1 - 1; j++)
    {
      for (E_Int i = 0; i <= ni1 - 1; i++)
      {
        ind = i + j * ni1;
        result += ratio[ind] * surf[ind] * field[ind];
      }
    }
  }
}

// ============================================================================
// Compute linear integral of field F, node/center case, 1D
// ============================================================================
void K_POST::integStructCellCenter1D(
  const E_Int ni,
  const E_Float* ratio, const E_Float* length, const E_Float* field,
  E_Float& result)
{
  result = 0.0;
  for (E_Int i = 0; i <= ni - 1; i++)
  {
    result += ratio[i] * length[i] * field[i];
  }
}

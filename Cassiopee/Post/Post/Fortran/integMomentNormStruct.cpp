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
// Integre les grandeurs de M = OM^F.vect(n)
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integMomentNormStruct2D(E_Int ni, E_Int nj, E_Int nk, 
                                      E_Int center2node, 
                                      E_Int posx, E_Int posy, E_Int posz,
                                      E_Float cx, E_Float cy, E_Float cz, 
                                      FldArrayF& coord, 
                                      FldArrayF& F, FldArrayF& ratio, 
                                      FldArrayF& resultat)
{
  E_Int NI, NJ;
  FldArrayF result(3);
  result.setAllValuesAtNull();

  E_Int numberOfVariables = F.getNfld();
  E_Float* resultat1 = resultat.begin(1);
  E_Float* resultat2 = resultat.begin(2);
  E_Float* resultat3 = resultat.begin(3);

  if      (nk == 1) {NI = ni; NJ = nj;}
  else if (nj == 1) {NI = ni; NJ = nk;}
  else if (ni == 1) {NI = nj; NJ = nk;}
  else return 0;
 
  // Compute surface of each "block" i cell, with coordinates coord
  E_Int ncells = (NI-1) * (NJ-1); 
  FldArrayF nsurf(ncells, 3);  

  K_METRIC::compNormStructSurf(
    NI, NJ, coord.begin(posx), coord.begin(posy), coord.begin(posz), 
    nsurf.begin(1), nsurf.begin(2), nsurf.begin(3));

  if (center2node == 1) 
  { 
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {  
      // Compute integral, coordinates defined in node and field F in center 
      integMomentNormStructCellCenter2D(
        NI, NJ, cx, cy, cz, ratio.begin(), 
        coord.begin(posx), coord.begin(posy),coord.begin(posz),
        nsurf.begin(1),nsurf.begin(2), nsurf.begin(3),   
        F.begin(), result.begin()
      );
      
      resultat1[n-1] += result[0];   
      resultat2[n-1] += result[1];
      resultat3[n-1] += result[2];
    }
  }
  else
  {
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      // Compute integral, coordinates and field have the same size
      integMomentNormStructNodeCenter2D(
        NI, NJ, cx, cy, cz, ratio.begin(), 
        coord.begin(posx), coord.begin(posy), coord.begin(posz),
        nsurf.begin(1), nsurf.begin(2), nsurf.begin(3), F.begin(),  
        result.begin()
      );

      resultat1[n-1] += result[0];   
      resultat2[n-1] += result[1];
      resultat3[n-1] += result[2];
    }
  }
  return 1;
}

// ============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
//       and F have the same size
// ============================================================================
void K_POST::integMomentNormStructNodeCenter2D(
  const E_Int ni, const E_Int nj, const E_Float cx,
  const E_Float cy, const E_Float cz, const E_Float* ratio, const E_Float* xt,
  const E_Float* yt, const E_Float* zt, const E_Float* sx, const E_Float* sy,
  const E_Float* sz, const E_Float* field, E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  E_Int ni1 = ni - 1;

  #pragma omp parallel
  {
    E_Int ind, ind1, ind2, ind3, ind4;
    E_Float f, f1, f2, f3, f4;
    E_Float mx, my, mz;
    E_Float sx0, sy0, sz0;
    E_Float centerx, centery, centerz;

    #pragma omp for collapse(2) reduction(+:res1,res2,res3)
    for (E_Int j = 0; j < nj - 1; j++)
    {
      for (E_Int i = 0; i < ni - 1; i++)
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

        f = 0.25 * (f1 + f2 + f3 + f4);

        sx0 = sx[ind];
        sy0 = sy[ind];
        sz0 = sz[ind];

        centerx = 0.25 * (xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4]);
        centerx = centerx - cx;

        centery = 0.25 * (yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4]);
        centery = centery - cy;

        centerz = 0.25 * (zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4]);
        centerz = centerz - cz;

        mx = centery * sz0 - centerz * sy0;
        my = centerz * sx0 - centerx * sz0;
        mz = centerx * sy0 - centery * sx0;

        res1 += f * mx;
        res2 += f * my;
        res3 += f * mz;
      }
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

//=============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
//       are defined in nodes and F is defined in center
//=============================================================================
void K_POST::integMomentNormStructCellCenter2D(
  const E_Int ni, const E_Int nj,
  const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* sx,
  const E_Float* sy, const E_Float* sz, const E_Float* field, E_Float* result)
{
  E_Float res1 = 0.0;
  E_Float res2 = 0.0;
  E_Float res3 = 0.0;

  E_Int ni1 = ni - 1;

  #pragma omp parallel
  {
    E_Int ind, ind1, ind2, ind3, ind4;
    E_Float f;
    E_Float mx, my, mz;
    E_Float sx0, sy0, sz0;
    E_Float centerx, centery, centerz;

    #pragma omp for collapse(2) reduction(+:res1,res2,res3)
    for (E_Int j = 0; j < nj - 1; j++)
    {
      for (E_Int i = 0; i < ni - 1; i++)
      {
        ind1 = i + j * ni;
        ind2 = ind1 + ni;
        ind3 = ind1 + 1;
        ind4 = ind3 + ni;
        ind = i + j * ni1;

        f = ratio[ind] * field[ind];

        sx0 = sx[ind];
        sy0 = sy[ind];
        sz0 = sz[ind];

        centerx = xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4];
        centerx = 0.25 * centerx - cx;

        centery = yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4];
        centery = 0.25 * centery - cy;

        centerz = zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4];
        centerz = 0.25 * centerz - cz;

        mx = centery * sz0 - centerz * sy0;
        my = centerz * sx0 - centerx * sz0;
        mz = centerx * sy0 - centery * sx0;

        res1 += f * mx;
        res2 += f * my;
        res3 += f * mz;
      }
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}
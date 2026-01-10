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
// Integre les grandeurs de M = OM^F.vect(n)
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integMomentNormUnstruct2D(E_Int center2node,
                                      E_Int posx, E_Int posy, E_Int posz,
                                      E_Float cx, E_Float cy, E_Float cz, 
                                      FldArrayI& cn, const char* eltType, FldArrayF& coord, 
                                      FldArrayF& F, FldArrayF& ratio, 
                                      FldArrayF& resultat)
{
  FldArrayF result(3);
  result.setAllValuesAtNull();

  E_Int numberOfVariables = F.getNfld();
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
    nsurf1, nsurf2, nsurf3);

  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node 
    // and field F in center 
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      integMomentNormUnstructCellCenter(
        cn,
        cx, cy, cz, ratio.begin(), 
        coord.begin(posx), coord.begin(posy), coord.begin(posz),
        nsurf1, nsurf2, nsurf3, F.begin(n),
        result.begin()
      );

      res1[n-1] += result[0];   
      res2[n-1] += result[1];
      res3[n-1] += result[2];
    }
  }
  else
  {
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      // Compute integral, coordinates and field have the same size
      integMomentNormUnstructNodeCenter(
        cn,
        cx, cy, cz, ratio.begin(), 
        coord.begin(posx), coord.begin(posy), coord.begin(posz), 
        nsurf1, nsurf2, nsurf3, F.begin(n), 
        result.begin()
      );
    
      res1[n-1] += result[0];   
      res2[n-1] += result[1];
      res3[n-1] += result[2];        
    }
  }
  return 1;
}


// ============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
// and F have the same size
// ============================================================================
void K_POST::integMomentNormUnstructNodeCenter(
  FldArrayI& cn,
  const E_Float cx, const E_Float cy, const E_Float cz, 
  const E_Float* ratio, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* sx, const E_Float* sy, const E_Float* sz, 
  const E_Float* field, E_Float* result)
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
    E_Int ind;
    E_Float f, fi;
    E_Float mx, my, mz, sx0, sy0, sz0;
    E_Float centerx, centery, centerz;
  
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int elOffset = nepc[ic];
    E_Int nvpe = cm.getNfld();
    E_Float nvpeinv = 1./nvpe;

    #pragma omp for reduction(+:res1,res2,res3)
    for (E_Int i = 0; i < nelts; i++)
    {
      f = 0.0;
      centerx = 0.0;
      centery = 0.0;
      centerz = 0.0;

      for (E_Int j = 1; j <= nvpe; j++)
      {
        ind = cm(i, j) - 1;

        fi = ratio[ind] * field[ind];
        f += fi;

        centerx += xt[ind];
        centery += yt[ind];
        centerz += zt[ind];
      }

      f *= nvpeinv;
      centerx = nvpeinv * centerx - cx;
      centery = nvpeinv * centery - cy;
      centerz = nvpeinv * centerz - cz;

      sx0 = sx[i+elOffset];
      sy0 = sy[i+elOffset];
      sz0 = sz[i+elOffset];

      mx = centery * sz0 - centerz * sy0;
      my = centerz * sx0 - centerx * sz0;
      mz = centerx * sy0 - centery * sx0;

      res1 += f * mx;
      res2 += f * my;
      res3 += f * mz;
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

// ============================================================================
// Compute linear integral of the moment.norm (OM^F.n), coordinates 
// are defined in nodes and F is defined in center, unstructured case
// ============================================================================
void K_POST::integMomentNormUnstructCellCenter(
  FldArrayI& cn,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* sx, const E_Float* sy, const E_Float* sz, 
  const E_Float* field, E_Float* result)
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
    E_Int ind;
    E_Float f;
    E_Float mx, my, mz, sx0, sy0, sz0;
    E_Float centerx, centery, centerz;
  
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int elOffset = nepc[ic];
    E_Int nvpe = cm.getNfld();
    E_Float nvpeinv = 1./nvpe;

    #pragma omp for reduction(+:res1,res2,res3)
    for (E_Int i = 0; i < nelts; i++)
    {
      f = ratio[i+elOffset] * field[i+elOffset];
      centerx = 0.0;
      centery = 0.0;
      centerz = 0.0;

      for (E_Int j = 1; j <= nvpe; j++)
      {
        ind = cm(i, j) - 1;

        centerx += xt[ind];
        centery += yt[ind];
        centerz += zt[ind];
      }

      centerx = nvpeinv * centerx - cx;
      centery = nvpeinv * centery - cy;
      centerz = nvpeinv * centerz - cz;

      sx0 = sx[i+elOffset];
      sy0 = sy[i+elOffset];
      sz0 = sz[i+elOffset];

      mx = centery * sz0 - centerz * sy0;
      my = centerz * sx0 - centerx * sz0;
      mz = centerx * sy0 - centery * sx0;

      res1 += f * mx;
      res2 += f * my;
      res3 += f * mz;
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}

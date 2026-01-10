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
// Integre les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integMomentUnstruct2D(E_Int center2node,
                                    E_Int posx, E_Int posy, E_Int posz,
                                    E_Float cx, E_Float cy, E_Float cz, 
                                    FldArrayI& cn, const char* eltType, FldArrayF& coord, 
                                    FldArrayF& F, FldArrayF& ratio, 
                                    FldArrayF& resultat)
{
  FldArrayF res(3);
  E_Int ntotElts = 0;
  E_Int nc = cn.getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }
 
  // Compute surface of each "block" i cell, with coordinates coord
  FldArrayF snx(ntotElts); // normale a la surface  
  FldArrayF sny(ntotElts);
  FldArrayF snz(ntotElts);
  FldArrayF surf(ntotElts);

  K_METRIC::compSurfUnstruct(
    cn, eltType,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    snx.begin(), sny.begin(), snz.begin(), surf.begin()
  );

  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node and field F in center 
    integMomentUnstructCellCenter(
      cn,
      cx, cy, cz, ratio.begin(),
      coord.begin(posx), coord.begin(posy), coord.begin(posz),
      surf.begin(), F.begin(1), F.begin(2), F.begin(3),
      res.begin()
    );
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    integMomentUnstructNodeCenter(
      cn,
      cx, cy, cz, ratio.begin(),
      coord.begin(posx), coord.begin(posy), coord.begin(posz),
      surf.begin(), F.begin(1), F.begin(2), F.begin(3),
      res.begin()
    );
  }

  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];
   
  return 1;
}

//=============================================================================
// Integre les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integMomentUnstruct1D(E_Int center2node,
                                    E_Int posx, E_Int posy, E_Int posz,
                                    E_Float cx, E_Float cy, E_Float cz, 
                                    FldArrayI& cn, const char* eltType, FldArrayF& coord, 
                                    FldArrayF& F, FldArrayF& ratio, 
                                    FldArrayF& resultat)
{
  FldArrayF res(3);
  E_Int ntotElts = 0;
  E_Int nc = cn.getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }

  // Compute surface of each "block" i cell, with coordinates coord
  FldArrayF length(ntotElts);
  K_METRIC::compUnstructSurf1d(
    cn, eltType,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    length.begin()
  );

  if (center2node == 1) 
  { 
    // Compute integral, coordinates defined in node and field F in center 
    integMomentUnstructCellCenter(
      cn,
      cx, cy, cz, ratio.begin(),
      coord.begin(posx), coord.begin(posy), coord.begin(posz), 
      length.begin(), F.begin(1), F.begin(2), F.begin(3),
      res.begin()
    );
  }
  else
  {
    integMomentUnstructNodeCenter(
      cn,
      cx, cy, cz, ratio.begin(),
      coord.begin(posx), coord.begin(posy), coord.begin(posz), 
      length.begin(), F.begin(1), F.begin(2), F.begin(3),
      res.begin()
    );
  }

  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];   

  return 1;
}

// ============================================================================
// Compute surface integral of the moment M (OM^F), coordinates 
// and field have the same size
// I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
// Aire(ABCD) = ||AB^AC||/2 + ||DB^DC||/2
// ============================================================================
void K_POST::integMomentUnstructNodeCenter(
  FldArrayI& cn,
  const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz, E_Float* result)
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
    E_Float fx, fy, fz;
    E_Float f1, f2, f3;
    E_Float mx, my, mz;
    E_Float dx, dy, dz, ri, si;

    #pragma omp for reduction(+:res1,res2,res3)
    for (E_Int i = 0; i < nelts; i++)
    {
      f1 = 0.0;
      f2 = 0.0;
      f3 = 0.0;

      for (E_Int j = 1; j <= nvpe; j++)
      {
        ind = cm(i, j) - 1;
        
        ri = ratio[ind];
        dx = xt[ind] - cx;
        dy = yt[ind] - cy;
        dz = zt[ind] - cz;
        fx = vx[ind];
        fy = vy[ind];
        fz = vz[ind];

        mx = dy * fz - dz * fy;
        my = dz * fx - dx * fz;
        mz = dx * fy - dy * fx;

        f1 += ri*mx;
        f2 += ri*my;
        f3 += ri*mz;
      }

      si = surf[i+elOffset];
      res1 += nvpeinv*si*f1;
      res2 += nvpeinv*si*f2;
      res3 += nvpeinv*si*f3;
    }
  }

  result[0] = res1; 
  result[1] = res2;
  result[2] = res3;
}

// ============================================================================
// Compute surface integral of the moment M (OM^F)
// coordinates are defined in nodes and F is defined in center (unstructured)
// ============================================================================
void K_POST::integMomentUnstructCellCenter(
  FldArrayI& cn,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result)
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
    E_Float mx, my, mz, sri;
    E_Float centerx, centery, centerz;

    #pragma omp for reduction(+:res1,res2,res3)
    for (E_Int i = 0; i < nelts; i++)
    {
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

      mx = centery * vz[i] - centerz * vy[i];
      my = centerz * vx[i] - centerx * vz[i];
      mz = centerx * vy[i] - centery * vx[i];

      sri = surf[i+elOffset] * ratio[i];
      res1 += sri * mx;
      res2 += sri * my;
      res3 += sri * mz;
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}
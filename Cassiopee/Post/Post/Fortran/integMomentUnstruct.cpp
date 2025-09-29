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
    integMomentUnstructNodeCenter(
      cn, eltType,
      cx, cy, cz, ratio.begin(),
      coord.begin(posx), coord.begin(posy), coord.begin(posz),
      surf.begin(), F.begin(1), F.begin(2), F.begin(3),
      res.begin()
    );
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    integMomentUnstructCellCenter(
      cn, eltType,
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
    integMomentUnstructNodeCenter(
      cn, eltType,
      cx, cy, cz, ratio.begin(),
      coord.begin(posx), coord.begin(posy), coord.begin(posz), 
      length.begin(), F.begin(1), F.begin(2), F.begin(3),
      res.begin()
    );
  }
  else
  {
    integMomentUnstructCellCenter(
      cn, eltType,
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
void K_POST::integMomentUnstructCellCenter(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz, E_Float* result)
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

  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int elOffset = nepc[ic];
    E_Int nfpe=1;
    E_Float nfpeinv;
    
    E_Int ind;
    E_Float fx, fy, fz;
    E_Float f1, f2, f3;
    E_Float mx, my, mz;
    E_Float dx, dy, dz, ri, si;
    E_Float res1 = 0.0;
    E_Float res2 = 0.0;
    E_Float res3 = 0.0;

    if (strcmp(eltTypes[ic], "BAR") == 0) nfpe = 2;
    else if (strcmp(eltTypes[ic], "TRI") == 0) nfpe = 3;
    else if (strcmp(eltTypes[ic], "QUAD") == 0) nfpe = 4;

    nfpeinv = 1./nfpe;

    for (E_Int i = 0; i < nelts; i++)
    {
      f1 = 0.0;
      f2 = 0.0;
      f3 = 0.0;

      for (E_Int j = 1; j <= nfpe; j++)
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
      res1 += si*f1;
      res2 += si*f2;
      res3 += si*f3;
    }

    result[0] += nfpeinv*res1; 
    result[1] += nfpeinv*res2;
    result[2] += nfpeinv*res3;
  }
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute surface integral of the moment M (OM^F)
// coordinates are defined in nodes and F is defined in center (unstructured)
// ============================================================================
void K_POST::integMomentUnstructNodeCenter(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result)
{
  E_Float res1, res2, res3;
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

  result[0] = 0.0;
  result[1] = 0.0;
  result[2] = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int ind1, ind2, ind3, ind4;
    E_Float mx, my, mz, sri;
    E_Float centerx, centery, centerz;
  
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

        centerx = K_CONST::ONE_THIRD * (xt[ind1] + xt[ind2] + xt[ind3]) - cx;
        centery = K_CONST::ONE_THIRD * (yt[ind1] + yt[ind2] + yt[ind3]) - cy;
        centerz = K_CONST::ONE_THIRD * (zt[ind1] + zt[ind2] + zt[ind3]) - cz;

        mx = centery * vz[i] - centerz * vy[i];
        my = centerz * vx[i] - centerx * vz[i];
        mz = centerx * vy[i] - centery * vx[i];

        sri = surf[i+elOffset] * ratio[i];
        res1 = sri * mx;
        res2 = sri * my;
        res3 = sri * mz;

        result[0] += res1;
        result[1] += res2;
        result[2] += res3;
      
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

        centerx = K_CONST::ONE_FOURTH * (xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4]) - cx;
        centery = K_CONST::ONE_FOURTH * (yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4]) - cy;
        centerz = K_CONST::ONE_FOURTH * (zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4]) - cz;

        mx = centery * vz[i] - centerz * vy[i];
        my = centerz * vx[i] - centerx * vz[i];
        mz = centerx * vy[i] - centery * vx[i];

        sri = surf[i+elOffset] * ratio[i];
        res1 = sri * mx;
        res2 = sri * my;
        res3 = sri * mz;

        result[0] += res1;
        result[1] += res2;
        result[2] += res3;
      
      }
    }
    else if (strcmp(eltTypes[ic], "BAR") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;

        centerx = K_CONST::ONE_HALF * (xt[ind1] + xt[ind2]) - cx;
        centery = K_CONST::ONE_HALF * (yt[ind1] + yt[ind2]) - cy;
        centerz = K_CONST::ONE_HALF * (zt[ind1] + zt[ind2]) - cz;

        mx = centery * vz[i] - centerz * vy[i];
        my = centerz * vx[i] - centerx * vz[i];
        mz = centerx * vy[i] - centery * vx[i];

        sri = surf[i+elOffset] * ratio[i];
        res1 = sri * mx;
        res2 = sri * my;
        res3 = sri * mz;

        result[0] += res1;
        result[1] += res2;
        result[2] += res3;
      }
    }
  }
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}
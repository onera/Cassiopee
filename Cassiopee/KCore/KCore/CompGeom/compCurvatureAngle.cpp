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

# include "CompGeom/compGeom.h"
# include "Connect/connect.h"
# include <math.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   Calcul de l'angle entre les segments pour un "i-array" dans le plan (x,y)
   IN: npts: nbre de pts du i-array
   IN: xt, yt, zt: coord. des pts du i-array
   IN: dirVect is (approximatively) the direction vector orthogonal to the 
   mean plane in which lies the i-array
   OUT: angle: angle entre 0 et 360 degres.
*/
//=============================================================================
void K_COMPGEOM::compCurvatureAngle(E_Int npts, 
                                    E_Float* xt, E_Float* yt, E_Float* zt, 
                                    E_Float* dirVect,
                                    FldArrayF& angle)
{
  E_Int indm, indp, ind;
  E_Float x, xp, xm;
  E_Float y, yp, ym;
  E_Float z, zp, zm;
  E_Float dx1, dx2, dy1, dy2, dz1, dz2;
  E_Float n1, n2;
  E_Float inv, pv1, pv2, pv3, pvn;
  E_Float pv, ps;
  E_Float pi = 4*atan(1.);
  E_Float alp = 180./pi;

  FldArrayI node(npts+4);
  E_Int* nodep = node.begin();
 
  // traitement des bords : contour ferme ou ouvert ?
  E_Int N = npts-1;
  E_Float dx = xt[0]-xt[N];
  E_Float dy = yt[0]-yt[N];
  E_Float dz = zt[0]-zt[N];
  E_Float dist = dx*dx+dy*dy+dz*dz;
  if (dist < 1.e-12)
  {
    for (E_Int i = 0; i <= N; i++) nodep[i+2] = i;
    nodep[0] = N-2; nodep[1] = N-1;
    nodep[N+3] = 1; nodep[N+4] = 2;
  }
  else  // ouvert 
  {
    for (E_Int i = 0; i <= N; i++) nodep[i+2] = i;
    nodep[0] = 0; nodep[1] = 0;
    nodep[N+3] = N; nodep[N+4] = N;
  }
  for (E_Int i = 2; i <= N+2; i++)
  {
    indm = nodep[i-1]; indp = nodep[i+1]; ind = i-2;
    x = xt[ind]; xp = xt[indp]; xm = xt[indm];
    y = yt[ind]; yp = yt[indp]; ym = yt[indm];
    z = zt[ind]; zp = zt[indp]; zm = zt[indm];
    dx1 = xm-x; dx2 = xp-x;
    dy1 = ym-y; dy2 = yp-y;
    dz1 = zm-z; dz2 = zp-z;
    n1 = dx1*dx1+dy1*dy1+dz1*dz1;
    n2 = dx2*dx2+dy2*dy2+dz2*dz2;
    n1 = K_FUNC::E_max(n1, 1.e-12);
    n2 = K_FUNC::E_max(n2, 1.e-12);
    
    inv = 1./sqrt(n1*n2);
    pv1 = dy1*dz2 - dz1*dy2;
    pv2 = dz1*dx2 - dx1*dz2;
    pv3 = dx1*dy2 - dy1*dx2;
    
    //E_Float normdir = K_FUNC::normalize<3>(dirVect);
    pv = pv1*dirVect[0] + pv2*dirVect[1] + pv3*dirVect[2];
    pvn = pv1*pv1+pv2*pv2+pv3*pv3;
    pvn = sqrt(pvn) * inv;
    
    //printf("%f %f\n", pv, pvn);
    pv = pv * inv;
    ps = dx2*dx1 + dy2*dy1 + dz1*dz2;
    ps = ps * inv;
    
    if (pv >= 0. && ps >= 0.)
      angle[ind] = asin(pvn)* alp;
    else if (pv <= 0. && ps >= 0.)
      angle[ind] = 360. - asin(pvn) * alp;
    else if (pv >= 0. && ps <= 0.)
      angle[ind] = 180. - asin(pvn) * alp;
    else
      angle[ind] = 180. + asin(pvn) * alp;
  }
}

//=============================================================================
/* Calcul de l'angle entre les segments pour une "BAR" dans le plan (x,y)
   IN: npts: nbre de pts de la BAR
   IN: xt, yt, zt: coord. des pts de la BAR
   IN: cn: connectivite vertex-elts
   IN: dirVect is (approximatively) the direction vector orthogonal to the 
   mean plane in which lies the i-array
   OUT: angle entre 0 et 360 degres.
*/
//=============================================================================
void K_COMPGEOM::compCurvatureAngleForBar(
  E_Int npts, 
  E_Float* xt, E_Float* yt, E_Float* zt,
  FldArrayI& cn,  E_Float* dirVect,
  FldArrayF& angle)
{
  E_Int ind;
  E_Int i;
  E_Int indp, indm, s;
  E_Int e1, e2, ind1, ind2, ind3, ind4;
  E_Float x, xp, xm;
  E_Float y, yp, ym;
  E_Float z, zp, zm;
  E_Float dx1, dx2, dy1, dy2, dz1, dz2;
  E_Float n1, n2;
  E_Float inv, pv1, pv2, pv3;
  E_Float pv, ps;
  E_Float pi = 4*atan(1.);
  E_Float alp = 180./pi;
  
  vector< vector<E_Int> > cVE(npts);
  K_CONNECT::connectEV2VE(cn, cVE);
  for (i = 0; i < npts; i++)
  {
    ind = i;
    s = cVE[i].size();
    if (s == 2) // 2 elts
    {
      e1 = cVE[i][0];
      e2 = cVE[i][1];
      ind1 = cn(e1,1)-1;
      ind2 = cn(e1,2)-1;
      ind3 = cn(e2,1)-1;
      ind4 = cn(e2,2)-1;
      if (ind == ind1 && ind == ind4)
      {
        indp = ind2;
        indm = ind3;
      }
      else if (ind == ind2 && ind == ind3)
      {
        indp = ind4;
        indm = ind1;
      }
      else if (ind == ind1 && ind == ind3)
      {
        // reverse case : ambiguous
        indp = ind2;
        indm = ind4;
      }
      else
      {
        // reverse case: ambiguous
        indp = ind1;
        indm = ind3;
      }           
    }
    else if (s == 1)
    {
      e1 = cVE[i][0];
      ind1 = cn(e1,1)-1;
      ind2 = cn(e1,2)-1;
      indp = ind1;
      indm = ind2;
    }
    else
    {
      indp = ind;
      indm = ind;
    }
    
    // Compute angle
    x = xt[ind]; xp = xt[indp]; xm = xt[indm];
    y = yt[ind]; yp = yt[indp]; ym = yt[indm];
    z = zt[ind]; zp = zt[indp]; zm = zt[indm];
    dx1 = xm-x; dx2 = xp-x;
    dy1 = ym-y; dy2 = yp-y;
    dz1 = zm-z; dz2 = zp-z;
    n1 = dx1*dx1+dy1*dy1+dz1*dz1;
    n2 = dx2*dx2+dy2*dy2+dz2*dz2;
    n1 = K_FUNC::E_max(n1, 1.e-12);
    n2 = K_FUNC::E_max(n2, 1.e-12);
    inv = 1./sqrt(n1*n2);
    pv1 = dy1*dz2 - dz1*dy2;
    pv2 = dz1*dx2 - dx1*dz2;
    pv3 = dx1*dy2 - dy1*dx2;
    //E_Float normdir = K_FUNC::normalize<3>(dirVect);
    pv = pv1*dirVect[0] + pv2*dirVect[1] + pv3*dirVect[2]; 
    pv = pv * inv;
    ps = dx2*dx1 + dy2*dy1 + dz1*dz2;
    ps = ps * inv;
    pv = K_FUNC::E_min(pv, 1.);
    pv = K_FUNC::E_max(pv,-1.);
    ps = K_FUNC::E_min(ps, 1.);
    ps = K_FUNC::E_max(ps,-1.);
    
    if (pv >= 0. && ps >= 0.)
      angle[i] = asin(pv)* alp; 
    else if (pv <= 0. && ps >= 0.)
      angle[i] = 360 + asin(pv) * alp;
    else if (pv >= 0. && ps <= 0.)
      angle[i] = 180. - asin(pv) * alp;
    else
      angle[i] = 180. - asin(pv) * alp;
  }
}

//=============================================================================
/* Calcul de l'angle pour un "TRI"
   IN: npts: nbre de vertex dans le maillage TRI
   IN: xt, yt, zt: coord. des pts du TRI
   IN: cn: connectivite du TRI
   OUT: angle pour chaque element et pour chaque arete. 
   Retourne 1: calcul angle correct
   Retourne 0: echec calcul
*/  
//=============================================================================
E_Int K_COMPGEOM::compCurvatureAngleForTri(
  E_Int npts, 
  E_Float* xt, E_Float* yt, E_Float* zt,
  FldArrayI& cn, FldArrayF& angle)
{
  E_Int ind;
  E_Int indp, indm, inda;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6;
  E_Float x, xp, xm, xa;
  E_Float y, yp, ym, ya;
  E_Float z, zp, zm, za;
  E_Float dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3;
  E_Float n1, n2, n3;
  E_Float inv1, inv2;
  E_Float pvx, pvy, pvz;
  E_Float pv, ps1, ps2;
  E_Float pi = 4*atan(1.);
  E_Float alp = 180./pi;

  E_Int nelts = cn.getSize();
  vector< vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs("TRI", npts, cn, cEEN);

  for (E_Int elt1 = 0; elt1 < nelts; elt1++)
  {
    vector<E_Int>& eltVoisins = cEEN[elt1];//liste des tri voisins de elts
    E_Int nvoisins = eltVoisins.size();
    ind1 = cn(elt1,1)-1;
    ind2 = cn(elt1,2)-1;
    ind3 = cn(elt1,3)-1;
    for (E_Int n = 0; n < nvoisins; n++ )
    {
      E_Int elt2 = eltVoisins[n];//indice de l element voisin 
      ind4 = cn(elt2,1)-1;
      ind5 = cn(elt2,2)-1;
      ind6 = cn(elt2,3)-1;
      ind = ind1;
      inda = ind2;
      indp = ind5;
      indm = ind3;
      if (ind4 == ind1)
      {
        if (ind6 == ind2)
        {
          indm = ind3;
          indp = ind5;
          inda = ind2;
        }
        else if (ind5 == ind3)
        {
          indm = ind2;
          indp = ind6;
          inda = ind3;
        }
      }
      else if (ind5 == ind1)
      {
        if (ind4 == ind2)
        {
          indm = ind3;
          indp = ind6;
          inda = ind2;
        }
        else if (ind6 == ind3)
        {
          indm = ind2;
          indp = ind4;
          inda = ind3;
        }
      }
      else if (ind6 == ind1)
      {
        if (ind5 == ind2)
        {
          indm = ind3;
          indp = ind4; 
          inda = ind2;
        }
        else if (ind4 == ind3)
        {
          indm = ind2;
          indp = ind5;
          inda = ind3;
        }
      }
      else if (ind5 == ind2)
      {
        if (ind4 == ind3)
        {
          indm = ind1;
          indp = ind6;
          inda = ind3;
        }
      }
      else if (ind6 == ind2)
      {
        if (ind5 == ind3)
        {
          indm = ind1;
          indp = ind4;
          inda = ind3;
        }
      }
      else if (ind4 == ind2)
      {
        if (ind6 == ind3)
        {
          indm = ind1;
          indp = ind5;
          inda = ind3;
        }
      }
      else
      {
        return 0;
      }
      
      // Compute angle
      x = xt[ind]; xp = xt[indp]; xm = xt[indm]; xa = xt[inda];
      y = yt[ind]; yp = yt[indp]; ym = yt[indm]; ya = yt[inda];
      z = zt[ind]; zp = zt[indp]; zm = zt[indm]; za = zt[inda];
      dx1 = xm-x; dx2 = xp-x; dx3 = xa-x;
      dy1 = ym-y; dy2 = yp-y; dy3 = ya-x;
      dz1 = zm-z; dz2 = zp-z; dz3 = za-x;
      n1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
      n1 = K_FUNC::E_max(n1, 1.e-12); 
      n2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
      n2 = K_FUNC::E_max(n2, 1.e-12);
      n3 = sqrt(dx3*dx3+dy3*dy3+dz3*dz3);
      n3 = K_FUNC::E_max(n3, 1.e-12);
      inv1 = 1./(n1*n2);
      pvx = dy1*dz2 - dz1*dy2;
      pvy = dz1*dx2 - dx1*dz2;
      pvz = dx1*dy2 - dx2*dy1;
      pv = sqrt(pvx*pvx+pvy*pvy+pvz*pvz);
      pv = K_FUNC::E_max(pv, 1.e-12); 
      inv2 = 1./(pv*n3);
      ps1 = dx2*dx1 + dy2*dy1 + dz2*dz1;
      ps1 = ps1 * inv1;
      ps2 = dx3*pvx + dy3*pvy + dz3*pvz;
      ps2 = ps2 * inv2;
      ps2 = K_FUNC::E_min(ps2, 1.);
      ps2 = K_FUNC::E_max(ps2,-1.);
      ps1 = K_FUNC::E_min(ps1, 1.);
      ps1 = K_FUNC::E_max(ps1,-1.);
           
      if (ps2 >= 0.)
        angle(elt1,n+1) = acos(ps1)* alp; 
      else
        angle(elt1,n+1) = 360 - acos(ps1) * alp;
    }
  }
  return 1;
}

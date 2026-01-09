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
# include <math.h>
using namespace std;//dbx
using namespace K_FLD;
//=============================================================================
/* 
   Calcul de la courbure d'un i-array dans le plan (x,y)
   IN: npts: nbre de pts du i-array
   IN: xt, yt, zt: coord. du i-array
   OUT: curv: courbure > 0 dans le sens oppose a la normale
   courbure < 0 sinon.
*/
//=============================================================================
void K_COMPGEOM::compCurvature(E_Int npts,
                               E_Float* xt, E_Float* yt, E_Float* zt,
                               FldArrayF& curv)
{
  FldArrayI node(npts+4);
  E_Int* nodep = node.begin();

  // Traitement des bords: contour ferme ou ouvert ?
  E_Int N = npts-1;
  E_Float dx = xt[0]-xt[npts-1];
  E_Float dy = yt[0]-yt[npts-1];
  E_Float dz = zt[0]-zt[npts-1];
  E_Float dist = dx*dx+dy*dy+dz*dz;
  E_Int open = 0;
  dx = xt[1]-xt[npts-2];
  dy = yt[1]-yt[npts-2];
  dz = zt[1]-zt[npts-2];
  E_Float dist2 = dx*dx+dy*dy+dz*dz;
  if (dist < 1.e-12 && dist2 > 1.e-12) // evite les contours en C 
  {
    for (E_Int i = 0; i <= N; i++) nodep[i+2] = i;
    nodep[0] = N-2; nodep[1] = N-1;
    nodep[N+3] = 1; nodep[N+4] = 2;
  }
  else  // ouvert 
  {
    open = 1;
    for (E_Int i = 0; i <= N; i++) nodep[i+2] = i;
    nodep[0] = 0; nodep[1] = 0;
    nodep[N+3] = N; nodep[N+4] = N;
  }
  E_Int im1, ip1, ii;
  E_Float xp, yp, zp, xpp, ypp, zpp, n1, n2, n3, denom;

  for (E_Int i = 2; i <= N+2; i++)
  {
    im1 = nodep[i-1];
    ip1 = nodep[i+1];
    ii = i-2;
    xp = (xt[ip1]-xt[im1])/2.;
    yp = (yt[ip1]-yt[im1])/2.;
    zp = (zt[ip1]-zt[im1])/2.;
    xpp = (xt[ip1]-2*xt[ii]+xt[im1]);
    ypp = (yt[ip1]-2*yt[ii]+yt[im1]);
    zpp = (zt[ip1]-2*zt[ii]+zt[im1]);

    n1 = (zpp*yp-ypp*zp);
    n2 = (xpp*zp-zpp*xp);
    n3 = (ypp*xp-xpp*yp);
    denom = xp*xp+yp*yp+zp*zp;
    if (K_FUNC::fEqualZero(denom) == true) curv[ii] = 0.;
    else
    {
      denom = denom*denom*denom;
      curv[ii] = sqrt((n1*n1+n2*n2+n3*n3)/denom);
    } 
  }

  if (open == 1)
  {
    curv[0] = curv[1];
    curv[N] = curv[N-1];
  }
  else curv[N] = curv[0];
}

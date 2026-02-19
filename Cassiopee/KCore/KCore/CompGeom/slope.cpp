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
# include <math.h>
# include <stdio.h>
# include "CompGeom/compGeom.h"

// Compute the slope
void K_COMPGEOM::slope(
  const E_Int m,
  const E_Float* x0, const E_Float* y0, const E_Float* z0,
  E_Float* dx, E_Float* dy, E_Float* dz
)
{
  // 2 pts case
  if (m == 2)
  {
    E_Float dxp = x0[1] - x0[0];
    E_Float dyp = y0[1] - y0[0];
    E_Float dzp = z0[1] - z0[0];
    E_Float den = 1.0 / sqrt(dxp*dxp + dyp*dyp + dzp*dzp);
    dx[0] = dxp * den;
    dy[0] = dyp * den;
    dz[0] = dzp * den;
    dx[1] = dx[0];
    dy[1] = dy[0];
    dz[1] = dz[0];
    return;
  }

  // first vertex
  E_Float dxm = x0[1] - x0[0];
  E_Float dxp = x0[2] - x0[1];
  E_Float dym = y0[1] - y0[0];
  E_Float dyp = y0[2] - y0[1];
  E_Float dzm = z0[1] - z0[0];  
  E_Float dzp = z0[2] - z0[1];
  E_Float d2m = dxm*dxm + dym*dym + dzm*dzm;
  E_Float d2p = dxp*dxp + dyp*dyp + dzp*dzp;
  E_Float dm = sqrt(d2m);
  E_Float dp = sqrt(d2p);
  E_Float ambda = (1.0 + dp/dm) * (1.0 + dp/dm);
  E_Float cxm = ambda*dxm - dxm - dxp;
  E_Float cym = ambda*dym - dym - dyp;
  E_Float czm = ambda*dzm - dzm - dzp;
  E_Float denm = 1.0 / sqrt(cxm*cxm + cym*cym + czm*czm);
  dx[0] = cxm * denm;
  dy[0] = cym * denm;
  dz[0] = czm * denm;
 
  // last vertex
  E_Int m1 = m-1;
  dxm = x0[m1-1] - x0[m1-2];
  dxp = x0[m1] - x0[m1-1];
  dym = y0[m1-1] - y0[m1-2];
  dyp = y0[m1] - y0[m1-1];
  dzm = z0[m1-1] - z0[m1-2];
  dzp = z0[m1] - z0[m1-1];
  d2m = dxm*dxm + dym*dym + dzm*dzm;
  d2p = dxp*dxp + dyp*dyp + dzp*dzp;
  dm = sqrt(d2m);
  dp = sqrt(d2p);
  ambda = (1.0 + dm/dp) * (1.0 + dm/dp);
  E_Float cxp = ambda*dxp - dxm - dxp;
  E_Float cyp = ambda*dyp - dym - dyp;
  E_Float czp = ambda*dzp - dzm - dzp;
  E_Float denp = 1.0 / sqrt(cxp*cxp + cyp*cyp + czp*czp);
  dx[m1] = cxp * denp;
  dy[m1] = cyp * denp;
  dz[m1] = czp * denp;

  // other points
  for (E_Int k = 1; k < m-1; k++)
  {
    E_Int kp = k+1;
    E_Int km = k-1;
    dxm = x0[k] - x0[km];
    dxp = x0[kp] - x0[k];
    dym = y0[k] - y0[km];
    dyp = y0[kp] - y0[k];
    dzm = z0[k] - z0[km];
    dzp = z0[kp] - z0[k];
    d2m = dxm*dxm + dym*dym + dzm*dzm;
    d2p = dxp*dxp + dyp*dyp + dzp*dzp;
    ambda = d2p / (d2m + d2p);
    E_Float cx = (1.0 - ambda)*dxm + ambda*dxp;
    E_Float cy = (1.0 - ambda)*dym + ambda*dyp;
    E_Float cz = (1.0 - ambda)*dzm + ambda*dzp;
    E_Float den = 1.0 / sqrt(cx*cx + cy*cy + cz*cz);
    dx[k] = cx * den;
    dy[k] = cy * den;
    dz[k] = cz * den;
  }
}

//=============================================================================
void K_COMPGEOM::slopebar(
  const E_Int npts, const E_Int net,
  const E_Int* cn1, const E_Int* cn2,
  const E_Float* x0, const E_Float* y0, const E_Float* z0,
  E_Float* dx, E_Float* dy, E_Float* dz
)
{
  if (net == 1)
  {
    E_Int ind1 = cn1[0] - 1;
    E_Int ind2 = cn2[0] - 1;
    E_Float dxp = x0[ind2] - x0[ind1];
    E_Float dyp = y0[ind2] - y0[ind1];
    E_Float dzp = z0[ind2] - z0[ind1];
    E_Float den = 1.0 / sqrt(dxp*dxp + dyp*dyp + dzp*dzp);
    dx[0] = dxp * den;
    dy[0] = dyp * den;
    dz[0] = dzp * den;
    dx[1] = dx[0];
    dy[1] = dy[0];
    dz[1] = dz[0];
    return;
  }

  E_Int net1 = net - 1;

  // first vertex
  E_Int ind1 = cn1[0] - 1;
  E_Int ind2 = cn2[0] - 1;
  E_Int ind3 = cn2[1] - 1;
  E_Float dxm = x0[ind2] - x0[ind1];
  E_Float dxp = x0[ind3] - x0[ind2];
  E_Float dym = y0[ind2] - y0[ind1];
  E_Float dyp = y0[ind3] - y0[ind2];
  E_Float dzm = z0[ind2] - z0[ind1];
  E_Float dzp = z0[ind3] - z0[ind2];
  E_Float d2m = dxm*dxm + dym*dym + dzm*dzm;
  E_Float d2p = dxp*dxp + dyp*dyp + dzp*dzp;
  E_Float dm = sqrt(d2m);
  E_Float dp = sqrt(d2p);
  E_Float ambda = (1.0 + dp/dm) * (1.0 + dp/dm);
  E_Float cxm = ambda*dxm - dxm - dxp;
  E_Float cym = ambda*dym - dym - dyp;
  E_Float czm = ambda*dzm - dzm - dzp;
  E_Float denm = 1.0 / sqrt(cxm*cxm + cym*cym + czm*czm);
  dx[ind1] = cxm * denm;
  dy[ind1] = cym * denm;
  dz[ind1] = czm * denm;
 
  // last vertex
  ind1 = cn1[net1-1] - 1;
  ind2 = cn1[net1] - 1;
  ind3 = cn2[net1] - 1;
  dxm = x0[ind2] - x0[ind1];
  dxp = x0[ind3] - x0[ind2];
  dym = y0[ind2] - y0[ind1];
  dyp = y0[ind3] - y0[ind2];
  dzm = z0[ind2] - z0[ind1];
  dzp = z0[ind3] - z0[ind2];
  d2m = dxm*dxm + dym*dym + dzm*dzm;
  d2p = dxp*dxp + dyp*dyp + dzp*dzp;
  dm = sqrt(d2m);
  dp = sqrt(d2p);
  ambda = (1.0 + dm/dp) * (1.0 + dm/dp);
  E_Float cxp = ambda*dxp - dxm - dxp;
  E_Float cyp = ambda*dyp - dym - dyp;
  E_Float czp = ambda*dzp - dzm - dzp;
  E_Float denp = 1.0 / sqrt(cxp*cxp + cyp*cyp + czp*czp);
  dx[ind3] = cxp * denp;
  dy[ind3] = cyp * denp;
  dz[ind3] = czp * denp;


  // other points
  for (E_Int et = 0; et <= net1-1; et++)
  {
    E_Int km = cn1[et] - 1;
    E_Int k  = cn2[et] - 1;
    E_Int kp = cn2[et+1] - 1;
    dxm = x0[k] - x0[km];
    dxp = x0[kp] - x0[k];
    dym = y0[k] - y0[km];
    dyp = y0[kp] - y0[k];
    dzm = z0[k] - z0[km];
    dzp = z0[kp] - z0[k];
    d2m = dxm*dxm + dym*dym + dzm*dzm;
    d2p = dxp*dxp + dyp*dyp + dzp*dzp;
    ambda = d2p / (d2m + d2p);
    E_Float cx = (1.0 - ambda)*dxm + ambda*dxp;
    E_Float cy = (1.0 - ambda)*dym + ambda*dyp;
    E_Float cz = (1.0 - ambda)*dzm + ambda*dzp;
    E_Float den = 1.0 / sqrt(cx*cx + cy*cy + cz*cz);
    dx[k] = cx * den;
    dy[k] = cy * den;
    dz[k] = cz * den;
  }
}

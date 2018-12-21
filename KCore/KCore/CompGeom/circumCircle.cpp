/*    
    Copyright 2013-2019 Onera.

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

//=============================================================================
/* Calcul le rayon du cercle circonscrit */
//=============================================================================
E_Float K_COMPGEOM::circumCircleRadius(E_Float* p1, E_Float* p2, E_Float* p3)
{
  // AB
  E_Float xa = p2[0]-p1[0]; E_Float ya = p2[1]-p1[1]; E_Float za = p2[2]-p1[2];
  // AC
  E_Float xb = p3[0]-p1[0]; E_Float yb = p3[1]-p1[1]; E_Float zb = p3[2]-p1[2];
  // BC
  E_Float xc = p3[0]-p2[0]; E_Float yc = p3[1]-p2[1]; E_Float zc = p3[2]-p2[2];
  
  E_Float a2 = xa*xa + ya*ya + za*za;
  E_Float b2 = xb*xb + yb*yb + zb*zb;
  E_Float c2 = xc*xc + yc*yc + zc*zc;

  // Aire du triangle
  E_Float A = 2*b2*c2 + 2*c2*a2 + 2*a2*b2 - a2*a2 - b2*b2 - c2*c2;
  if (K_FUNC::fEqualZero(A, 1.e-24) == true) return 0.;
  E_Float R = sqrt(a2*b2*c2 / A);
  return R;
}
//=============================================================================
/* Calcul le rayon du cercle inscrit */
//=============================================================================
E_Float K_COMPGEOM::inscribedCircleRadius(
  E_Float* p1, E_Float* p2, E_Float* p3)
{
  E_Float xa = p2[0]-p3[0]; E_Float ya = p2[1]-p3[1]; E_Float za = p2[2]-p3[2];
  E_Float xb = p3[0]-p1[0]; E_Float yb = p3[1]-p1[1]; E_Float zb = p3[2]-p1[2];
  E_Float xc = p1[0]-p2[0]; E_Float yc = p1[1]-p2[1]; E_Float zc = p1[2]-p2[2];
  E_Float a = sqrt(xa*xa + ya*ya + za*za); //dBC
  E_Float b = sqrt(xb*xb + yb*yb + zb*zb); //dAC
  E_Float c = sqrt(xc*xc + yc*yc + zc*zc); //dAB
  
  E_Float peri = a+b+c;
  E_Float aire = compTriangleArea(a, b, c);

  if (K_FUNC::fEqualZero(peri) == true) return 0.;
  else return 2 * aire / peri;
}

//=============================================================================
/* Calcul du cercle circonscrit a un triangle P1P2P3 :
   retourne le centre et le rayon du cercle */
//=============================================================================
E_Int K_COMPGEOM::circumCircle(E_Float* p1, E_Float* p2, E_Float* p3,
                               E_Float* pc, E_Float& R)
{
  // AB
  E_Float xa = p2[0]-p1[0]; E_Float ya = p2[1]-p1[1]; E_Float za = p2[2]-p1[2];
  // AC
  E_Float xb = p3[0]-p1[0]; E_Float yb = p3[1]-p1[1]; E_Float zb = p3[2]-p1[2];
  // BC
  E_Float xc = p3[0]-p2[0]; E_Float yc = p3[1]-p2[1]; E_Float zc = p3[2]-p2[2];

  E_Float a2 = xa*xa + ya*ya + za*za;
  E_Float b2 = xb*xb + yb*yb + zb*zb;
  E_Float c2 = xc*xc + yc*yc + zc*zc;

  E_Float A = 2*b2*c2 + 2*c2*a2 + 2*a2*b2 - a2*a2 - b2*b2 - c2*c2;
  if ( K_FUNC::fEqualZero(A, 1.e-24) == true) 
  { R = 0.; pc[0] = 0; pc[1] = 0; pc[2] = 0; return -1; }

  R = sqrt(a2*b2*c2 / A);

  E_Float nx = ya*zb - za*yb;
  E_Float ny = za*xb - xa*zb;
  E_Float nz = xa*yb - ya*xb;
  E_Float tx = ya*nz - za*ny;
  E_Float ty = za*nx - xa*nz;
  E_Float tz = xa*ny - ya*nx;
  E_Float norm = tx*tx + ty*ty + tz*tz;
  E_Float normi;
  if (K_FUNC::fEqualZero(norm) != true) normi = 1./sqrt(norm);
  else normi = 1.e10;
  tx = tx*normi; ty = ty*normi; tz = tz*normi;
  E_Float alpha = R*R - (xa*xa+ya*ya+za*za)*0.25;
  alpha = sqrt(alpha);
  pc[0] = 0.5*(p1[0]+p2[0]) + alpha*tx;
  pc[1] = 0.5*(p1[1]+p2[1]) + alpha*ty;
  pc[2] = 0.5*(p1[2]+p2[2]) + alpha*tz;
 
  E_Float l = (p3[0]-pc[0])*(p3[0]-pc[0]) + (p3[1]-pc[1])*(p3[1]-pc[1]) +
    (p3[2]-pc[2])*(p3[2]-pc[2]);
  if (K_FUNC::fEqualZero(l - R*R, 1.e-10) == true) return 0;
  
  pc[0] = 0.5*(p1[0]+p2[0]) - alpha*tx;
  pc[1] = 0.5*(p1[1]+p2[1]) - alpha*ty;
  pc[2] = 0.5*(p1[2]+p2[2]) - alpha*tz;
  l = (p3[0]-pc[0])*(p3[0]-pc[0]) + (p3[1]-pc[1])*(p3[1]-pc[1]) +
    (p3[2]-pc[2])*(p3[2]-pc[2]);
  if (K_FUNC::fEqualZero(l - R*R, 1.e-10) == true) return 0;
  return -1;
}

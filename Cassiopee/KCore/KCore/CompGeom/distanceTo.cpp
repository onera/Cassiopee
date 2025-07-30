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

# include "kcore.h"

//=============================================================================
/* 
   Calcul le projete de p sur le plan defini par le triangle.
   Retourne la distance entre p et son projete.
   
   IN: p0, p1, p2: triangle coordinates
   IN: p: tested point
   IN: treatment: traitement effectue si p ne se projete pas dans le triangle
   treatment=0, retourne le projete sur le plan du triangle
   treatment=1, retourne le sommet du triangle le plus proche du projete
   treatment=2, retourne le point de l'edge le plus proche
   OUT: dist2: square distance to projected point
   OUT: in: true, projection falls in the triangle.
            false, projection falls out of triangle.
   OUT: xp, yp, zp: coordinates of projected point.
   OUT: sigma0, sigma1: coord. parametrique de P dans le triangle
   return:
   0: OK
   -1: Failed.
*/
//=============================================================================
E_Int K_COMPGEOM::distanceToTriangle(
  E_Float* p0, E_Float* p1, E_Float* p2,
  E_Float* p, E_Int treatment,
  E_Float& dist2, E_Bool& in, 
  E_Float& xp, E_Float& yp, E_Float& zp,
  E_Float& sigma0, E_Float& sigma1)
{
  E_Float e0x, e0y, e0z, e1x, e1y, e1z;
  E_Float nx, ny, nz, t, px, py, pz;
  E_Float det, deti;

  in = false;
  dist2 = 0.;
  xp = p[0]; yp = p[1]; zp = p[2];
  e0x = p1[0]-p0[0]; e0y = p1[1]-p0[1]; e0z = p1[2]-p0[2];
  e1x = p2[0]-p0[0]; e1y = p2[1]-p0[1]; e1z = p2[2]-p0[2];
  nx = e0y*e1z-e0z*e1y;
  ny = e0z*e1x-e0x*e1z;
  nz = e0x*e1y-e0y*e1x;

  det = e0x*e1y*nz - e0x*ny*e1z - e0y*e1x*nz 
    + e0y*nx*e1z + e0z*e1x*ny - e0z*e1y*nx;

  if (K_FUNC::fEqualZero(det, 1.e-24) == true) return -1;

  px = p[0]-p0[0]; py = p[1]-p0[1]; pz = p[2]-p0[2];

  deti = 1./det;
  sigma1 = (e1y*nz - ny*e1z)*px - (e1x*nz-nx*e1z)*py + (e1x*ny-nx*e1y)*pz;
  sigma1 = sigma1 * deti;
  sigma0 = -(e0y*nz-ny*e0z)*px + (e0x*nz-nx*e0z)*py - (e0x*ny-nx*e0y)*pz;
  sigma0 = sigma0 * deti;
  t = -(e0y*e1z-e1y*e0z)*px + (e0x*e1z-e1x*e0z)*py - (e0x*e1y-e1x*e0y)*pz; 
  t = t * deti;

  xp = p[0] + t*nx;
  yp = p[1] + t*ny;
  zp = p[2] + t*nz;
  dist2 = (xp-p[0])*(xp-p[0]) + (yp-p[1])*(yp-p[1]) + (zp-p[2])*(zp-p[2]);
  E_Float eps = K_CONST::E_GEOM_CUTOFF;

  if (sigma0 >= -eps && 
      sigma1 >= -eps && sigma0+sigma1 <= 1.+2*eps)
  {
    in = true;
  }
  else if (treatment == 1)
  {
    // retourne dist = min des distances aux sommets des triangles
    xp = p0[0]; yp = p0[1]; zp = p0[2];
    E_Float d0 = (p[0]-p0[0])*(p[0]-p0[0]) + 
    (p[1]-p0[1])*(p[1]-p0[1])+(p[2]-p0[2])*(p[2]-p0[2]);  
    E_Float d1 =  (p[0]-p1[0])*(p[0]-p1[0]) + 
    (p[1]-p1[1])*(p[1]-p1[1])+(p[2]-p1[2])*(p[2]-p1[2]);  
    E_Float d2 = (p[0]-p2[0])*(p[0]-p2[0]) + 
    (p[1]-p2[1])*(p[1]-p2[1])+(p[2]-p2[2])*(p[2]-p2[2]);  
    if (d1 < d0) 
    {
      if (d1 < d2)
      {xp = p1[0]; yp = p1[1]; zp = p1[2]; dist2 = d1;}
      else //d1 > d2
      {xp = p2[0]; yp = p2[1]; zp = p2[2]; dist2 = d2;}
    }
    else //d1 > d0
    {
      if (d0 < d2)
      {xp = p0[0]; yp = p0[1]; zp = p0[2]; dist2 = d0;}
      else 
      {xp = p0[0]; yp = p0[1]; zp = p0[2]; dist2 = d2;}
    }
  }
  else if (treatment == 2)
  {
    E_Float d0, d1, d2, d3;
    E_Float xp1, yp1, zp1, xp2, yp2, zp2, xp3, yp3, zp3;
    E_Bool in1, in2, in3;
    distanceToBar(p0, p1, p, 0, xp1, yp1, zp1, in1, d1);
    distanceToBar(p0, p2, p, 0, xp2, yp2, zp2, in2, d2);
    distanceToBar(p1, p2, p, 0, xp3, yp3, zp3, in3, d3);
    if (in1 == false) d1 = 1.e6; 
    if (in2 == false) d2 = 1.e6;
    if (in3 == false) d3 = 1.e6;
    if (d1 < d2 && d1 < d3)
    {
      xp = xp1; yp = yp1; zp = zp1; dist2 = d1;
    }
    else if (d2 < d1 && d2 < d3)
    {
      xp = xp2; yp = yp2; zp = zp2; dist2 = d2;
    }
    else
    {
      xp = xp3; yp = yp3; zp = zp3; dist2 = d3;
    }
    if (in1 == false && in2 == false && in3 == false)
    {
      xp = p0[0]; yp = p0[1]; zp = p0[2];
      d0 = (p[0]-p0[0])*(p[0]-p0[0]) + 
        (p[1]-p0[1])*(p[1]-p0[1])+(p[2]-p0[2])*(p[2]-p0[2]);  
      d1 = (p[0]-p1[0])*(p[0]-p1[0]) + 
        (p[1]-p1[1])*(p[1]-p1[1])+(p[2]-p1[2])*(p[2]-p1[2]);  
      d2 = (p[0]-p2[0])*(p[0]-p2[0]) + 
        (p[1]-p2[1])*(p[1]-p2[1])+(p[2]-p2[2])*(p[2]-p2[2]);  
      if (d1 < d0) 
      {
        if (d1 < d2)
        {xp = p1[0]; yp = p1[1]; zp = p1[2]; dist2 = d1;}
        else //d1 > d2
        {xp = p2[0]; yp = p2[1]; zp = p2[2]; dist2 = d2;}
      }
      else //d1 > d0
      {
        if (d0 < d2)
        {xp = p0[0]; yp = p0[1]; zp = p0[2]; dist2 = d0;}
        else 
        {xp = p2[0]; yp = p2[1]; zp = p2[2]; dist2 = d2;}
      }
    }
  }
  return 0;
}

//=============================================================================
/* Distance d'un point a une BAR (AB)
   IN: p: tested point
   IN: pA, pB: BAR (segment)
   IN: treatement: traitement effectue si p ne se projete pas dans la BAR
   treatment=0, retourne le projete sur la BAR prolongee a l'infini
   treatment=1, retourne le sommet de la BAR le plus proche du projete
   OUT: xp, yp, zp: coord du point projete
   OUT: in: true, le pt projete est dans la BAR,
            false, le pt est dehors
   OUT: dist2: distance au carre du pt p au pt projete.
   Retourne 0: OK
   -1: Failed.
*/
//=============================================================================
E_Int K_COMPGEOM::distanceToBar(E_Float* pA, E_Float* pB,
                                E_Float* p, E_Int treatment,
                                E_Float& xp, E_Float& yp, E_Float& zp,
                                E_Bool& in, E_Float& dist2)
{
  in = false;

  E_Float dxAP = p[0]-pA[0]; E_Float dxAB = pB[0]-pA[0]; 
  E_Float dyAP = p[1]-pA[1]; E_Float dyAB = pB[1]-pA[1]; 
  E_Float dzAP = p[2]-pA[2]; E_Float dzAB = pB[2]-pA[2]; 
  E_Float scal = dxAP*dxAB + dyAP*dyAB + dzAP*dzAB;
  
  E_Float dAB2 = dxAB*dxAB+dyAB*dyAB+dzAB*dzAB;
  if (K_FUNC::fEqualZero(dAB2, 1.e-24) == true)
  {
    if (K_FUNC::fEqualZero(dxAP*dxAP + dyAP*dyAP + dzAP*dzAP, 1.e-24) == true)
    {
      in = true; xp = pA[0]; yp = pA[1]; zp = pA[2]; return 0;
    }
    else
    {
      in = false; xp = p[0]; yp = p[1]; zp = p[2]; return -1;
    }
  }
  
  E_Float alpha = scal / dAB2;
  xp = pA[0] + alpha*dxAB;
  yp = pA[1] + alpha*dyAB;
  zp = pA[2] + alpha*dzAB;
  dist2 = (p[0]-xp)*(p[0]-xp) + (p[1]-yp)*(p[1]-yp) + (p[2]-zp)*(p[2]-zp);
  if (alpha < 0. || alpha > 1.) 
  {
    in = false;
    if (treatment == 1)
    {
      if (alpha < 0.5) { xp = pA[0]; yp = pA[1]; zp = pA[2]; }
      else { xp = pB[0]; yp = pB[1]; zp = pB[2]; }
      dist2 = (p[0]-xp)*(p[0]-xp) + (p[1]-yp)*(p[1]-yp) + (p[2]-zp)*(p[2]-zp);
    }
  }
  else in = true;
  return 0;
}

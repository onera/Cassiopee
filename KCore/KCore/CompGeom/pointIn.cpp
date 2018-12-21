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

#include "kcore.h"

#define isLeft(p0,p1,p2) (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1])
//=============================================================================
/*
  IN: p0, p1: segment
  IN: p: pt a tester
  IN: eps: tolerance
  OUT: 1: p est le pt p0
  OUT: 2: p est le pt p1
  OUT: 3: p est dans p0-p1
  OUT: 0: sinon
*/
//=============================================================================
E_Int K_COMPGEOM::pointInSegment(E_Float* p0, E_Float* p1, 
                                 E_Float* p, E_Float eps)
{
  E_Float dx = p1[0]-p0[0];
  E_Float dy = p1[1]-p0[1];
  E_Float dz = p1[2]-p0[2];
  E_Float dx0 = p[0]-p0[0];
  E_Float dy0 = p[1]-p0[1];
  E_Float dz0 = p[2]-p0[2];
  E_Float dx1 = p1[0]-p[0];
  E_Float dy1 = p1[1]-p[1];
  E_Float dz1 = p1[2]-p[2];
  
  E_Float scal = dx0*dx + dy0*dy + dz0*dz;
  E_Float norm = dx*dx + dy*dy + dz*dz; 
  E_Float norm0 = dx0*dx0 + dy0*dy0 + dz0*dz0;
  E_Float norm1 = dx1*dx1 + dy1*dy1 + dz1*dz1;

  eps = eps*sqrt(norm);
  E_Float eps2 = eps*eps;

  if (norm0 < eps2) return 1;
  if (norm1 < eps2) return 2;
  if (K_FUNC::E_abs(scal*scal - norm*norm0) > eps2) return 0; // non colineaires

  E_Float a;
  if (K_FUNC::E_abs(dx) > eps) a = dx0/dx;
  else if (K_FUNC::E_abs(dy) > eps) a = dy0/dy;
  else if (K_FUNC::E_abs(dz) > eps) a = dz0/dz;
  else
  { // p0 = p1
    if (norm0 > eps2) return 0;
    else return 1;
  }

  if (a <= eps && a >= -eps) return 1;
  else if (a >= 1-eps && a <= 1+eps) return 2;
  else if (a >= eps && a <= 1-eps) return 3; // IN
  return 0;
}
//=============================================================================
/* 
   IN: p1, p2, p3: triangle coordinates
   IN: p: tested point
   IN: Ball: null or circum-circle of triangle (xR,yR,R)
   IN: BB: null or bounding box of triangle (xmin,ymin,xmax,ymax)

   return:
   0: out
   1: in   
*/
//=============================================================================
E_Int K_COMPGEOM::pointInTriangle2D(E_Float* p0, E_Float* p1, E_Float* p2,
                                    E_Float* p,
                                    E_Float* Ball, E_Float* BB)
{
  E_Float x = p[0];
  E_Float y = p[1];

  if (BB != NULL)
  {
    E_Float xmin = BB[0];
    E_Float ymin = BB[1];
    E_Float xmax = BB[2];
    E_Float ymax = BB[3];
    if (x < xmin || x > xmax || y < ymin || y > ymax) return 0;
  }

  if (Ball != NULL)
  {
    E_Float xR = BB[0];
    E_Float yR = BB[1];
    E_Float R  = BB[2];
    E_Float d = (x - xR)*(x - xR)+(y - yR)*(y - yR);
    if (d > R*R) return 0;
  }
  
  E_Int wn = 0;

  // Boucle sur les aretes du triangle
  if (p0[1] <= y)
  {
    if (p1[1] > y)
      if (isLeft(p0, p1, p) > 0) wn++;
  }
  else
  {
    if (p1[1] <= y)
      if (isLeft( p0, p1, p) < 0) wn--;
  }
  
  if (p1[1] <= y)
  {
    if (p2[1] > y)
      if (isLeft(p1, p2, p) > 0) wn++;
  }
  else
  {
    if (p2[1] <= y)
      if (isLeft( p1, p2, p) < 0) wn--;
  }
  
  if (p2[1] <= y)
  {
    if (p0[1] > y)
      if (isLeft(p2, p0, p) > 0) wn++;
  }
  else
  {
    if (p0[1] <= y)
      if (isLeft( p2, p0, p) < 0) wn--;
  }
  
  return (wn == 0 ? 0 : 1);
}
//=============================================================================
/* 
   IN: p1, p2, p3: triangle coordinates
   IN: p: tested point
   IN: Ball: null or circum-circle of triangle (xR,yR,R)
   IN: BB: null or bounding box of triangle (xmin,ymin,xmax,ymax)

   return:
   0: out
   1: in 
*/
//=============================================================================
E_Int K_COMPGEOM::pointInPolygon2D(E_Float* xt, E_Float* yt, E_Float* zt,
                                   K_FLD::FldArrayI& connect,
                                   E_Float* p,
                                   E_Float* Ball, E_Float* BB)
{
  E_Int ind1, ind2;
  E_Float x = p[0];
  E_Float y = p[1];
  E_Float p0[3];
  E_Float p1[3]; 

  if (BB != NULL)
  {
    E_Float xmin = BB[0];
    E_Float ymin = BB[1];
    E_Float xmax = BB[2];
    E_Float ymax = BB[3];
    if (x < xmin || x > xmax || y < ymin || y > ymax) return 0;
  }

  if (Ball != NULL)
  {
    E_Float xR = BB[0];
    E_Float yR = BB[1];
    E_Float R  = BB[2];
    E_Float d = (x - xR)*(x - xR)+(y - yR)*(y - yR);
    if (d > R*R) return 0;
  }
  
  E_Int wn = 0;

  // Boucle sur les aretes du polygone
  E_Int na = connect.getSize();

  for (E_Int i = 0; i < na; i++)
  {
    ind1 = connect(i,1)-1;
    ind2 = connect(i,2)-1;
    
    p0[0] = xt[ind1]; p0[1] = yt[ind1]; p0[2] = zt[ind1];
    p1[0] = xt[ind2]; p1[1] = yt[ind2]; p1[2] = zt[ind2];
      
    if (p0[1] <= y) 
    {
      if (p1[1] > y)
      {
        if (isLeft(p0, p1, p) > 0) wn++;
      }
    }
    else
    {
      if (p1[1] <= y)
      {
        if (isLeft( p0, p1, p) < 0) wn--;          
      }
    }
  }
  
  return (wn == 0 ? 0 : 1);
}
//=============================================================================
/* 
   IN: p1, p2, p3: triangle coordinates
   IN: p: tested point
   IN: epsilon: tolerance
   IN: Ball: null or circum-circle of triangle (xR,yR,R)
   IN: BB: null or bounding box of triangle (xmin,ymin,xmax,ymax)

   return:
   0: out
   1: in 
*/
//=============================================================================
E_Int K_COMPGEOM::pointInPolygon2D(E_Float* xt, E_Float* yt, E_Float* zt,
                                   K_FLD::FldArrayI& connect,
                                   E_Float* p, E_Float epsilon,
                                   E_Float* Ball, E_Float* BB)
{
  E_Int ret = 0;
  E_Float pt[3]; 
  pt[0] = p[0]; pt[1] = p[1]; pt[2] = p[2];
  ret = pointInPolygon2D(xt, yt, zt, connect,
                         pt,
                         Ball, BB);
  if (ret == 1) return 1;
  pt[0] = p[0] + epsilon;
  ret = pointInPolygon2D(xt, yt, zt, connect,
                         pt,
                         Ball, BB);
  if (ret == 1) return 1;
  pt[0] = p[0] - epsilon;
  ret = pointInPolygon2D(xt, yt, zt, connect,
                         pt,
                         Ball, BB);
  if (ret == 1) return 1;
  pt[0] = p[0]; pt[1] = p[1] + epsilon;
  ret = pointInPolygon2D(xt, yt, zt, connect,
                         pt,
                         Ball, BB);
  if (ret == 1) return 1;
  pt[1] = p[1] - epsilon;
  ret = pointInPolygon2D(xt, yt, zt, connect,
                         pt,
                         Ball, BB);
  if (ret == 1) return 1;
  return 0;
}
//=============================================================================
E_Int K_COMPGEOM::pointInBB(E_Float xmin, E_Float ymin, E_Float zmin,
                            E_Float xmax, E_Float ymax, E_Float zmax,
                            E_Float* p, E_Float tol)
{
  if (p[0] < xmin - tol) return 0;
  if (p[0] > xmax + tol) return 0;
  if (p[1] < ymin - tol) return 0;
  if (p[1] > ymax + tol) return 0;
  if (p[2] < zmin - tol) return 0;
  if (p[2] > zmax + tol) return 0;
  return 1;
}

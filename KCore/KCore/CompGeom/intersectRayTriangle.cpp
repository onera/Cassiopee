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

#include "compGeom.h"

using namespace K_FLD;
using namespace K_FUNC;
using namespace std;

#define dot(u,v) ( u[0]*v[0] + u[1]*v[1] + u[2]*v[2] )

//=============================================================================
E_Int K_COMPGEOM::intersectRayTriangle(E_Float* p1, E_Float* p2, E_Float* p3,
                                       E_Float* pR0, E_Float* pR1,
                                       E_Float* pi)
{
  E_Float tol = 1.e-6;
  E_Float small = 1.e-18;
    
  E_Float u[3];
  E_Float v[3];
  E_Float n[3];
  E_Float dir[3];
  E_Float w[3];
  E_Float w0[3];
  E_Float a, b;
  
  // Aretes du triangle + normal au plan du triangle
  u[0] = p2[0] - p1[0];
  u[1] = p2[1] - p1[1];
  u[2] = p2[2] - p1[2];
  v[0] = p3[0] - p1[0];
  v[1] = p3[1] - p1[1];
  v[2] = p3[2] - p1[2];
  n[0] = u[1]*v[2]-u[2]*v[1];
  n[1] = u[2]*v[0]-u[0]*v[2];
  n[2] = u[0]*v[1]-u[1]*v[0];
  if (fEqualZero(n[0]*n[0] + n[1]*n[1] + n[2]*n[2], 1.e-28) == true)
    return -1;

  dir[0] = pR1[0] - pR0[0];
  dir[1] = pR1[1] - pR0[1]; 
  dir[2] = pR1[2] - pR0[2];
  
  w0[0] = pR0[0] - p1[0];
  w0[1] = pR0[1] - p1[1];
  w0[2] = pR0[2] - p1[2];
  
  a = -dot(n, w0);
  b = dot(n, dir);

  if (K_FUNC::E_abs(b) < small)
  {
    if (fEqualZero(a) == true) return 2;
    else return 0;
  }

  // Intersection
  E_Float r, s, t;
  r = a / b;
  
  pi[0] = pR0[0] + r * dir[0];
  pi[1] = pR0[1] + r * dir[1];
  pi[2] = pR0[2] + r * dir[2];
  
  // pi est dans le triangle?
  E_Float uu, uv, vv, wu, wv, D;
  uu = dot(u,u);
  uv = dot(u,v);
  vv = dot(v,v);
  w[0] = pi[0] - p1[0];
  w[1] = pi[1] - p1[1];
  w[2] = pi[2] - p1[2];
  wu = dot(w,u);
  wv = dot(w,v);
  D = uv * uv - uu * vv;
  
  // Test les coordonnees parametriques
  s = (uv * wv - vv * wu) / D;
  if (s < -tol || s > 1.+tol) return 0;
  t = (uv * wu - uu * wv) / D;
  if (t < -tol || (s+t) > 1.+tol) return 0;
  
  return 1;  
}

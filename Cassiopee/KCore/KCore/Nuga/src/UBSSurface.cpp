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
//Authors : Sam Landier (sam.landier@onera.fr)

#include "Nuga/include/UBSSurface.h"
#include <algorithm>

UBSSurface* UBSSurface::buildUBSSurface
(const K_FLD::FloatArray& pos, E_Int nj)
{
  if ((pos.cols() == 0) || (nj == 0)) return 0; // Error
  return new UBSSurface(pos, nj);
}

UBSSurface::UBSSurface(const K_FLD::FloatArray& pos, E_Int nj)
:_pos(pos), _Umin(0.), _Umax(nj+1), _Vmin(0.), _Vmax((pos.cols() / nj) + 1)
{
  // Set the control points.
  __setControlPoints(_pos, nj);

  // Store pointers to the base functions.
  _base.push_back(&__N0);
  _base.push_back(&__N1);
  _base.push_back(&__N2);
  _base.push_back(&__N3);

  // Store pointers to the first derivative of the base functions.
  _baseD1.push_back(&__DN0);
  _baseD1.push_back(&__DN1);
  _baseD1.push_back(&__DN2);
  _baseD1.push_back(&__DN3);

  // Store pointers to the second derivative of the base functions.
  _baseD2.push_back(&__D2N0);
  _baseD2.push_back(&__D2N1);
  _baseD2.push_back(&__D2N2);
  _baseD2.push_back(&__D2N3);
}

UBSSurface::~UBSSurface(void)
{
}

void UBSSurface::point (E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _base, v, _base, P);
}

void UBSSurface::DU1(E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _baseD1, v, _base, P);
}


void UBSSurface::DU2(E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _baseD2, v, _base, P);
}


void UBSSurface::DV1(E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _base, v, _baseD1, P);
}


void UBSSurface::DV2(E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _base, v, _baseD2, P);
}


void UBSSurface::DUV(E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _baseD1, v, _baseD1, P);
}

void UBSSurface::DAUV1(E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _baseD1, v, _base, P); // DU1
  __eval (u, _base, v, _baseD1, P+3); // DV1
}

void UBSSurface::DAUV2(E_Float u, E_Float v, E_Float* P) const
{
  __eval (u, _baseD1, v, _base,   P); // DU1
  __eval (u,  _base,  v, _baseD1, P+3); // DV1
  __eval (u, _baseD2, v, _base,   P+6); // DU2
  __eval (u, _base,   v, _baseD2, P+9); // DV2
  __eval (u, _baseD1, v, _baseD1, P+12); // DUV
}

void UBSSurface::__eval (E_Float u, const std::vector<pBaseFunc>& FUs, 
                         E_Float v, const std::vector<pBaseFunc>& FVs, 
                         E_Float* P) const
{
  P[0] = P[1] = P[2] = 0.;

  // Ensure u and v are in the correct interval.
  u = std::min(_Umax, u);
  u = std::max(_Umin, u);
  v = std::min(_Vmax, v);
  v = std::max(_Vmin, v);
 
  E_Int i = (v == _Vmax)? E_Int(_Vmax - 1.) : E_Int(v);
  E_Int j = (u == _Umax)? E_Int(_Umax - 1.) : E_Int(u);
  
  // Normalize.
  u = u-j;
  v = v-i;

  ///
  for (E_Int k = j; k < j+4; ++k)
  {
    pBaseFunc Nku = FUs[k-j];
    E_Float xu = Nku(u);

    if (::fabs(xu) < EPSILON)
      continue;

    for (E_Int l = i; l < i+4; ++l)
    {
      pBaseFunc Nlv = FVs[l-i];
      E_Float    yv = Nlv(v);

      if (::fabs(yv) < EPSILON)
        continue;

      E_Int     nkl = _ctrlPts(l,k);

      for (E_Int i = 0; i < 3; ++i)
        P[i] += xu*yv*_pos(i,nkl);
    }
  }
}

void UBSSurface::__setControlPoints(const K_FLD::FloatArray& pos, E_Int nj)
{
  if (pos.cols() == 0) return;

  E_Int ni = pos.cols() / nj;

  /* Extend the control points to force the spline
     to interpolate the boundary points.*/

  E_Int m = ni;//ctrlPoints.rows();
  E_Int n = nj;//ctrlPoints.cols();
  E_Int none = IDX_NONE;

  E_Int r = 2;
  E_Int s = 2*r;

  _ctrlPts.resize(m+s, n+s, &none);

  // Set the input points.
  for (E_Int i = 0; i < m; ++i)
    for (E_Int j = 0; j < n; ++j)
      _ctrlPts(i+r,j+r) = j + i * nj;

  for (E_Int i = 0; i < r; ++i)
    for (E_Int j = 0; j < r; ++j)
      _ctrlPts(i,j) = 0;//j + i * nj ctrlPoints(0,0);

  for (E_Int i = 0; i < r; ++i)
    for (E_Int j = 0; j < r; ++j)
      _ctrlPts(i,j+n+r) = n-1;//j + i * nj ctrlPoints(0,n-1);

  for (E_Int i = 0; i < r; ++i)
    for (E_Int j = 0; j < r; ++j)
      _ctrlPts(i+m+r,j) = (m-1)*nj;//j + i * nj ctrlPoints(m-1,0);

  for (E_Int i = 0; i < r; ++i)
    for (E_Int j = 0; j < r; ++j)
      _ctrlPts(i+m+r,j+n+r) = (n-1) + (m-1) * nj;//j + i * nj ctrlPoints(m-1,n-1);


  // Add the extra points (4 blocks).
  E_Int P0;
  // Left
  for (E_Int i = 0; i < m; ++i)
  {
     P0 = i*nj;//j + i * nj ctrlPoints(i,0);
    for (E_Int j = 0; j < r; ++j)
      _ctrlPts(i+r,j) = P0;
  }
  // Right
  for (E_Int i = 0; i < m; ++i)
  {
     P0 =  (n-1) + i*nj;//j + i * nj ctrlPoints(i,n-1);
    for (E_Int j = n; j < n+r; ++j)
      _ctrlPts(i+r,j+r) = P0;
  }
  //Top
  for (E_Int j = 0; j < n; ++j)
  {
    P0 = j;//j + i * nj ctrlPoints(0,j);
    for (E_Int i = 0; i < r; ++i)
        _ctrlPts(i,j+r) = P0;
  }
  //Bottom
  for (E_Int j = 0; j < n; ++j)
  {
    P0 = j + (m-1)*nj;//j + i * nj ctrlPoints(m-1,j);
    for (E_Int i = m; i < m+r; ++i)
        _ctrlPts(i+r,j+r) = P0;
  }
}

void UBSSurface::triangulate(K_FLD::FloatArray& pos, 
                             K_FLD::IntArray& connectT3)
{
  E_Float P0[3];
  E_Int   Ti[3];

  connectT3.clear();
  pos.clear();
  K_FLD::IntArray net;
  net.resize(E_Int(_Vmax)+1, E_Int(_Umax)+1);

  // Create net points.
  for (E_Int v = 0; v <= _Vmax; ++v)
  {
    for (E_Int u =0; u <= _Umax; ++u)
    {
      point(u,v, P0);
      pos.pushBack(P0, P0+3);
      net(v,u) = pos.cols()-1;
    }
  }

  // Create the connectivity.
  for (E_Int v = 0; v < _Vmax; ++v)
  {
    for (E_Int u =0; u < _Umax; ++u)
    {
      Ti[0] = net(v,u);
      Ti[1] = net(v,u+1);
      Ti[2] = net(v+1,u);

      connectT3.pushBack(Ti, Ti+3);

      Ti[0] = net(v+1,u+1);
      Ti[1] = net(v+1,u);
      Ti[2] = net(v,u+1);

      connectT3.pushBack(Ti, Ti+3);
    }
  }
}
/*
void 
UBSSurface::quadrangulate
(K_FLD::FloatArray& pos, K_FLD::IntArray& connectQ4)
{
  E_Float P0[3];
  E_Int   Ti[3];

  connectT3.clear();
  pos.clear();
  K_FLD::IntArray net;
  net.resize(E_Int(_Vmax)+1, E_Int(_Umax)+1);

  // Create net points.
  for (E_Int v = 0; v <= _Vmax; ++v)
  {
    for (E_Int u =0; u <= _Umax; ++u)
    {
      point(u,v, P0);
      pos.pushBack(P0, P0+3);
      net(v,u) = pos.cols()-1;
    }
  }

  // Create the connectivity.
  for (E_Int v = 0; v < _Vmax; ++v)
  {
    for (E_Int u =0; u < _Umax; ++u)
    {
      Ti[0] = net(v,u);
      Ti[1] = net(v,u+1);
      Ti[2] = net(v+1,u);

      connectT3.pushBack(Ti, Ti+3);

      Ti[0] = net(v+1,u+1);
      Ti[1] = net(v+1,u);
      Ti[2] = net(v,u+1);

      connectT3.pushBack(Ti, Ti+3);
    }
  }
}
*/

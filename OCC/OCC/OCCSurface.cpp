/*    
    Copyright 2013-2018 Onera.

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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "OCCSurface.h"
#include "BRep_Tool.hxx"
#include "GeomLib_Tool.hxx"
#include "Fld/ArrayAccessor.h"
#include "Search/KdTree.h"
#include "MeshElement/Edge.h"
#include "Connect/merge.h"
#include "Connect/IdTool.h"

///
K_OCC::OCCSurface::OCCSurface(const TopoDS_Face& F, E_Int id) :_F(F), _surface(BRep_Tool::Surface (F)), _parent(id), _normalize_domain(true){}

///
K_OCC::OCCSurface::OCCSurface(const OCCSurface& rhs)
:_F(rhs._F), _surface(rhs._surface), _parent(rhs._parent), _edges(rhs._edges), _U0(rhs._U0), _U1(rhs._U1), _V0(rhs._V0), _V1(rhs._V1), _normalize_domain(rhs._normalize_domain){}

///
K_OCC::OCCSurface::~OCCSurface() {}

///
E_Int
K_OCC::OCCSurface::parameters
(const K_FLD::FloatArray&coord3D, K_FLD::FloatArray& UVs)
{
  UVs.clear();
    
  E_Int err(0), sz(coord3D.cols());
  E_Float UV[2];
  UVs.resize(2, sz, 0.);
  
  for (size_t i=0; (i < sz) && !err; ++i)
  {
    err = parameters(coord3D.col(i), UV[0], UV[1]);
    UVs(0,i)=UV[0];
    UVs(1,i)=UV[1];
  }
  return err;
}

///
E_Int
K_OCC::OCCSurface::parameters(const E_Float* pt, E_Float & u, E_Float& v)
{
  u=v=-1;
    
  gp_Pnt Point;
  Point.SetCoord(pt[0], pt[1], pt[2]);
  
  E_Int ok = GeomLib_Tool::Parameters(_surface, Point, 1.e+100/*dummytol*/, u, v); 
  
  if (!ok || u==-1 || v== -1)
    return 1;
  return 0;
}

E_Int
K_OCC::OCCSurface::parametersSample(const K_FLD::FloatArray&coord3D, K_FLD::FloatArray& UVs)
{
  E_Float UV[2];
  E_Int err(0), sz(coord3D.cols());
  
  UVs.clear();
  UVs.resize(2, sz, 0.);
  
  //1 do a sampling in both real and parametric space.
  K_FLD::FloatArray pos3D, pos2D;
  __sample_contour(coord3D.cols()/*sampling rate (ASSUME that only contour points are inside...*/, pos3D, pos2D);
  
  //2 for each point in coord3D, find the closest point in the sampling pos3D and Interpolate parameters from pos2D
  // build the KdTree.
  K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(pos3D);
  K_SEARCH::KdTree<K_FLD::FloatArray> tree(crdA);
  size_t szSmpl(pos3D.cols());
  E_Int Nc, N[2]; // N[0] its predecessor, N[1] its successor
  E_Float d2[2], lambda[2]; //d2[0] for predecessor, d2[1] for successor
  
  for (size_t i = 0; i < sz; ++i)
  {
    const E_Float* Pt = coord3D.col(i);
    Nc = tree.getClosest(Pt);
    
    //init
    N[0] = (Nc > 0) ? Nc-1 : szSmpl-1; //prior
    N[1] = (Nc+1)%szSmpl;                // next
      
    //
    d2[0]=K_MESH::Edge::edgePointMinDistance2<3>(pos3D.col(N[0]), pos3D.col(Nc), Pt, lambda[0]);
    d2[1]=K_MESH::Edge::edgePointMinDistance2<3>(pos3D.col(Nc), pos3D.col(N[1]), Pt, lambda[1]);
    
    // Linear interpolation
    if (d2[0] < d2[1])
      K_FUNC::sum<2>(1.-lambda[0], pos2D.col(N[0]), lambda[0], pos2D.col(Nc), UVs.col(i));
    else
      K_FUNC::sum<2>(1.-lambda[1], pos2D.col(Nc), lambda[1], pos2D.col(N[1]), UVs.col(i)); 
  }
  
  return err;
}

E_Int K_OCC::OCCSurface::__sample_contour(E_Int Nsample, K_FLD::FloatArray& pos3D, K_FLD::FloatArray& pos2D)
{
  E_Float U[2], Pt[3];
  
  //
  U[1]=_V0;
  for (size_t a=0; a < Nsample-1; ++a) // up to Nsample-1 to avoid to create duplicates
  {
    U[0]=_U0 + (_U1-_U0)*a/(Nsample-1);
    point(U[0], _V0, Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  //
  U[0]=_U1;
  for (size_t a=0; a < Nsample-1; ++a)
  {
    U[1]=_V0 + (_V1-_V0)*a/(Nsample-1);
    point(_U1, U[1], Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  //
  U[1]=_V1;
  for (size_t a=0; a < Nsample-1; ++a)
  {
    U[0]=_U1 - (_U1-_U0)*a/(Nsample-1);
    point(U[0], _V1, Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  //
  U[0]=_U0;
  for (size_t a=0; a < Nsample-1; ++a)
  {
    U[1]=_V1 - (_V1-_V0)*a/(Nsample-1);
    point(_U0, U[1], Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  std::vector<E_Int > nids;
  K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(pos3D);
  ::merge(crdA, E_EPSILON, nids);
  
  K_CONNECT::unchanged pred(nids);
  K_CONNECT::IdTool::compress(pos3D, pred);
  K_CONNECT::IdTool::compress(pos2D, pred);
  

  
  return 0;
}

/// Computes the surface point P for the input (u,v) parameters.
void K_OCC::OCCSurface::point(E_Float u, E_Float v, E_Float* P) const{
  gp_Pnt Pt;
  _surface->D0(u,v, Pt);
  P[0]=Pt.X();P[1]=Pt.Y();P[2]=Pt.Z();
}

/// Computes the first U-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DU1(E_Float u, E_Float v, E_Float* P) const{
  gp_Pnt Pt;
  gp_Vec            DU1, DV1;
  
  _surface->D1(u,v, Pt, DU1, DV1);
  
  P[0]=DU1.X();P[1]=DU1.Y();P[2]=DU1.Z();
}

/// Computes the second U-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DU2(E_Float u, E_Float v, E_Float* P) const{
  gp_Pnt            Pt;
  gp_Vec            DU1, DV1, DU2, DV2, DUV;
  
  _surface->D2 (u, v, Pt, DU1, DV1, DU2, DV2, DUV);
  
  P[0]=DU2.X();P[1]=DU2.Y();P[2]=DU2.Z();
}

/// Computes the first V-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DV1(E_Float u, E_Float v, E_Float* P) const{
  gp_Pnt Pt;
  gp_Vec            DU1, DV1;
  _surface->D1(u,v, Pt, DU1, DV1);
  
  P[0]=DV1.X();P[1]=DV1.Y();P[2]=DV1.Z();
}

/// Computes the second U-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DV2(E_Float u, E_Float v, E_Float* P) const{
  gp_Pnt            Pt;
  gp_Vec            DU1, DV1, DU2, DV2, DUV;
  
  _surface->D2 (u, v, Pt, DU1, DV1, DU2, DV2, DUV);
  
  P[0]=DV2.X();P[1]=DV2.Y();P[2]=DV2.Z();
}

/// Computes the first crossed UV-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DUV(E_Float u, E_Float v, E_Float* P) const{
  gp_Pnt            Pt;
  gp_Vec            DU1, DV1, DU2, DV2, DUV;
  
  _surface->D2 (u, v, Pt, DU1, DV1, DU2, DV2, DUV);
  
  P[0]=DUV.X();P[1]=DUV.Y();P[2]=DUV.Z();
}

/// Checks whether input parameters are in the bounds for this surface
bool K_OCC::OCCSurface::in_bounds(E_Float u, E_Float v) const
{
  if (u < _U0) return false;
  if (u > _U1) return false;
  if (v < _V0) return false;
  if (v > _V1) return false;
    
  return true;
}

void K_OCC::OCCSurface::__normalize(E_Float& u, E_Float& v) const 
{
  if (!_normalize_domain) return;
  
  // normalize in [0,1]
  u = (u - _U0) / (_U1 - _U0);
  v = (v - _V0) / (_V1 - _V0);
}


void K_OCC::OCCSurface::__denormalize(E_Float& u, E_Float& v) const 
{

  if (!_normalize_domain) return;
  
  u = u*(_U1-_U0) + _U0; // put in [U0,U1]
  v = v*(_V1-_V0) + _V0; // put in [V0,V1]
  
  u = std::max(_U0, u);
  u = std::min(_U1, u);
  v = std::max(_V0, v);
  v = std::min(_V1, v);
  
  //assert(u <= _U1 && u >= _U0);
  //assert(v <= _V1 && v >= _V0);
  
  // put in [eps,1-eps] : to avoid problems at seam with revolution surface
  
  //u = std::max(u, E_EPSILON);
  //u = std::min(u, _U1 - E_EPSILON);
  //v = std::max(v, E_EPSILON);
  // v = std::min(v, _V1 - E_EPSILON);
}

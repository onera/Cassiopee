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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "OCCSurface.h"
#include "BRep_Tool.hxx"
#include "GeomLib_Tool.hxx"
#include "GeomAPI_ProjectPointOnSurf.hxx"
# include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/KdTree.h"
#include "Nuga/include/Edge.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/IdTool.h"
#include "TopExp_Explorer.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS.hxx"
#include "ShapeAnalysis.hxx"
#include "ShapeAnalysis_Surface.hxx"
#include "StdFail_NotDone.hxx"
#include "BRepAdaptor_Surface.hxx"

///
K_OCC::OCCSurface::OCCSurface(const TopoDS_Face& F, TopTools_IndexedMapOfShape& occ_edges, E_Int pid) 
:_F(F), _surface(BRep_Tool::Surface(F)), _parent(pid), _normalize_domain(true)
{
  // check if the surface is of revolution
  _isUClosed = _isVClosed = false;
  _isRevol = false;
  ShapeAnalysis::GetFaceUVBounds(F, _U0, _U1, _V0, _V1);
  //printf("bounds U=%f %f - V=%f %f\n",_U0,_U1,_V0,_V1);
  // Analyse par OCC
  ShapeAnalysis_Surface o(_surface);
  _isUClosed = o.IsUClosed();
  _isVClosed = o.IsVClosed();
  BRepAdaptor_Surface h(_F);
  _isUPeriodic = h.IsUPeriodic();
  if (_isUPeriodic) _uPeriod = h.UPeriod();
  else _uPeriod = _U1-_U0;
  _isVPeriodic = h.IsVPeriodic();
  if (_isVPeriodic) _vPeriod = h.VPeriod();
  else _vPeriod = _V1-_V0;
  // Comme on dessine uniquement sur la premiere periode:
  if (_uPeriod > _U1-_U0+1.e-2) _isUClosed = false;
  if (_vPeriod > _V1-_V0+1.e-2) _isVClosed = false;
  
  //printf("closed %d %d\n", isUClosed, isVClosed);

  //__get_params_and_type(F, _U0, _U1, _V0, _V1, _isUClosed, _isVClosed);
      
  // Traverse the edges
  TopExp_Explorer edge_expl;
  __traverse_face_edges(F, edge_expl, occ_edges, _edges);
}

///
K_OCC::OCCSurface::OCCSurface(const OCCSurface& rhs)
:_F(rhs._F), _surface(rhs._surface), _parent(rhs._parent), _edges(rhs._edges), _U0(rhs._U0), _U1(rhs._U1), _V0(rhs._V0), _V1(rhs._V1), _normalize_domain(rhs._normalize_domain){}

///
K_OCC::OCCSurface::~OCCSurface() {}

// homothetie par rapport au centre
void K_OCC::OCCSurface::shrink(K_FLD::FloatArray& coord3D, E_Float factor)
{
  E_Float xc=0.,yc=0.,zc=0.;
  E_Float x,y,z;
  E_Int npts = coord3D.getSize();
  for (E_Int i = 0; i < npts; i++)
  {  
    x = coord3D(0,i); y = coord3D(1,i); z = coord3D(2,i);
    xc += x; yc += y; zc += z;
  }
  xc = xc/npts; yc = yc/npts; zc = zc/npts;
  for (E_Int i = 0; i < npts; i++)
  {  
    x = coord3D(0,i); y = coord3D(1,i); z = coord3D(2,i);
    coord3D(0,i) = xc + factor*(x-xc);
    coord3D(1,i) = yc + factor*(y-yc);
    coord3D(2,i) = zc + factor*(z-zc);
  }
}





// parametre le contour (coord3D, connectB) avec les parametres de la surface
// version originale
E_Int
K_OCC::OCCSurface::parameters
(const K_FLD::FloatArray& coord3D, const K_FLD::IntArray& connectB, K_FLD::FloatArray& UVs) const 
{
  UVs.clear();
  
  E_Int err(0), sz(coord3D.cols());
  UVs.resize(2, sz, NUGA::FLOAT_MAX);

  //traverse the edge to check 
  for (E_Int i=0; (i < connectB.cols()) && !err; ++i)
  {
    E_Int Ni = connectB(0,i);
    E_Int Nj = connectB(1,i);
    
    if (UVs(0, Ni) == NUGA::FLOAT_MAX)
    {
      err = parameters(coord3D.col(Ni), UVs(0,Ni), UVs(1,Ni), Ni);
    }
    if (!err && UVs(0, Nj) == NUGA::FLOAT_MAX)
    {
      err = parameters(coord3D.col(Nj), UVs(0,Nj), UVs(1,Nj), Nj);
    }

    if (_isUClosed && ::fabs(UVs(0,Ni) - UVs(0,Nj)) > NUGA::PI) err = 1;
    if (_isVClosed && ::fabs(UVs(1,Ni) - UVs(1,Nj)) > NUGA::PI) err = 1;
  }
  
  for (E_Int k=0; k < UVs.cols(); ++k)
    if (UVs(0,k) == NUGA::FLOAT_MAX)
      UVs(0,k) = UVs(1,k) = 0.;

  return err;
}

// version originale
E_Int
K_OCC::OCCSurface::parameters(const E_Float* pt, E_Float& u, E_Float& v, E_Int index) const 
{
  u=v=-1;
    
  gp_Pnt Point;
  Point.SetCoord(pt[0], pt[1], pt[2]);
  
  E_Int ok = GeomLib_Tool::Parameters(_surface, Point, 1.e+100/*dummytol*/, u, v);
  if (!ok || u==-1 || v==-1) return 1;
  
  // normalize in [0,1]
  __normalize(u,v);
  
  return 0;
}

// Parametre le contour coord3D par rapport a la surface (edges discretises)
E_Int
K_OCC::OCCSurface::parametersSample(const K_FLD::FloatArray&coord3D, K_FLD::FloatArray& UVs)
{
#ifdef DEBUG_CAD_READER
  std::cout << "sampling contour..." << std::endl;
#endif  
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
  
  for (E_Int i = 0; i < sz; ++i)
  {
    const E_Float* Pt = coord3D.col(i);
    Nc = tree.getClosest(Pt);
    
    //init
    N[0] = (Nc > 0) ? Nc-1 : szSmpl-1; //prior
    N[1] = (Nc+1)%szSmpl;              // next
      
    //
    d2[0]=K_MESH::Edge::edgePointMinDistance2<3>(pos3D.col(N[0]), pos3D.col(Nc), Pt, lambda[0]);
    d2[1]=K_MESH::Edge::edgePointMinDistance2<3>(pos3D.col(Nc), pos3D.col(N[1]), Pt, lambda[1]);
    
    // Linear interpolation
    if (d2[0] < d2[1])
      NUGA::sum<2>(1.-lambda[0], pos2D.col(N[0]), lambda[0], pos2D.col(Nc), UVs.col(i));
    else
      NUGA::sum<2>(1.-lambda[1], pos2D.col(Nc), lambda[1], pos2D.col(N[1]), UVs.col(i)); 
  }
  
  return err;
}

E_Int K_OCC::OCCSurface::__sample_contour(E_Int Nsample, K_FLD::FloatArray& pos3D, K_FLD::FloatArray& pos2D)
{
  E_Float V0n(0.), V1n(1.), U0n(0.), U1n(1.);
  
  if (!_normalize_domain)
  {
    V0n = _V0;
    V1n = _V1;
    U0n = _U0;
    U1n = _U1;
  }
  
  E_Float U[2], Pt[3];
  // edge 1
  U[1] = V0n;
  for (E_Int a=0; a < Nsample-1; ++a) // up to Nsample-1 to avoid to create duplicates
  {
    U[0]=U0n + (U1n-U0n)*E_Float(a)/E_Float(Nsample-1);
    point(U[0], V0n, Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  // edge 2
  U[0]=U1n;
  for (E_Int a=0; a < Nsample-1; ++a)
  {
    U[1]=V0n + (V1n-V0n)*E_Float(a)/E_Float(Nsample-1);
    point(U1n, U[1], Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  // edge 3
  U[1]=V1n;
  for (E_Int a=0; a < Nsample-1; ++a)
  {
    U[0]=U1n - (U1n-U0n)*E_Float(a)/E_Float(Nsample-1);
    point(U[0], V1n, Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  // edge 4
  U[0]=U0n;
  for (E_Int a=0; a < Nsample-1; ++a)
  {
    U[1]=V1n - (V1n-V0n)*E_Float(a)/E_Float(Nsample-1);
    point(U0n, U[1], Pt);
    
    pos2D.pushBack(U, U+2);
    pos3D.pushBack(Pt, Pt+3);
  }
  
  std::vector<E_Int > nids;
  K_FLD::ArrayAccessor<K_FLD::FloatArray > crdA(pos3D);
  ::merge(crdA, EPSILON, nids);
  
  K_CONNECT::unchanged pred(nids);
  K_CONNECT::IdTool::compress(pos3D, pred);
  K_CONNECT::IdTool::compress(pos2D, pred);
  
//  _U0 = U0tmp;
//  _U1 = U1tmp;
//  _V0 = V0tmp;
//  _V1 = V1tmp;
  
  return 0;
}

/// Computes the surface point P for the input (u,v) parameters.
void K_OCC::OCCSurface::point(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt Pt;
  
  __denormalize(u,v);
  
  _surface->D0(u, v, Pt);
  P[0]=Pt.X(); P[1]=Pt.Y(); P[2]=Pt.Z();
  
//  std::cout << "coordinated passed to D0 : " << u << "/" << v << std::endl;
//  std::cout << "returned PT : " << P[0] << "/" << P[1] << "/" << P[2] << std::endl;
}

/// Computes the first U-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DU1(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt Pt;
  gp_Vec DU1, DV1;
  
  __denormalize(u,v);
  
  _surface->D1(u, v, Pt, DU1, DV1);
  
  P[0]=DU1.X(); P[1]=DU1.Y(); P[2]=DU1.Z();
  
  if (!_normalize_domain) return;
  
  P[0] *= (_U1 - _U0);
  P[1] *= (_U1 - _U0);
  P[2] *= (_U1 - _U0);
}

/// Computes the second U-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DU2(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt  Pt;
  gp_Vec  DU1, DV1, DU2, DV2, DUV;
  
  __denormalize(u,v);
  
  _surface->D2(u, v, Pt, DU1, DV1, DU2, DV2, DUV);
  
  P[0]=DU2.X(); P[1]=DU2.Y(); P[2]=DU2.Z();
  
  if (!_normalize_domain) return;
  
  P[0] *= (_U1 - _U0)*(_U1 - _U0);
  P[1] *= (_U1 - _U0)*(_U1 - _U0);
  P[2] *= (_U1 - _U0)*(_U1 - _U0);

}

/// Computes the first V-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DV1(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt Pt;
  gp_Vec DU1, DV1;
  
  __denormalize(u,v);
  
  _surface->D1(u, v, Pt, DU1, DV1);
  
  P[0]=DV1.X(); P[1]=DV1.Y(); P[2]=DV1.Z();
  
  if (!_normalize_domain) return;
  
  P[0] *= (_V1 - _V0);
  P[1] *= (_V1 - _V0);
  P[2] *= (_V1 - _V0);
}

/// Computes the second U-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DV2(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt Pt;
  gp_Vec DU1, DV1, DU2, DV2, DUV;
  
  __denormalize(u,v);
  
  _surface->D2(u, v, Pt, DU1, DV1, DU2, DV2, DUV);
  
  P[0]=DV2.X(); P[1]=DV2.Y(); P[2]=DV2.Z();
  
  if (!_normalize_domain) return;
  
  P[0] *= (_V1 - _V0)*(_V1 - _V0);
  P[1] *= (_V1 - _V0)*(_V1 - _V0);
  P[2] *= (_V1 - _V0)*(_V1 - _V0);

}

/// Computes the first crossed UV-derivative at P(u,v) on the surface.
void K_OCC::OCCSurface::DUV(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt Pt;
  gp_Vec DU1, DV1, DU2, DV2, DUV;
  
  __denormalize(u,v);
  
  _surface->D2(u, v, Pt, DU1, DV1, DU2, DV2, DUV);
  
  P[0]=DUV.X(); P[1]=DUV.Y(); P[2]=DUV.Z();
  
  if (!_normalize_domain) return;

  P[0] *= (_U1 - _U0)*(_V1 - _V0);
  P[1] *= (_U1 - _U0)*(_V1 - _V0);
  P[2] *= (_U1 - _U0)*(_V1 - _V0);
}

void K_OCC::OCCSurface::DAUV1(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt Pt;
  gp_Vec DU1, DV1;
  
  __denormalize(u,v);
  
  _surface->D1(u, v, Pt, DU1, DV1);

  P[0]=DU1.X(); P[1]=DU1.Y(); P[2]=DU1.Z();
  P[3]=DV1.X(); P[4]=DV1.Y(); P[5]=DV1.Z();
  
  if (!_normalize_domain) return;

  P[0] *= (_U1 - _U0);
  P[1] *= (_U1 - _U0);
  P[2] *= (_U1 - _U0);
  P[3] *= (_V1 - _V0);
  P[4] *= (_V1 - _V0);
  P[5] *= (_V1 - _V0);
}

void K_OCC::OCCSurface::DAUV2(E_Float u, E_Float v, E_Float* P) const
{
  gp_Pnt Pt;
  gp_Vec DU1, DV1, DU2, DV2, DUV;
  
  __denormalize(u,v);
  
  _surface->D2(u, v, Pt, DU1, DV1, DU2, DV2, DUV);
  
  P[0] =DU1.X(); P[1] =DU1.Y(); P[2] =DU1.Z();
  P[3] =DV1.X(); P[4] =DV1.Y(); P[5] =DV1.Z();
  P[6] =DU2.X(); P[7] =DU2.Y(); P[8] =DU2.Z();
  P[9] =DV2.X(); P[10]=DV2.Y(); P[11]=DV2.Z();
  P[12]=DUV.X(); P[13]=DUV.Y(); P[14]=DUV.Z();

  if (!_normalize_domain) return;

  P[0] *= (_U1 - _U0);
  P[1] *= (_U1 - _U0);
  P[2] *= (_U1 - _U0);
  P[3] *= (_V1 - _V0);
  P[4] *= (_V1 - _V0);
  P[5] *= (_V1 - _V0);

  P[6]  *= (_U1 - _U0)*(_U1 - _U0);
  P[7]  *= (_U1 - _U0)*(_U1 - _U0);
  P[8]  *= (_U1 - _U0)*(_U1 - _U0);
  P[9]  *= (_V1 - _V0)*(_V1 - _V0);
  P[10] *= (_V1 - _V0)*(_V1 - _V0);
  P[11] *= (_V1 - _V0)*(_V1 - _V0);
  P[12] *= (_U1 - _U0)*(_V1 - _V0);
  P[13] *= (_U1 - _U0)*(_V1 - _V0);
  P[14] *= (_U1 - _U0)*(_V1 - _V0);
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
}


// Trouve les edges de la face F
void K_OCC::OCCSurface::__traverse_face_edges(const TopoDS_Face& F, TopExp_Explorer& edge_expl, TopTools_IndexedMapOfShape& occ_edges, std::vector<E_Int>& edges)
{
  for (edge_expl.Init(F, TopAbs_EDGE); edge_expl.More(); edge_expl.Next())
  {
    const TopoDS_Edge& E = TopoDS::Edge(edge_expl.Current());
     
    if (BRep_Tool::Degenerated (E)) continue;

    // Get edge id in the flat list
    E_Int edge_idx = occ_edges.FindIndex(E);
    if (edge_idx == 0) //doesn' exist so add it (due to surface of revolution process)
    {
      occ_edges.Add(E);
      edge_idx = occ_edges.FindIndex(E);
    }

#ifdef DEBUG_CAD_READER
    //assert (E.IsSame (occ_edges(edge_idx)));
#endif

    // Take orientation into account
    if (E.Orientation() != occ_edges(edge_idx).Orientation())
      edge_idx = -edge_idx;
        
    edges.push_back(edge_idx);
  }
}

void K_OCC::OCCSurface::__get_params_and_type
(const TopoDS_Face& F, E_Float& U0, E_Float& U1, E_Float& V0, E_Float& V1, bool& isUClosed, bool& isVClosed)
{
  E_Float paramtol = 1.e-6;
  
  ShapeAnalysis::GetFaceUVBounds(F, U0, U1, V0, V1);
  
  //printf("bounds U=%f %f - V=%f %f\n",U0,U1,V0,V1); 
  
  // min must be a strictly positive epsilon to avoid seam issue for revol surface (and should not hurt for others..)
  U0 = (::fabs(U0) < paramtol ) ? paramtol : U0;
  V0 = (::fabs(V0) < paramtol ) ? paramtol : V0;
  
  E_Float twoPI = 2.*NUGA::PI;
  
  // if max is 2Pi (revolution) , we similarily put a value smaller to avoid the seam issu
  U1 = (::fabs(U1 - twoPI) < paramtol ) ? twoPI - paramtol : U1;
  V1 = (::fabs(V1 - twoPI) < paramtol ) ? twoPI - paramtol : V1;
  
  isUClosed = (U0 == paramtol) & (U1 == (twoPI - paramtol));
  isVClosed = (V0 == paramtol) & (V1 == (twoPI - paramtol));
  
}

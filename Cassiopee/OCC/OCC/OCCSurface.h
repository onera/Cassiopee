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

#ifndef OCCSURFACE_H
#define	OCCSURFACE_H

#include "TopoDS_Face.hxx"
#include "Geom_Surface.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "Nuga/include/DynArray.h"
#include <map>

namespace K_OCC
{

class OCCSurface {
public:
  OCCSurface(const TopoDS_Face& F, TopTools_IndexedMapOfShape& occ_edges, E_Int pid);
  OCCSurface (const OCCSurface& rhs);
  ~OCCSurface();
  
  // Shrink curve
  void shrink(K_FLD::FloatArray& coord3D, E_Float factor);
  
  // Projete coord3D sur la surface
  void project(K_FLD::FloatArray& coord3D) const;

  // Discretise la surface dans coord3D
  void discretize(K_FLD::FloatArray& coord3D, K_FLD::IntArray& connect, E_Int ni, E_Int nj);

  // Order BAR  
  E_Int findNextPoint(K_FLD::IntArray& found, std::vector< std::vector<E_Int> >& node2Elt) const;
  void orderBAR(E_Int npts, K_FLD::FloatArray& coord3D, K_FLD::IntArray& connectB, K_FLD::IntArray& ind, K_FLD::IntArray& start) const;
  void parcoursBAR(K_FLD::FloatArray& pos3D, K_FLD::IntArray& connectB);
  E_Int findNextElement(E_Int e, K_FLD::IntArray& found, K_FLD::IntArray& connectB,
                        std::vector< std::vector<E_Int> >& node2Elt) const;
  void dupBAR(K_FLD::FloatArray& pos3D, K_FLD::IntArray& connectB, K_FLD::IntArray& switcha, std::map< E_Int, E_Int >& mirror);
  E_Int findNonAmbStart(E_Int npts, K_FLD::FloatArray& coord3D) const;

  // version CB
  E_Int parameters2(K_FLD::FloatArray& coord3D, K_FLD::IntArray& connectB, K_FLD::FloatArray& UVs) const ;
  E_Int parameters2(const E_Float* pt, E_Float & u, E_Float& v, E_Int index, E_Float& up, E_Float& vp, E_Float& upp, E_Float& vpp) const ;
  
  // Calcule les parametres UV de coord3D sur la surface
  E_Int parameters(const K_FLD::FloatArray&coord3D, const K_FLD::IntArray& connectB, K_FLD::FloatArray& UVs) const ;
  E_Int parametersSample(const K_FLD::FloatArray&coord3D, K_FLD::FloatArray& UVs);
  
  // Calcule les parametres UV du point P sur la surface
  E_Int parameters(const E_Float* pt, E_Float & u, E_Float& v, E_Int index=-1) const ;
  
  /// Computes the surface point P for the input (u,v) parameters.
  void point(E_Float u, E_Float v, E_Float* P) const;

  /// Computes the first U-derivative at P(u,v) on the surface.
  void DU1(E_Float u, E_Float v, E_Float* P) const;

  /// Computes the second U-derivative at P(u,v) on the surface.
  void DU2(E_Float u, E_Float v, E_Float* P) const;

  /// Computes the first V-derivative at P(u,v) on the surface.
  void DV1(E_Float u, E_Float v, E_Float* P) const;

  /// Computes the second U-derivative at P(u,v) on the surface.
  void DV2(E_Float u, E_Float v, E_Float* P) const;

  /// Computes the first crossed UV-derivative at P(u,v) on the surface.
  void DUV(E_Float u, E_Float v, E_Float* P) const;
  
  /// Compute all derivatives in one go
  void DAUV1(E_Float u, E_Float v, E_Float* P) const;
  void DAUV2(E_Float u, E_Float v, E_Float* P) const;

  /// Checks whether input parameters are in the bounds for this surface
  bool in_bounds(E_Float u, E_Float v) const ;
  
private:
  void __get_params_and_type(const TopoDS_Face& F, E_Float& U0, E_Float& U1, E_Float& V0, E_Float& V1, bool& isUClosed, bool& isVClosed);
  
  E_Int __sample_contour(E_Int Nsample, K_FLD::FloatArray& pos3D, K_FLD::FloatArray& pos2D);
  
  void __traverse_face_edges(const TopoDS_Face& F, TopExp_Explorer& edge_expl, TopTools_IndexedMapOfShape& occ_edges, std::vector<E_Int>& edges);
  
  void __normalize(E_Float& u, E_Float& v) const ; 
  void __denormalize(E_Float & u, E_Float& v) const ;
  
public:

  const TopoDS_Face& _F;
  Handle(Geom_Surface) _surface;
  
  E_Int _parent; //solid
  std::vector<E_Int> _edges;
  
  E_Float _U0, _U1, _V0, _V1;
  bool _isUClosed, _isVClosed, _isRevol;
  bool _isUPeriodic, _isVPeriodic;
  E_Float _uPeriod, _vPeriod;
  bool _normalize_domain;
};

}

#endif	/* OCCSURFACE_H */


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

#ifndef OCCSURFACE_H
#define	OCCSURFACE_H

#include "Fld/DynArray.h"
#include "TopoDS_Face.hxx"
#include "Geom_Surface.hxx"

namespace K_OCC
{

class OCCSurface {
public:
  OCCSurface(const TopoDS_Face&, E_Int id=0);
  OCCSurface (const OCCSurface& rhs);
  ~OCCSurface();
  
  
  ///
  E_Int parameters(const K_FLD::FloatArray&coord3D, K_FLD::FloatArray& UVs);
  E_Int parametersSample(const K_FLD::FloatArray&coord3D, K_FLD::FloatArray& UVs);
  
  ///
  E_Int parameters(const E_Float* pt, E_Float & u, E_Float& v);
  
  
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
  
  /// Checks whether input parameters are in the bounds for this surface
  bool in_bounds(E_Float u, E_Float v) const ;
  
private:
  E_Int __sample_contour(E_Int Nsample, K_FLD::FloatArray& pos3D, K_FLD::FloatArray& pos2D);
  
  void __normalize(E_Float& u, E_Float& v) const ; 
  void __denormalize(E_Float & u, E_Float& v) const ;
  
public:

  const TopoDS_Face& _F;
  Handle(Geom_Surface) _surface;
  
  E_Int _parent; //solid
  std::vector<E_Int> _edges;
  
  E_Float _U0, _U1, _V0, _V1;
  
  bool _normalize_domain;
};

}

#endif	/* OCCSURFACE_H */


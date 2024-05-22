/*    
    Copyright 2013-2024 Onera.

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

#include "occ.h"
#include "TopoDS_Edge.hxx"

E_Boolean areEdgeIdentical(const TopoDS_Edge& E1, const TopoDS_Edge& E2)
{
  BRepAdaptor_Curve C1(E1);
  GeomAdaptor_Curve geomAdap(C1.Curve());
  Standard_Real u0 = geomAdap.FirstParameter();
  Standard_Real u1 = geomAdap.LastParameter();
  
  BRepAdaptor_Curve C2(E2);
  GeomAdaptor_Curve geomAdap(C2.Curve());
  Standard_Real v0 = geomAdap.FirstParameter();
  Standard_Real v1 = geomAdap.LastParameter();
  
  // topological check

  // Get extreme points
  gp_Pnt Pt1A; gp_Pnt Pt1B; gp_Pnt Pt2A; gp_Pnt Pt2B;
  C1.D0(u0, Pt1A);
  C1.D0(u1, Pt1B);
  C2.D0(v0, Pt2A);
  C2.D0(v1, Pt2B);
}

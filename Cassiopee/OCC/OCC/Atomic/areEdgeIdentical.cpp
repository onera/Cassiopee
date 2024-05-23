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
#include "BRepAdaptor_Curve.hxx"
#include "GeomAdaptor_Curve.hxx"

E_Boolean areEdgeIdentical(const TopoDS_Edge& E1, const TopoDS_Edge& E2)
{
  BRepAdaptor_Curve C1(E1);
  GeomAdaptor_Curve geomAdap1(C1.Curve());
  Standard_Real u0 = geomAdap1.FirstParameter();
  Standard_Real u1 = geomAdap1.LastParameter();
  
  BRepAdaptor_Curve C2(E2);
  GeomAdaptor_Curve geomAdap2(C2.Curve());
  Standard_Real v0 = geomAdap2.FirstParameter();
  Standard_Real v1 = geomAdap2.LastParameter();
  
  // topological check
  E_Int sens = 0;
  if (u0 == v0 && u1 == v1) sens = 1;
  else if (u0 == v1 && u1 == v0) sens = 2;
  else return false;

  // Get extreme points
  gp_Pnt Pt1A; gp_Pnt Pt1B; gp_Pnt Pt2A; gp_Pnt Pt2B;
  C1.D0(u0, Pt1A);
  C1.D0(u1, Pt1B);
  if (sens == 1)
  {
    C2.D0(v0, Pt2A);
    C2.D0(v1, Pt2B);
  }
  else
  {
    C2.D0(v0, Pt2B);
    C2.D0(v1, Pt2A);
  }

  // Check distances
  E_Float tol = 1.e-10;
  E_Float dx, dy,dz, dist1, dist2;
  dx = Pt1A.X()-Pt2A.X();
  dy = Pt1A.Y()-Pt2A.Y();
  dz = Pt1A.Z()-Pt2A.Z();
  dist1 = dx*dx+dy*dy+dz*dz;
  dx = Pt1B.X()-Pt2B.X();
  dy = Pt1B.Y()-Pt2B.Y();
  dz = Pt1B.Z()-Pt2B.Z();
  dist2 = dx*dx+dy*dy+dz*dz;

  if (dist1 < tol*tol && dist2 < tol*tol) return true;

  return false;
}

// Python function
PyObject* K_OCC::areEdgeIdentical(PyObject* self, PyObject* args)
{
  return NULL;
}
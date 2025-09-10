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

#include "occ.h"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Ax2.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Shape.hxx>

//=====================================================================
// Add a part cylinder to CAD hook
//=====================================================================
PyObject* K_OCC::addCylinder(PyObject* self, PyObject* args)
{
  PyObject* hook; 
  E_Float xc, yc, zc, xaxis, yaxis, zaxis, R, H;
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ R_ R_, 
    &hook, &xc, &yc, &zc, &xaxis, &yaxis, &zaxis, &R, &H)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  GETMAPEDGES;
  
  // Define the radius, height, and angle of the cylinder
  //Standard_Real angle = 2*M_PI;
  // Create the cylinder with closing faces 
  //TopoDS_Shape part = BRepPrimAPI_MakeCylinder(R, H, angle).Shape();

  // Create only the side of cylinder
  gp_Ax2 axis(gp_Pnt(xc, yc, zc), gp_Dir(xaxis, yaxis, zaxis)); // Axis of the cylinder
  gp_Cylinder cylinder(axis, R); // Radius of 10.0
  Standard_Real uMin = 0.0;
  Standard_Real uMax = 2 * M_PI; // Full circle
  Standard_Real vMin = 0.0;
  Standard_Real vMax = H; // Height of the cylinder
  TopoDS_Face face = BRepBuilderAPI_MakeFace(cylinder, uMin, uMax, vMin, vMax);

  // Rebuild a single compound
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
    
  for (E_Int i = 1; i <= surfaces.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    builder.Add(compound, F);
  }
  for (E_Int i = 1; i <= edges.Extent(); i++)
  {
    TopoDS_Edge E = TopoDS::Edge(edges(i));
    builder.Add(compound, E);
  }
  builder.Add(compound, face);

  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after addPartCylinder: Nb edges=%d\n", se->Extent());
  printf("INFO: after addPartCylinder: Nb faces=%d\n", sf->Extent());
  
  Py_INCREF(Py_None);
  return Py_None;
}

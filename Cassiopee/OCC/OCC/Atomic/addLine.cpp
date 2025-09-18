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

//=====================================================================
// Add a line to CAD hook
//=====================================================================
PyObject* K_OCC::addLine(PyObject* self, PyObject* args)
{
  PyObject* hook; 
  E_Float x1, y1, z1, x2, y2, z2;
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_, &hook, &x1, &y1, &z1, &x2, &y2, &z2)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  GETMAPEDGES;
  
  gp_Pnt p1(x1, y1, z1); // Bottom left
  gp_Pnt p2(x2, y2, z2); // Bottom right

  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(p1, p2);

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
  builder.Add(compound, edge);
  
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
  
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after addLine: Nb edges=%d\n", se->Extent());
  printf("INFO: after addLine: Nb faces=%d\n", sf->Extent());
  
  Py_INCREF(Py_None);
  return Py_None;
}

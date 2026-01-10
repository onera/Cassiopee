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
// Add a box to CAD hook
//=====================================================================
PyObject* K_OCC::addBox(PyObject* self, PyObject* args)
{
  PyObject* hook; 
  E_Float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
  E_Float x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8;
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ TRRR_ TRRR_ TRRR_ TRRR_ TRRR_ TRRR_, 
    &hook, 
    &x1, &y1, &z1, &x2, &y2, &z2, &x3, &y3, &z3, &x4, &y4, &z4,
    &x5, &y5, &z5, &x6, &y6, &z6, &x7, &y7, &z7, &x8, &y8, &z8)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  GETMAPEDGES;

  /* new square */
  gp_Pnt p1(x1, y1, z1); // Bottom left
  gp_Pnt p2(x2, y2, z2); // Bottom right
  gp_Pnt p3(x3, y3, z3); // Top right
  gp_Pnt p4(x4, y4, z4); // Top left
  gp_Pnt p5(x5, y5, z5); // Bottom left
  gp_Pnt p6(x6, y6, z6); // Bottom right
  gp_Pnt p7(x7, y7, z7); // Top right
  gp_Pnt p8(x8, y8, z8); // Top left

  TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(p1, p2);
  TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(p2, p3);
  TopoDS_Edge edge3 = BRepBuilderAPI_MakeEdge(p3, p4);
  TopoDS_Edge edge4 = BRepBuilderAPI_MakeEdge(p4, p1);

  TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge1, edge2, edge3, edge4);
  TopoDS_Face face1 = BRepBuilderAPI_MakeFace(wire);

  edge1 = BRepBuilderAPI_MakeEdge(p5, p8);
  edge2 = BRepBuilderAPI_MakeEdge(p8, p7);
  edge3 = BRepBuilderAPI_MakeEdge(p7, p6);
  edge4 = BRepBuilderAPI_MakeEdge(p6, p5);

  wire = BRepBuilderAPI_MakeWire(edge1, edge2, edge3, edge4);
  TopoDS_Face face2 = BRepBuilderAPI_MakeFace(wire);

  edge1 = BRepBuilderAPI_MakeEdge(p1, p4);
  edge2 = BRepBuilderAPI_MakeEdge(p4, p8);
  edge3 = BRepBuilderAPI_MakeEdge(p8, p5);
  edge4 = BRepBuilderAPI_MakeEdge(p5, p1);

  wire = BRepBuilderAPI_MakeWire(edge1, edge2, edge3, edge4);
  TopoDS_Face face3 = BRepBuilderAPI_MakeFace(wire);

  edge1 = BRepBuilderAPI_MakeEdge(p2, p6);
  edge2 = BRepBuilderAPI_MakeEdge(p6, p7);
  edge3 = BRepBuilderAPI_MakeEdge(p7, p3);
  edge4 = BRepBuilderAPI_MakeEdge(p3, p2);

  wire = BRepBuilderAPI_MakeWire(edge1, edge2, edge3, edge4);
  TopoDS_Face face4 = BRepBuilderAPI_MakeFace(wire);

  edge1 = BRepBuilderAPI_MakeEdge(p4, p3);
  edge2 = BRepBuilderAPI_MakeEdge(p3, p7);
  edge3 = BRepBuilderAPI_MakeEdge(p7, p8);
  edge4 = BRepBuilderAPI_MakeEdge(p8, p4);

  wire = BRepBuilderAPI_MakeWire(edge1, edge2, edge3, edge4);
  TopoDS_Face face5 = BRepBuilderAPI_MakeFace(wire);

  edge1 = BRepBuilderAPI_MakeEdge(p1, p5);
  edge2 = BRepBuilderAPI_MakeEdge(p5, p6);
  edge3 = BRepBuilderAPI_MakeEdge(p6, p2);
  edge4 = BRepBuilderAPI_MakeEdge(p2, p1);

  wire = BRepBuilderAPI_MakeWire(edge1, edge2, edge3, edge4);
  TopoDS_Face face6 = BRepBuilderAPI_MakeFace(wire);

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
  builder.Add(compound, face1);
  builder.Add(compound, face2);
  builder.Add(compound, face3);
  builder.Add(compound, face4);
  builder.Add(compound, face5);
  builder.Add(compound, face6);

  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after addBox: Nb edges=%d\n", se->Extent());
  printf("INFO: after addBox: Nb faces=%d\n", sf->Extent());
  
  Py_INCREF(Py_None);
  return Py_None;
}

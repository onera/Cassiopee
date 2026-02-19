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
// Add a square to CAD hook from P0, width, height, depth
//=====================================================================
PyObject* K_OCC::addSquare(PyObject* self, PyObject* args)
{
  PyObject* hook; 
  E_Float x0, y0, z0, width, height;
  E_Int makeFace;
  char* name;
  if (!PYPARSETUPLE_(args, O_ TRRR_ RR_ I_ S_, &hook, &x0, &y0, &z0, 
    &width, &height, &makeFace, &name)) return NULL;

  GETSHAPE;

  /* new square */
  gp_Pnt p1(x0, y0, z0); // Bottom left
  gp_Pnt p2(x0+width, y0, z0); // Bottom right
  gp_Pnt p3(x0+width, y0+height, z0); // Top right
  gp_Pnt p4(x0, y0+height, z0); // Top left

  TopoDS_Edge edge1 = BRepBuilderAPI_MakeEdge(p1, p2);
  TopoDS_Edge edge2 = BRepBuilderAPI_MakeEdge(p2, p3);
  TopoDS_Edge edge3 = BRepBuilderAPI_MakeEdge(p3, p4);
  TopoDS_Edge edge4 = BRepBuilderAPI_MakeEdge(p4, p1);

  TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge1, edge2, edge3, edge4);
  TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);

#ifdef USEXCAF

  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
  builder.Add(compound, wire);
  if (makeFace == 1) builder.Add(compound, face);
  
  TDocStd_Document* doc = (TDocStd_Document*)packet[5];
  addShape2OCAF(compound, name, *doc);
  TopoDS_Shape* newshp = copyOCAF2TopShape(*doc);
  delete shape;
  SETSHAPE(newshp);
  Py_INCREF(Py_None);
  return Py_None;

#else
  GETMAPSURFACES;
  GETMAPEDGES;

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
  // Add the square face or wire
  if (makeFace == 1) builder.Add(compound, face);
  else builder.Add(compound, wire);
  
  delete shape;
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
  SETSHAPE(newshp);

  printf("INFO: after addSquare: Nb edges=%d\n", se->Extent());
  printf("INFO: after addSquare: Nb faces=%d\n", sf->Extent());
  
  Py_INCREF(Py_None);
  return Py_None;
#endif
}

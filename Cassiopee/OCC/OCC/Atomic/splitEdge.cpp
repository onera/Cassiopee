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
// split edge
#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "TopoDS_Wire.hxx"
#include "BRepCheck_Wire.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "ShapeAnalysis.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "Geom2d_Curve.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "TColGeom_SequenceOfSurface.hxx"
#include "ShapeExtend_CompositeSurface.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "BRep_Builder.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"

//=====================================================================
// Split edge by point or parameter
//=====================================================================
PyObject* K_OCC::splitEdge(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int no; E_Float param; E_Float x, y, z;
  if (!PYPARSETUPLE_(args, O_ I_ R_ TRRR_, &hook, &no, &param, &x, &y, &z)) return NULL;

  GETSHAPE;
  GETMAPEDGES;
  GETMAPSURFACES;

  TopoDS_Edge edge = TopoDS::Edge(edges(no)); // existing edge
  E_Float f, l;
  Handle(Geom_Curve) curve = BRep_Tool::Curve(edge, f, l);

  // Find parameter of the split point
  if (param == -999.)
  {
    gp_Pnt splitPoint(x, y, z);
    GeomAdaptor_Curve adaptor(curve);
    param = GeomAPI_ProjectPointOnCurve(splitPoint, curve).Parameter(1);
  }

  // Build two new edges
  TopoDS_Edge e1 = BRepBuilderAPI_MakeEdge(curve, f, param);
  TopoDS_Edge e2 = BRepBuilderAPI_MakeEdge(curve, param, l);

  // rebuild compound
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
    if (i != no)
    {
      TopoDS_Edge E = TopoDS::Edge(edges(i));
      builder.Add(compound, E);
    }
    else
    {
      builder.Add(compound, e1);
      builder.Add(compound, e2);
    }
  }
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);

  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after split: Nb edges=%d\n", se->Extent());
  printf("INFO: after split: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

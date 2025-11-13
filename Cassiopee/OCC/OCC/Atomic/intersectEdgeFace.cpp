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
// check intersection between edges and faces

#include "occ.h"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "Precision.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepOffsetAPI_ThruSections.hxx"
#include "GeomFill_BSplineCurves.hxx"
#include "Geom_BSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "BRepAlgoAPI_Section.hxx"

//=====================================================================
// Loft
//=====================================================================
PyObject* K_OCC::intersectEdgeFace(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges; PyObject* listFaces; 
  if (!PYPARSETUPLE_(args, OOO_ , &hook, &listEdges, &listFaces)) return NULL;

  GETSHAPE;
  GETMAPEDGES;
  GETMAPSURFACES;

  E_Int nedges = PyList_Size(listEdges);
  E_Int nfaces = PyList_Size(listFaces);
  
  // Get edges and add to compound1
  BRep_Builder builder1;
  TopoDS_Compound compound1;
  builder1.MakeCompound(compound1);
  for (E_Int i = 0; i < nedges; i++)
  {
    PyObject* noO = PyList_GetItem(listEdges, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    builder1.Add(compound1, E);
  }
  
  // Get faces and add to compund2
  BRep_Builder builder2;
  TopoDS_Compound compound2;
  builder2.MakeCompound(compound2);
  for (E_Int i = 0; i < nfaces; i++)
  {
    PyObject* noO = PyList_GetItem(listFaces, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(no));
    builder2.Add(compound2, F);
  }

  BRepAlgoAPI_Section section(compound1, compound2);
  section.ComputePCurveOn1(Standard_True);
  section.Approximation(Standard_True);
  section.Build();
  PyObject* l = PyList_New(0);

  if (section.IsDone())
  {
    TopoDS_Shape result = section.Shape();

    // extract vertices if possible
    TopTools_IndexedMapOfShape v = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(result, TopAbs_VERTEX, v);
    for (E_Int i = 1; i <= v.Extent(); i++)
    {
      const TopoDS_Vertex& vertex = TopoDS::Vertex(v(i));
      gp_Pnt point = BRep_Tool::Pnt(vertex);
      PyObject* val = Py_BuildValue("[d,d,d]", point.X(), point.Y(), point.Z());
      PyList_Append(l, val); Py_DECREF(val);
    }
  }

  // Rebuild the hook
  //delete shape;
  //SETSHAPE(newshp);

  //printf("INFO: after intersect: Nb edges=%d\n", se->Extent());
  //printf("INFO: after intersect: Nb faces=%d\n", sf->Extent());

  return l;
}

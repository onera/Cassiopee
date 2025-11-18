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
// build a gordon surface from a net of edges

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
#include "BRepOffsetAPI_ThruSections.hxx"
#include "GeomFill_BSplineCurves.hxx"
#include "Geom_BSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "occ_gordon.h"

//=====================================================================
// Loft
//=====================================================================
PyObject* K_OCC::addGordonSurface(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listSet1; PyObject* listSet2; 
  if (!PYPARSETUPLE_(args, OOO_ , &hook, &listSet1, &listSet2)) return NULL;

  GETSHAPE;
  GETMAPEDGES;

  E_Int n1 = PyList_Size(listSet1);
  E_Int n2 = PyList_Size(listSet2);
  TopoDS_Shape* newshp = NULL;
  Standard_Real firstParam, lastParam;

  std::vector<Handle (Geom_Curve)> ucurves;
  std::vector<Handle (Geom_Curve)> vcurves;
                                                      
  // Get curves from 1
  for (E_Int i = 0; i < n1; i++)
  {
    PyObject* noO = PyList_GetItem(listSet1, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    //TopoDS_Wire W = BRepBuilderAPI_MakeWire(E);
    Handle(Geom_Curve) baseCurve = BRep_Tool::Curve(E, firstParam, lastParam);
    ucurves.push_back(baseCurve);
  }

  // Get curves from 2
  for (E_Int i = 0; i < n2; i++)
  {
    PyObject* noO = PyList_GetItem(listSet2, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    //TopoDS_Wire W = BRepBuilderAPI_MakeWire(E);
    Handle(Geom_Curve) baseCurve = BRep_Tool::Curve(E, firstParam, lastParam);
    vcurves.push_back(baseCurve);
  }

  Handle(Geom_BSplineSurface) surf;
  surf = occ_gordon::interpolate_curve_network(ucurves, vcurves, 1.e-4);
  TopoDS_Face F = BRepBuilderAPI_MakeFace(surf, Precision::Confusion());
  newshp = new TopoDS_Shape(F);

  // Rebuild the hook
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after addGordonSurface: Nb edges=%d\n", se->Extent());
  printf("INFO: after addGordonSurface: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

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
// loft a list of edges

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
#include "occ_gordon.h"

//=====================================================================
// Loft
//=====================================================================
PyObject* K_OCC::loft(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listProfiles; PyObject* listGuides; 
  if (!PYPARSETUPLE_(args, OOO_ , &hook, &listProfiles, &listGuides)) return NULL;

  GETSHAPE;
  GETMAPEDGES;

  E_Int nprofiles = PyList_Size(listProfiles);
  E_Int nguides = PyList_Size(listGuides);
  TopoDS_Shape* newshp = NULL;

  if (nguides == 0) // loft without guides
  { 
    // Use opencascade BRepOffsetAPI_ThruSections
    BRepOffsetAPI_ThruSections loftBuilder(/*isSolid=*/false, /*is ruled=*/true, /*preserveOrientation=*/1);
    loftBuilder.SetContinuity(GeomAbs_C2);
    //loftBuilder.SetSmoothing(True);

    // Get Profiles and add to builder
    for (E_Int i = 0; i < nprofiles; i++)
    {
      PyObject* noO = PyList_GetItem(listProfiles, i);
      E_Int no = PyInt_AsLong(noO);
      const TopoDS_Edge& E = TopoDS::Edge(edges(no));
      TopoDS_Wire W = BRepBuilderAPI_MakeWire(E);
      loftBuilder.AddWire(W);
    }
    loftBuilder.Build();
    const TopoDS_Shape& loftedSurface = loftBuilder.Shape();
    newshp = new TopoDS_Shape(loftedSurface);
  }
  else // loft with guides
  {
    Standard_Real firstParam, lastParam;
    std::vector<Handle (Geom_Curve)> ucurves;
    std::vector<Handle (Geom_Curve)> vcurves;

    // Get curves from profiles
    for (E_Int i = 0; i < nprofiles; i++)
    {
      PyObject* noO = PyList_GetItem(listProfiles, i);
      E_Int no = PyInt_AsLong(noO);
      const TopoDS_Edge& E = TopoDS::Edge(edges(no));
      Handle(Geom_Curve) baseCurve = BRep_Tool::Curve(E, firstParam, lastParam);
      ucurves.push_back(baseCurve);
    }

    // Get curves from guides
    for (E_Int i = 0; i < nguides; i++)
    {
      PyObject* noO = PyList_GetItem(listGuides, i);
      E_Int no = PyInt_AsLong(noO);
      const TopoDS_Edge& E = TopoDS::Edge(edges(no));
      Handle(Geom_Curve) baseCurve = BRep_Tool::Curve(E, firstParam, lastParam);
      vcurves.push_back(baseCurve);
    }

    Handle(Geom_BSplineSurface) surf;
    surf = occ_gordon::interpolate_curve_network(ucurves, vcurves, 1.e-4);
    TopoDS_Face F = BRepBuilderAPI_MakeFace(surf, Precision::Confusion());
    newshp = new TopoDS_Shape(F);
  }

  // Rebuild the hook
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after loft: Nb edges=%d\n", se->Extent());
  printf("INFO: after loft: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

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
#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepOffsetAPI_ThruSections.hxx"

//=====================================================================
// Loft
//=====================================================================
PyObject* K_OCC::loft(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listProfiles; PyObject* listGuides; 
  if (!PYPARSETUPLE_(args, OOO_ , &hook, &listProfiles, &listGuides)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  //TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];

  BRepOffsetAPI_ThruSections loftBuilder(/*isSolid=*/false, /*is ruled=*/true, /*preserveOrientation=*/1);
  loftBuilder.SetContinuity(GeomAbs_C2);
  //loftBuilder.SetSmoothing(True);

  // Get Profiles and add to builder
  E_Int nprofiles = PyList_Size(listProfiles);
  for (E_Int i = 0; i < nprofiles; i++)
  {
    PyObject* noO = PyList_GetItem(listProfiles, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    TopoDS_Wire W = BRepBuilderAPI_MakeWire(E);
    loftBuilder.AddWire(W);
  }
  loftBuilder.Build();
  TopoDS_Shape loftedSurface = loftBuilder.Shape();

  // Get guides

  // tigl
  //GeomFill_FillingStyle style = GeomFill_CoonsC2Style;
  //style = GeomFill_StretchStyle;
  //SurfMaker.Perform(_myTolerance, _mySameKnotTolerance, style, Standard_True);
  //_result = SurfMaker.Patches();

  // gordon surface

  
  // rebuild    
  TopoDS_Shape* newshp = new TopoDS_Shape(loftedSurface);
  
  // Rebuild the hook
  packet[0] = newshp;
  // Extract surfaces
  TopTools_IndexedMapOfShape* ptr = (TopTools_IndexedMapOfShape*)packet[1];
  delete ptr;
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf);
  packet[1] = sf;

  // Extract edges
  TopTools_IndexedMapOfShape* ptr2 = (TopTools_IndexedMapOfShape*)packet[2];
  delete ptr2;
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se);
  packet[2] = se;
  printf("INFO: after loft: Nb edges=%d\n", se->Extent());
  printf("INFO: after loft: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

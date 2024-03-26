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
// CAD trim a list of face with a list of edges

#include "occ.h"
#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "ShapeUpgrade_UnifySameDomain.hxx"
#include "BRepFeat_SplitShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"

//=====================================================================
// Fix the full shape
//=====================================================================
PyObject* K_OCC::trimFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listOfFaceNo; PyObject* listOfEdgeNo;
  if (!PYPARSETUPLE_(args, OOO_ , &hook, &listOfFaceNo, &listOfEdgeNo)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  
  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];
  TopoDS_Shape* newshp = new TopoDS_Shape();
  
  // Build a wire
  BRepBuilderAPI_MakeWire wireBuilder;
  for (E_Int i = 0; i < PyList_Size(listOfEdgeNo); i++)
  {
    PyObject* edgeNoO = PyList_GetItem(listOfEdgeNo, i);
    E_Int edgeNo = PyLong_AsLong(edgeNoO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(edgeNo));
    wireBuilder.Add(E);
  }
  TopoDS_Wire W = wireBuilder.Wire();

  BRepFeat_SplitShape splitter(*shp);
  // Add faces to splitter
  for (E_Int i = 0; i < PyList_Size(listOfFaceNo); i++)
  {
    PyObject* faceNoO = PyList_GetItem(listOfFaceNo, i);
    E_Int faceNo = PyLong_AsLong(faceNoO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(faceNo));
    splitter.Add(W, F);
  }
  splitter.Build();
  *newshp = splitter.Shape();

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
  printf("INFO: after fix: Nb edges=%d\n", se->Extent());
  printf("INFO: after fix: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

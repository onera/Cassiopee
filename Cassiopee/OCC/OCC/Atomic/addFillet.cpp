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
// Add fillet from edges
#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "BRep_Builder.hxx"
#include "BRepFilletAPI_MakeFillet.hxx"
#include "StdFail_NotDone.hxx"

void getEdgeMap(TopTools_IndexedMapOfShape& oldFaces, TopTools_IndexedMapOfShape& newFaces, PyObject*& faceMap);
void getFaceMap(TopTools_IndexedMapOfShape& oldFaces, TopTools_IndexedMapOfShape& newFaces, PyObject*& faceMap);

//=====================================================================
// Remove some faces and rebuild compound
//=====================================================================
PyObject* K_OCC::addFillet(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges; E_Float radius;
  PyObject* edgeMap; PyObject* faceMap;
  if (!PYPARSETUPLE_(args, OO_ R_ OO_, &hook, &listEdges, &radius, &edgeMap, &faceMap)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  // get top shape
  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];

  // get edges 
  BRepFilletAPI_MakeFillet mkFillet(*shp);
  bool fail = false;
  for (E_Int no = 0; no < PyList_Size(listEdges); no++)
  {
    PyObject* noEdgeO = PyList_GetItem(listEdges, no);
    E_Int noEdge = PyInt_AsLong(noEdgeO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(noEdge));
    try {
      mkFillet.Add(radius, E);
    } catch (StdFail_NotDone& e) { fail = true; }
  }
  if (fail) 
  {
    PyErr_SetString(PyExc_TypeError, "addFillet: failed.");  
    return NULL;
  }

  TopoDS_Shape shc;
  try {
    shc = mkFillet.Shape();
  } catch (StdFail_NotDone& e) { fail = true; }
  if (fail) 
  {
    PyErr_SetString(PyExc_TypeError, "addFillet: failed.");  
    return NULL;
  }
  
  // export
  TopoDS_Shape* newshp = new TopoDS_Shape(shc);
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf);
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se);
  getFaceMap(surfaces, *sf, faceMap);
  getEdgeMap(edges, *se, edgeMap);
  delete shp;
  TopTools_IndexedMapOfShape* ptr = (TopTools_IndexedMapOfShape*)packet[1];
  delete ptr;
  TopTools_IndexedMapOfShape* ptr2 = (TopTools_IndexedMapOfShape*)packet[2];
  delete ptr2;
  packet[0] = newshp;  
  packet[1] = sf;
  packet[2] = se;
  printf("INFO: after addFillet: Nb edges=%d\n", se->Extent());
  printf("INFO: after addFillet: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

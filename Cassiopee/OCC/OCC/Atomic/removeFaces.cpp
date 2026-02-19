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
// Remove CAD faces and rebuild compound
#include "occ.h"

#include "TopoDS.hxx"
#include "TopoDS_Edge.hxx"
#include "BRep_Tool.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "BRep_Builder.hxx"

//=============================================================
// Get new->old for edges
//=============================================================
void getEdgeMap(TopTools_IndexedMapOfShape& oldEdges, TopTools_IndexedMapOfShape& newEdges, PyObject*& edgeMap)
{
  bool found;
  E_Int neold = oldEdges.Extent();
  E_Int nenew = newEdges.Extent();
  PyList_SetSlice(edgeMap, 0, PyList_Size(edgeMap), NULL);
  for (E_Int i = 1; i <= nenew; i++)
  {
    const TopoDS_Edge& E = TopoDS::Edge(newEdges(i));

    found = false;
    for (E_Int j = 1; j <= neold; j++)
    {
      const TopoDS_Edge& EO = TopoDS::Edge(oldEdges(j));
      if (E.IsSame(EO)) 
      {
        //printf("edge identified %d %d\n", i, j);
        PyList_Append(edgeMap, PyLong_FromLong(j));
        found = true;
        break;
      }
    }
    if (not found) PyList_Append(edgeMap, PyLong_FromLong(-1));
    //printf("%d -> %p\n", i, (void*)&E);
  }
}

//===============================================================
// Get new->old for faces
//===============================================================
void getFaceMap(TopTools_IndexedMapOfShape& oldFaces, TopTools_IndexedMapOfShape& newFaces, PyObject*& faceMap)
{
  bool found;
  E_Int neold = oldFaces.Extent();
  E_Int nenew = newFaces.Extent();
  PyList_SetSlice(faceMap, 0, PyList_Size(faceMap), NULL);
  for (E_Int i = 1; i <= nenew; i++)
  {
    const TopoDS_Face& F = TopoDS::Face(newFaces(i));
    found = false;
    for (E_Int j = 1; j <= neold; j++)
    {
      const TopoDS_Face& FO = TopoDS::Face(oldFaces(j));
      if (F.IsSame(FO)) 
      {
        //printf("face identified %d %d\n", i, j);
        PyList_Append(faceMap, PyLong_FromLong(j));
        found = true;
        break;
      }
    }
    if (not found) PyList_Append(faceMap, PyLong_FromLong(-1));
  }
}

//=====================================================================
// Remove some faces and rebuild compound
// output edgeMap et faceMap 
//=====================================================================
PyObject* K_OCC::removeFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces; 
  PyObject* edgeMap; PyObject* faceMap;
  if (!PYPARSETUPLE_(args, OO_ OO_, &hook, &listFaces, &edgeMap, &faceMap)) return NULL;

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
  E_Int nbFaces = surfaces.Extent();

  //Handle(TDF_Data) data = new TDF_Data();
  //TDF_Label label = data->Root();
  //TDataStd_Name::Set(label, "MyFaceTag");
  //TDF_Tool::AddShape(label, F);

  ShapeBuild_ReShape reshaper;
  for (E_Int no = 0; no < PyList_Size(listFaces); no++)
  {
    PyObject* noFaceO = PyList_GetItem(listFaces, no);
    E_Int noFace = PyInt_AsLong(noFaceO);
    if (noFace >= 1 && noFace <= nbFaces)
    {
      const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
      printf("Info: removing face " SF_D_ "\n", noFace);
      reshaper.Remove(F);
    }
    else printf("Warning: removeFaces: invalid face number.\n");
  }
  TopoDS_Shape shc = reshaper.Apply(*shp);

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

  printf("INFO: after removeFaces: Nb edges=%d\n", se->Extent());
  printf("INFO: after removeFaces: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

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
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include "BRepAlgoAPI_Cut.hxx"

//=====================================================================
// Fix the full shape
//=====================================================================
PyObject* K_OCC::trimFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listOfFaceNo1; PyObject* listOfFaceNo2;
  if (!PYPARSETUPLE_(args, OOO_ , &hook, &listOfFaceNo1, &listOfFaceNo2)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  //TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  
  // Build remaining faces list
  std::list<E_Int> pl;
  E_Int nfaces = surfaces.Extent();
  for (E_Int i = 1; i <= nfaces; i++) pl.push_back(i);

  // Build compound from face1
  BRep_Builder builder1;
  TopoDS_Compound compound1;
  builder1.MakeCompound(compound1);
  for (E_Int i = 0; i < PyList_Size(listOfFaceNo1); i++)
  {
    PyObject* faceNoO = PyList_GetItem(listOfFaceNo1, i);
    E_Int faceNo = PyLong_AsLong(faceNoO);
    auto it = std::find(pl.begin(), pl.end(), faceNo);
    if (it != pl.end()) pl.erase(it);
    const TopoDS_Face& F = TopoDS::Face(surfaces(faceNo));
    builder1.Add(compound1, F);
  }

  // Build compound from face2
  BRep_Builder builder2;
  TopoDS_Compound compound2;
  builder1.MakeCompound(compound2);
  for (E_Int i = 0; i < PyList_Size(listOfFaceNo2); i++)
  {
    PyObject* faceNoO = PyList_GetItem(listOfFaceNo2, i);
    E_Int faceNo = PyLong_AsLong(faceNoO);
    auto it = std::find(pl.begin(), pl.end(), faceNo);
    if (it != pl.end()) pl.erase(it);
    const TopoDS_Face& F = TopoDS::Face(surfaces(faceNo));
    builder2.Add(compound2, F);
  }  

  // trim the compound
  TopoDS_Shape trimmedCompound1 = BRepAlgoAPI_Cut(compound1, compound2);
  TopoDS_Shape trimmedCompound2 = BRepAlgoAPI_Cut(compound2, compound1);
  
  // rebuild
  BRep_Builder builder3;
  TopoDS_Compound compound3;  
  builder3.MakeCompound(compound3);
    
  for (auto& i : pl)
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    builder3.Add(compound3, F);
  }

  TopExp_Explorer expl1(trimmedCompound1, TopAbs_FACE);
  while (expl1.More())
  {
    TopoDS_Shape shape = expl1.Current();
    TopoDS_Face face = TopoDS::Face(shape);
    builder3.Add(compound3, face);
    expl1.Next();
  }
  TopExp_Explorer expl2(trimmedCompound2, TopAbs_FACE);
  while (expl2.More()) 
  {
    TopoDS_Shape shape = expl2.Current();
    TopoDS_Face face = TopoDS::Face(shape);
    builder3.Add(compound3, face);
    expl2.Next();
  }
  
  TopoDS_Shape* newshp = new TopoDS_Shape(compound3);
  //TopoDS_Shape* newshp = new TopoDS_Shape(trimmedCompound2);
  

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
  printf("INFO: after trim: Nb edges=%d\n", se->Extent());
  printf("INFO: after trim: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

//===
/*
#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <Geom_Surface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

Handle(Geom_Surface) surface1 = ...; // Your first surface
Handle(Geom_Surface) surface2 = ...; // Your second surface
TopoDS_Face face1 = BRepBuilderAPI_MakeFace(surface1, Precision::Confusion());
TopoDS_Face face2 = BRepBuilderAPI_MakeFace(surface2, Precision::Confusion());

BRepAlgoAPI_Section section(face1, face2);
section.ComputePCurveOn1(Standard_True);
section.Approximation(Standard_True);
section.Build();

if (!section.IsDone()) {
    std::cerr << "Error: Intersection computation failed." << std::endl;
    return;
}

TopoDS_Shape intersection = section.Shape();
TopoDS_Wire wire = BRepBuilderAPI_MakeWire(TopoDS::Edge(intersection));
TopoDS_Face trimmedFace1 = BRepBuilderAPI_MakeFace(surface1, wire, Standard_True);
TopoDS_Face trimmedFace2 = BRepBuilderAPI_MakeFace(surface2, wire, Standard_True);

// on compounds == 
TopoDS_Shape trimmedCompound = BRepAlgoAPI_Cut(compound1, compound2);
*/
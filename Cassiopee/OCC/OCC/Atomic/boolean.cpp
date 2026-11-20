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
// boolean operations

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
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "BRepAlgoAPI_Fuse.hxx"
#include "ShapeUpgrade_UnifySameDomain.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_Common.hxx"

#include "XCAFDoc_DocumentTool.hxx"
#include "TDocStd_Document.hxx"
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"
#include "TDataStd_Name.hxx"

//=====================================================================
// boolean operations
// op=0 (fuse), 1 (cut), 2 (common)
// rev1, rev2: reverse solid
//=====================================================================
PyObject* K_OCC::boolean(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces1; PyObject* listFaces2;
  E_Int op, rev1, rev2;
  if (!PYPARSETUPLE_(args, OOO_ III_ , &hook, &listFaces1, &listFaces2, &op, &rev1, &rev2)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPSURFACES;
  
  // Get compound1
  E_Int nfaces1 = PyList_Size(listFaces1);
  BRep_Builder builder1;
  TopoDS_Compound compound1;
  builder1.MakeCompound(compound1);
  for (E_Int i = 0; i < nfaces1; i++)
  {
    PyObject* noO = PyList_GetItem(listFaces1, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(no));
    builder1.Add(compound1, F);
  }

  BRepBuilderAPI_Sewing sewingTool1;
  sewingTool1.Add(compound1);
  sewingTool1.Perform();
  TopoDS_Shape sewedShape1 = sewingTool1.SewedShape();
  if (sewedShape1.ShapeType() != TopAbs_SHELL)
  {
    PyErr_SetString(PyExc_TypeError,
                    "boolean: surface1 is not a shell.");
    return NULL;
  }

  TopoDS_Shell shell1 = TopoDS::Shell(sewedShape1);
  BRepBuilderAPI_MakeSolid builders1(shell1);
  if (!builders1.IsDone())
  {
    PyErr_SetString(PyExc_TypeError,
                    "boolean: surface1 is not a shell.");
    return NULL;
  }
  TopoDS_Solid solid1 = builders1.Solid();
  if (rev1 == 1) solid1 = TopoDS::Solid(solid1.Reversed());
  
  // Get compound2
  E_Int nfaces2 = PyList_Size(listFaces2);
  BRep_Builder builder2;
  TopoDS_Compound compound2;
  builder2.MakeCompound(compound2);
  for (E_Int i = 0; i < nfaces2; i++)
  {
    PyObject* noO = PyList_GetItem(listFaces2, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(no));
    builder2.Add(compound2, F);
  }

  BRepBuilderAPI_Sewing sewingTool2;
  sewingTool2.Add(compound2);
  sewingTool2.Perform();
  TopoDS_Shape sewedShape2 = sewingTool2.SewedShape();  
  if (sewedShape2.ShapeType() != TopAbs_SHELL)
  {
    PyErr_SetString(PyExc_TypeError,
                    "boolean: surface2 is not a shell.");
    return NULL;
  }

  TopoDS_Shell shell2 = TopoDS::Shell(sewedShape2);
  BRepBuilderAPI_MakeSolid builders2(shell2);
  if (!builders2.IsDone())
  {
    PyErr_SetString(PyExc_TypeError,
                    "boolean: surface2 is not a shell.");
    return NULL;
  }
  TopoDS_Solid solid2 = builders2.Solid();
  if (rev2 == 1) solid2 = TopoDS::Solid(solid2.Reversed());

  printf("before boolean\n"); fflush(stdout);

  TopoDS_Shape result;
  if (op == 0) 
  {
    BRepAlgoAPI_Fuse builder(solid1, solid2);
    if (!builder.IsDone())
    {
      PyErr_SetString(PyExc_TypeError,
                    "boolean: operation fails.");
      return NULL;
    }
    result = builder.Shape();
  }
  else if (op == 1) 
  { 
    BRepAlgoAPI_Cut builder(solid1, solid2);
    if (!builder.IsDone())
    {
      PyErr_SetString(PyExc_TypeError,
                    "boolean: operation fails.");
      return NULL;
    }
    result = builder.Shape();
  }
  else 
  {
    BRepAlgoAPI_Common builder(solid1, solid2);
    if (!builder.IsDone())
    {
      PyErr_SetString(PyExc_TypeError,
                    "boolean: operation fails.");
      return NULL;
    }
    result = builder.Shape();
  }

  //ShapeUpgrade_UnifySameDomain unify(result);
  //unify.Build();
  //TopoDS_Shape unified = unify.Shape();

#ifdef USEXCAF

  GETDOC;
  GETSHAPETOOL;

  // suppress somes faces
  std::map< E_Int, std::vector<E_Int> > label2Faces;
  std::map< E_Int, std::vector<E_Int> > label2Edges;
  getLabel2Edges(*doc, label2Edges);
  getLabel2Faces(*doc, label2Faces);

  for (E_Int i = 0; i < PyList_Size(listFaces1); i++)
  {
    PyObject* noFaceO = PyList_GetItem(listFaces1, i);
    E_Int noFace = PyInt_AsLong(noFaceO);
    for (size_t j = 0; j < label2Faces.size(); j++)
    {
      std::vector< E_Int >& f = label2Faces[j];
      f.erase(std::remove(f.begin(), f.end(), noFace), f.end());
    }
  }
  for (E_Int i = 0; i < PyList_Size(listFaces2); i++)
  {
    PyObject* noFaceO = PyList_GetItem(listFaces2, i);
    E_Int noFace = PyInt_AsLong(noFaceO);
    for (size_t j = 0; j < label2Faces.size(); j++)
    {
      std::vector< E_Int >& f = label2Faces[j];
      f.erase(std::remove(f.begin(), f.end(), noFace), f.end());
    }
  }
  copyTopShape2OCAF(*shape, label2Edges, label2Faces, *doc);

  // Add compound3 as a new shape
  TDF_Label label = shapeTool->AddShape(result);
  TDataStd_Name::Set(label, "boolean");

  // back copy
  TopoDS_Shape* newshp = copyOCAF2TopShape(*doc);
  delete shape;
  SETSHAPE(newshp);
  Py_INCREF(Py_None);
  return Py_None;

#else

  TopoDS_Shape* newshp = new TopoDS_Shape(result);
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after boolean: Nb edges=%d\n", se->Extent());
  printf("INFO: after boolean: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
#endif
}

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
#include "occ.h"

#include "TDF_Label.hxx"
#include "TDF_LabelSequence.hxx"
#include "TDF_Tool.hxx"
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "TDocStd_Document.hxx"
#include "TDataStd_Name.hxx"
#include "TDF_ChildIterator.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Face.hxx"
#include "Geom_Surface.hxx"
#include "BRep_Tool.hxx"

//=====================================================================
// Get face names
// Return [names, [face no]] 
//=====================================================================
PyObject* K_OCC::getFaceNameInOCAF(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TDocStd_Document* doc = (TDocStd_Document*)packet[5];
  if (doc == NULL) 
  {
    PyErr_SetString(PyExc_TypeError, "printOCAF: no OCAF document.");
    return NULL;
  }

  TopoDS_Shape* topShape = (TopoDS_Shape*)packet[0];

  TopTools_IndexedMapOfShape& allFaces = *(TopTools_IndexedMapOfShape*)packet[1];

  TDF_LabelSequence labels;
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
  shapeTool->GetShapes(labels);

  PyObject* out = NULL;
  out = PyList_New(0);

  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);
    Handle(TDataStd_Name) NAME = new TDataStd_Name();
    // Get shape associated with label
    TopoDS_Shape shape = shapeTool->GetShape(label);

    // Face list corresponding to that shape
    /*
    TopTools_IndexedMapOfShape faces = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(shape, TopAbs_FACE, faces);
    for (E_Int k = 1; k <= allFaces.Extent(); k++) // top shape
    {
      TopoDS_Face F1 = TopoDS::Face(allFaces(k));
      auto itr1 = allFaces(k);
      for (E_Int j = 1; j <= faces.Extent(); j++)
      {
        //printf("compound faces = %d\n", j);
        TopoDS_Face F2 = TopoDS::Face(faces(j));
        auto itr2 = faces(j);
        if (F1.IsEqual(F2) == true) printf("id detected=%d %d \n", j, k);
        if (itr1 == itr2) printf("id detected=%d %d\n", j, k);
      }
    }
    */

    /*
    TopExp_Explorer explorer2(shape, TopAbs_FACE);
    while (explorer2.More())
    {
      TopoDS_Face face1 = TopoDS::Face(explorer2.Current());
      
      TopExp_Explorer explorer(*topShape, TopAbs_FACE);
      while (explorer.More())
      {
        TopoDS_Face face2 = TopoDS::Face(explorer.Current());
        if (explorer.Current() == explorer2.Current()) printf("found\n");
        if (face1.IsEqual(face2)) printf("found\n");
        explorer.Next();
      }
      explorer2.Next();
    }
    */

    if (label.FindAttribute(TDataStd_Name::GetID(), NAME)) // retourne tous les attributs de type string
    { 
      TopExp_Explorer explorer2(shape, TopAbs_FACE);
      E_Int n = 0;
      while (explorer2.More())
      {
        TopoDS_Face face1 = TopoDS::Face(explorer2.Current());
        TopExp_Explorer explorer(*topShape, TopAbs_FACE);
        E_Int n1 = 0;
        while (explorer.More()) // topshape
        {
          TopoDS_Face face2 = TopoDS::Face(explorer.Current());
          if (explorer.Current() == explorer2.Current()) printf("found by ptr\n");
          if (face1.IsEqual(face2)) printf("found face\n");
          if (face1 == face2) printf("found face\n");
          if (face1.IsSame(face2)) printf("found face\n");
          Handle(TopoDS_TShape) tface1 = face1.TShape();
          Handle(TopoDS_TShape) tface2 = face2.TShape();
          if (tface1 == tface2) printf("found tface\n");


          Handle(Geom_Surface) surface1 = BRep_Tool::Surface(face1);
          Handle(Geom_Surface) surface2 = BRep_Tool::Surface(face2);
          //if (surface1->IsEqual(surface2)) printf("found geom\n");

          explorer.Next();
          n1++;
        }
        n++;
        explorer2.Next();
      }

      TCollection_ExtendedString labelName = NAME->Get();
      TCollection_AsciiString asciiStr(labelName);
      const char* name = asciiStr.ToCString(); // component name
        
      printf("has string attribute %s = %d\n", name, n);
      
      PyObject* pystring = PyUnicode_FromString(name);
      PyList_Append(out, pystring); Py_DECREF(pystring);
    }
  }
  return out;
}
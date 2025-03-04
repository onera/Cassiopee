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
#include "XCAFDoc_ShapeMapTool.hxx"

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

#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <BRepAlgoAPI_Common.hxx>

#include <TopLoc_Location.hxx>

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
    PyErr_SetString(PyExc_TypeError, "getFaceNameInOCAF: no OCAF document.");
    return NULL;
  }

  TopoDS_Shape* topShape = (TopoDS_Shape*)packet[0];
  
  // Get labels corresponding to shapes
  TDF_LabelSequence labels;
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
  shapeTool->GetShapes(labels);

  PyObject* out = NULL;
  out = PyList_New(0);

  TopTools_IndexedMapOfShape faces1 = TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*topShape, TopAbs_FACE, faces1);
  //TopExp_Explorer explorer1(*topShape, TopAbs_FACE);
  //for (; explorer1.More(); explorer1.Next()) { faces1.Add(explorer1.Current()); }

  printf("Top shape has %d faces \n", faces1.Extent());

  std::vector<E_Int> dejavu(faces1.Extent());
  for (size_t i = 0; i < dejavu.size(); i++) dejavu[i] = 0;

  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);
    Handle(TDataStd_Name) NAME = new TDataStd_Name();

    if (label.FindAttribute(TDataStd_Name::GetID(), NAME)) // retourne tous les attributs de type string
    { 
      TCollection_ExtendedString labelName = NAME->Get();
      std::cout << "Info: label name: " << labelName << std::endl;

      for (size_t i = 0; i < dejavu.size(); i++) dejavu[i] = 0;

      // Get shape associated with label
      TopoDS_Shape shape = shapeTool->GetShape(label);

      TopTools_IndexedMapOfShape faces2 = TopTools_IndexedMapOfShape();
      TopExp::MapShapes(shape, TopAbs_FACE, faces2);
      //TopExp_Explorer explorer2(shape, TopAbs_FACE);
      //for (; explorer2.More(); explorer2.Next()) { faces2.Add(explorer2.Current()); }

      printf("Info: this label shape has %d faces \n", faces2.Extent());

      /*
      E_Int n = 0;
      while (explorer2.More()) // label shape
      {
        TopoDS_Face face2 = TopoDS::Face(explorer2.Current()); // label face
        E_Int n1 = 0;
        while (explorer.More()) // topshape
        {
          TopoDS_Face face1 = TopoDS::Face(explorer.Current()); // topshape face
          
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
      */

      E_Int n = 0;
      std::vector<E_Int> number(faces2.Extent());

      for (E_Int i2 = 0; i2 < faces2.Extent(); i2++) // label face
      {
        n = 1;
        //TopoDS_Face face2 = TopoDS::Face(faces2(i2+1));
        TopoDS_Face face2 = TopoDS::Face(faces2(i2+1).Located(TopLoc_Location()));
        E_Float xmin2, ymin2, zmin2, xmax2, ymax2, zmax2;
        Bnd_Box bbox2;
        BRepBndLib::Add(face2, bbox2);
        bbox2.Get(xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);

        for (E_Int i1 = 0; i1 < faces1.Extent(); i1++) // top face
        {
          if (dejavu[i1] == 0)
          {
            //TopoDS_Face face1 = TopoDS::Face(faces1(i1+1));
            TopoDS_Face face1 = TopoDS::Face(faces1(i1+1).Located(TopLoc_Location()));

            //printf("testing %d vs %d ======================\n", i1, i2);
            //if (face1.IsEqual(face2)) printf("found face\n");
            //if (face1 == face2) printf("found face\n");
            //if (face1.IsSame(face2)) printf("found face\n");
            //if (face1.IsPartner(face2)) { printf("partner faces i1=%d i2=%d\n", i1+1, i2+1); n = i1+1; break; }

            //Handle(TopoDS_TShape) tface1 = face1.TShape();
            //Handle(TopoDS_TShape) tface2 = face2.TShape();
            //if (tface1 == tface2) printf("found tface\n");
          
            //Handle(Geom_Surface) surface1 = BRep_Tool::Surface(face1);
            //Handle(Geom_Surface) surface2 = BRep_Tool::Surface(face2);
            //if (surface1 == surface2) { printf("common surfaces i1=%d i2=%d\n", i1+1, i2+1); n = i1+1; break; }

            //GProp_GProps gprops;
            //BRepGProp::SurfaceProperties(face1, gprops);
            //double area1 = gprops.Mass();
            //BRepGProp::SurfaceProperties(face2, gprops);
            //double area2 = gprops.Mass();
            //printf("area=%g %g\n", area1, area2);

            E_Float xmin1, ymin1, zmin1, xmax1, ymax1, zmax1;
            Bnd_Box bbox1;
            BRepBndLib::Add(face1, bbox1);
            bbox1.Get(xmin1, ymin1, zmin1, xmax1, ymax1, zmax1);
            //if (i1 == 175)
            //{ printf("box1=%g %g %g %g %g %g\n", xmin1, ymin1, zmin1, xmax1, ymax1, zmax1); fflush(stdout); }
            //if ((xmin1 == xmin2) && (ymin1 == ymin2) && (zmin1 == zmin2) &&
            //   (xmax1 == xmax2) && (ymax1 == ymax2) && (zmax1 == zmax2))
            //  { /*printf("common bbox i1=%d i2=%d\n", i1+1, i2+1);*/ n = i1+1; dejavu[i1] = 1; break; }
            if (K_FUNC::fEqual(xmin1, xmin2) && K_FUNC::fEqual(ymin1, ymin2) && K_FUNC::fEqual(zmin1,zmin2) &&
              K_FUNC::fEqual(xmax1, xmax2) && K_FUNC::fEqual(ymax1, ymax2) && K_FUNC::fEqual(zmax1, zmax2))
             { /*printf("common bbox i1=%d i2=%d\n", i1+1, i2+1);*/ n = i1+1; dejavu[i1] = 1; break; }

            //if (i1 == 175) printf("box2=%g %g %g %g %g %g\n", xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);

            //BRepAlgoAPI_Common common(face1, face2);
            //TopoDS_Shape commonShape = common.Shape();
            //if (!commonShape.IsNull()) { printf("common faces i1=%d i2=%d\n", i1+1, i2+1); n = i2+1; break; }
          }
        }
        number[i2] = n;
      }

      // list of faces
      PyObject* fl = PyList_New(0);
      for (size_t i = 0; i < number.size(); i++)
      {
        PyObject* o = PyLong_FromLong(number[i]);
        PyList_Append(fl, o); Py_DECREF(o);
      }
      TCollection_AsciiString asciiStr(labelName);
      const char* name = asciiStr.ToCString(); // component name
      PyObject* pystring = PyUnicode_FromString(name);

      PyList_Append(out, pystring); Py_DECREF(pystring);
      PyList_Append(out, fl); Py_DECREF(fl);
    }
  }
  return out;
}

//====================================================================================
PyObject* K_OCC::getFaceNameInOCAF2(PyObject* self, PyObject* args)
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
    PyErr_SetString(PyExc_TypeError, "getFaceNameInOCAF: no OCAF document.");
    return NULL;
  }

  TopoDS_Shape* topShape = (TopoDS_Shape*)packet[0];
  
  // Get labels corresponding to shapes
  TDF_LabelSequence labels;
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
  shapeTool->GetShapes(labels);

  TopTools_IndexedMapOfShape faces = TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*topShape, TopAbs_FACE, faces);

  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i);
    Handle(TDataStd_Name) name = new TDataStd_Name();
    if (label.FindAttribute(TDataStd_Name::GetID(), name)) // retourne tous les attributs de type string
    { 
      TCollection_ExtendedString labelName = name->Get();
      std::cout << "Info: label name: " << labelName << std::endl;
    }

    Handle(XCAFDoc_ShapeMapTool) shapeMapTool;
    Handle(TDF_Attribute) attr;
    if (label.FindAttribute(XCAFDoc_ShapeMapTool::GetID(), attr))
    {
      //shapeMapTool = Handle(XCAFDoc_ShapeMapTool::DownCast(attr));

      //for (E_Int i = 1; i <= faces.Extent(); i++)
      //{
      //  TopoDS_Face F = TopoDS::Face(faces(i));
      //  bool ret = shapeMapTool->IsSubShape(faces(i));
      //  printf("found face %d = %d\n", i, ret);
      //}
    }
  }

  for (E_Int i = 1; i <= faces.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(faces(i));
    
    // Essai avec FindShape pour retrouver la label d'une face
    /*
    TDF_Label label = shapeTool->FindShape(F);
    Handle(TDataStd_Name) name = new TDataStd_Name();
    if (label.FindAttribute(TDataStd_Name::GetID(), name))
    { 
      TCollection_ExtendedString labelName = name->Get();
      std::cout << "Info: label name: " << labelName << std::endl;
    }*/

    // Essai avec findComponent
    /*
    TDF_LabelSequence labels;
    shapeTool->FindComponent(F, labels);
    for (Standard_Integer i = 1; i <= labels.Length(); i++)
    {
      TDF_Label label = labels.Value(i);
      Handle(TDataStd_Name) name = new TDataStd_Name();
      if (label.FindAttribute(TDataStd_Name::GetID(), name))
      { 
        TCollection_ExtendedString labelName = name->Get();
        std::cout << "Info: label name: " << labelName << std::endl;
      }
    }*/

    // Essai avec FindMainShape
    /*
    TDF_Label label = shapeTool->FindMainShape(F);
    Handle(TDataStd_Name) name = new TDataStd_Name();
    if (label.FindAttribute(TDataStd_Name::GetID(), name))
    { 
      TCollection_ExtendedString labelName = name->Get();
      std::cout << "Info: label name: " << labelName << std::endl;
    }
    */
  }

  // Essai avec isSubShape

  return Py_None;
}

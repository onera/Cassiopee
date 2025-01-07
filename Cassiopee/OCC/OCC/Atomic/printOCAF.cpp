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
#include "occ.h"

#include "STEPCAFControl_Reader.hxx"
#include "TDF_Label.hxx"
#include "TDF_LabelSequence.hxx"
#include "TDF_Tool.hxx"
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "TDocStd_Document.hxx"
#include "TDataStd_Name.hxx"
#include "TDF_ChildIterator.hxx"
#include <iostream>

// recursive iteration through labels
void iterateLabels(const TDF_Label& label) 
{
  // print label entry
  TCollection_AsciiString es;
  TDF_Tool::Entry(label, es); // retourne l'entry de la label (3:0:3 indiquant sa position dans la hierarchie)
  std::cout  <<  "Info: Label entry: " << es.ToCString() << std::endl;

  // Print the label string attributes
  Handle(TDataStd_Name) nameAttr;
  if (label.FindAttribute(TDataStd_Name::GetID(), nameAttr)) 
  {
    TCollection_ExtendedString name = nameAttr->Get();
    std::cout << "Info: Label string: " << name << std::endl;
  }

  // Iterate through child labels
  for (TDF_ChildIterator it(label); it.More(); it.Next()) 
  {
    iterateLabels(it.Value());
  }
}

//=====================================================================
// print all printable OCAF labels
//=====================================================================
PyObject* K_OCC::printOCAF(PyObject* self, PyObject* args)
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

  // parcours recursivement toutes les labels du document
  iterateLabels(doc->Main());

  Py_INCREF(Py_None);
  return Py_None;
}

//=====================================================================
// print shape printable OCAF labels
//=====================================================================
PyObject* K_OCC::printShapeOCAF(PyObject* self, PyObject* args)
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

  //============================================================
  // Retourne une collection de labels correspondants aux shapes
  //============================================================
  TDF_LabelSequence labels;
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc->Main());
  shapeTool->GetFreeShapes(labels);

  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i); // one label of collection

    TCollection_AsciiString es;
    TDF_Tool::Entry(label, es); // retourne l'entry de la label (3:0:3 indiquant sa position dans la hierarchie)
    std::cout  <<  "Info: Label entry: " << es.ToCString()  <<  std::endl;

    Handle(TDataStd_Name) NAME = new TDataStd_Name(); 
    if (label.FindAttribute(TDataStd_Name::GetID(), NAME)) // retourne tous les attributs de type string
    { 
      printf("has string attribute\n"); 
      TCollection_ExtendedString labelName = NAME->Get();
      std::cout << "Info: Label name: " << labelName << std::endl;
    } 
    else 
    { printf("no string attribute\n"); } 
  }

  Py_INCREF(Py_None);
  return Py_None;
}





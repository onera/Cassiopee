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
#include "TDataStd_AsciiString.hxx"
#include "TDF_ChildIterator.hxx"
#include <iostream>

#include "TNaming_UsedShapes.hxx"
#include "TNaming_NamedShape.hxx"
#include "TNaming_Evolution.hxx"
#include "TDF_Reference.hxx"
#include "TDataStd_Expression.hxx"
#include "TDataStd_Comment.hxx"
#include "TDataStd_ExtStringArray.hxx"
#include "TDataStd_ExtStringList.hxx"
#include "TDataStd_Integer.hxx"
#include "TDataStd_IntegerArray.hxx"
#include "TDataStd_IntegerList.hxx"
#include "TDataStd_Real.hxx"
#include "TDataStd_RealArray.hxx"
#include "TDataStd_RealList.hxx"


// recursive iteration through labels
void iterateLabels(const TDF_Label& label) 
{
  // print label entry
  TCollection_AsciiString es;
  TDF_Tool::Entry(label, es); // retourne l'entry de la label (3:0:3 indiquant sa position dans la hierarchie)
  std::cout  <<  "Info: Label entry: " << es.ToCString() << std::endl;

  // Number of attributes
  E_Int n = TDF_Tool::NbAttributes(label);
  std::cout  <<  "Info: number of attached attributes: " << n << std::endl;

  // Name String attributes
  Handle(TDataStd_Name) nameAttr;
  if (label.FindAttribute(TDataStd_Name::GetID(), nameAttr))
  {
    TCollection_ExtendedString name = nameAttr->Get();
    std::cout << ">>>>: Label string: " << name << std::endl;
  }

  // NamedShape attribute
  Handle(TNaming_NamedShape) tnnsAttr;
  if (label.FindAttribute(TNaming_NamedShape::GetID(), tnnsAttr))
  {
    std::cout << ">>>>: TNaming_NamedShape " << std::endl;
    std::cout << ">>>>: TNaming_NamedShape Version " << tnnsAttr->Version() << std::endl;
    TNaming_Evolution evo = tnnsAttr->Evolution();
    if (evo == TNaming_PRIMITIVE)
      std::cout << ">>>>: Evolution PRIMITIVE"<< std::endl;

  }

  // UsedShapes attribute
  Handle(TNaming_UsedShapes) tnusAttr;
  if (label.FindAttribute(TNaming_UsedShapes::GetID(), tnusAttr))
  {
    std::cout << ">>>>: TNaming_UsedShapes " << std::endl;
  }

  // Asciistring attribute
  Handle(TDataStd_AsciiString) tasAttr;
  if (label.FindAttribute(TDataStd_AsciiString::GetID(), tasAttr))
  {
    std::cout << ">>>>: Ascii string " << tasAttr->Get() << std::endl;
  }

  // Reference (parametre?)
  Handle(TDF_Reference) ref;
  if (label.FindAttribute(TDF_Reference::GetID(), ref))
  {
    std::cout << ">>>>: reference detected " << std::endl;
  }

  // Comment
  Handle(TDataStd_Comment) attComment;
  if (label.FindAttribute(TDataStd_Comment::GetID(), attComment))
  {
    std::cout << ">>>>: comment detected " << std::endl;
  }

  // Expression
  Handle(TDataStd_Expression) attExpression;
  if (label.FindAttribute(TDataStd_Expression::GetID(), attExpression))
  {
    std::cout << ">>>>: expression detected " << std::endl;
  }

  // ExtStringArray
  Handle(TDataStd_ExtStringArray) attExtStringArray;
  if (label.FindAttribute(TDataStd_Expression::GetID(), attExtStringArray))
  {
    std::cout << ">>>>: extStringArray detected "  << std::endl;
  }

  // ExtStringList
  Handle(TDataStd_ExtStringList) attExtStringList;
  if (label.FindAttribute(TDataStd_Expression::GetID(), attExtStringList))
  {
    std::cout << ">>>>: extStringList detected "  << std::endl;
  }

  // Integer
  Handle(TDataStd_Integer) attInteger;
  if (label.FindAttribute(TDataStd_Integer::GetID(), attInteger))
  {
    std::cout << ">>>>: integer detected "  << std::endl;
  }

  // IntegerArray
  Handle(TDataStd_IntegerArray) attIntegerArray;
  if (label.FindAttribute(TDataStd_IntegerArray::GetID(), attIntegerArray))
  {
    std::cout << ">>>>: integerArray detected "  << std::endl;
  }

  // IntegerList
  Handle(TDataStd_IntegerList) attIntegerList;
  if (label.FindAttribute(TDataStd_IntegerArray::GetID(), attIntegerList))
  {
    std::cout << ">>>>: integerList detected "  << std::endl;
  }

  // Real
  Handle(TDataStd_Real) attReal;
  if (label.FindAttribute(TDataStd_Real::GetID(), attReal))
  {
    std::cout << ">>>>: real detected "  << std::endl;
  }

  // RealArray
  Handle(TDataStd_RealArray) attRealArray;
  if (label.FindAttribute(TDataStd_RealArray::GetID(), attRealArray))
  {
    std::cout << ">>>>: realArray detected "  << std::endl;
  }

  // RealList
  Handle(TDataStd_RealList) attRealList;
  if (label.FindAttribute(TDataStd_RealArray::GetID(), attRealList))
  {
    std::cout << ">>>>: realList detected "  << std::endl;
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
  shapeTool->GetShapes(labels);

  for (Standard_Integer i = 1; i <= labels.Length(); i++)
  {
    TDF_Label label = labels.Value(i); // one label of collection

    TCollection_AsciiString es;
    TDF_Tool::Entry(label, es); // retourne l'entry de la label (3:0:3 indiquant sa position dans la hierarchie)
    std::cout  <<  "Info: Label entry: " << es.ToCString()  <<  std::endl;

    Handle(TDataStd_Name) NAME = new TDataStd_Name(); 
    if (label.FindAttribute(TDataStd_Name::GetID(), NAME)) // retourne tous les attributs de type string
    { 
      TCollection_ExtendedString labelName = NAME->Get();
      std::cout << "Info: Label name: " << labelName << std::endl;
    } 
    else 
    { /*printf("no string attribute\n");*/ } 
  }

  Py_INCREF(Py_None);
  return Py_None;
}

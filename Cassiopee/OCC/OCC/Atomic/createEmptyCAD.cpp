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
#include "TopoDS_Shape.hxx"
#include "TopoDS_Compound.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"

#include "TDocStd_Document.hxx"
#include "XCAFApp_Application.hxx"
#include "TColStd_SequenceOfAsciiString.hxx"

/*
#include <TDocStd_Application.hxx>
#include <TDocStd_Document.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelSequence.hxx>
#include <TNaming_NamedShape.hxx>
#include <XmlOcafDrivers.hxx>
#include <TopoDS_Shape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
*/

// ============================================================================
/* Convert CAD to OpenCascade hook */
// ============================================================================
PyObject* K_OCC::createEmptyCAD(PyObject* self, PyObject* args)
{
  char* fileName; char* fileFmt;
  if (!PyArg_ParseTuple(args, "ss", &fileName, &fileFmt)) return NULL;

  TopoDS_Shape* shp = new TopoDS_Shape(); // empty shape
  //TopoDS_Shape* shp = new TopoDS_Compound(); // empty shape
  //TopoDS_Shape* shp = NULL;

  // Extract surfaces
  TopTools_IndexedMapOfShape* surfs = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_FACE, *surfs);
  // Extract edges
  TopTools_IndexedMapOfShape* edges = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_EDGE, *edges);
  
  // copy the CAD file name and format
  E_Int l = strlen(fileName);
  char* fileNameC = new char [l+1];
  strcpy(fileNameC, fileName);
  l = strlen(fileFmt);
  char* fileFmtC = new char [l+1];
  strcpy(fileFmtC, fileFmt);

  // Document
  Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication(); // init app at first call
  TDocStd_Document* doc = NULL;
  //doc = new TDocStd_Document("MDTV-Standard");
  doc = new TDocStd_Document("XmlXCAF");
  app->InitDocument(doc);

  // capsule 
  PyObject* hook;
  E_Int sizePacket = 6;
  void** packet = new void* [sizePacket];
  packet[0] = shp; // the top shape
  packet[1] = surfs; // the face map
  packet[2] = edges; // the edge map
  packet[3] = fileNameC; // CAD file name
  packet[4] = fileFmtC; // CAD file format
  packet[5] = doc; // document

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif

  return hook;

}
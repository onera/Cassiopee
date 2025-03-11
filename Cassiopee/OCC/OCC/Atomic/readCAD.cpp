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

// IGES/STEP
#include "IGESControl_Reader.hxx" 
#include "STEPControl_Reader.hxx"

// Data structure
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"

// Document
#include "IGESCAFControl_Reader.hxx"
#include "STEPCAFControl_Reader.hxx"
#include "XCAFDoc_DocumentTool.hxx"
#include "TDocStd_Document.hxx"
#include "XCAFApp_Application.hxx"

#include "XmlDrivers.hxx"
#include "XmlXCAFDrivers.hxx"

#include "TDF_Label.hxx"
#include "TDF_LabelSequence.hxx"
#include "TDF_Tool.hxx"
#include "TDataStd_Name.hxx"

#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"
#include "Standard_Version.hxx"

// ============================================================================
/* Convert CAD to OpenCascade hook */
// ============================================================================
PyObject* K_OCC::readCAD(PyObject* self, PyObject* args)
{
  char* fileName; char* fileFmt;
  if (!PyArg_ParseTuple(args, "ss", &fileName, &fileFmt)) return NULL;

  FILE* ptrFile = fopen(fileName, "r");
  if (ptrFile == NULL)
  {
    PyErr_SetString(PyExc_TypeError, "readCAD: file not found.");
    return NULL;
  }
  fclose(ptrFile);

  TopoDS_Shape* shp = new TopoDS_Shape();
  
  //static Handle(TDocStd_Document) doc2 = new TDocStd_Document("MDTV-Standard");
  static Handle(TDocStd_Document) doc2 = new TDocStd_Document("XmlXCAF"); // static to avoid transcient
  Handle(XCAFApp_Application) app = XCAFApp_Application::GetApplication(); // init app at first call

  if (strcmp(fileFmt, "fmt_iges") == 0)
  {
    // simple read
    //IGESControl_Reader reader;
    //reader.ReadFile(fileName);
    //reader.ClearShapes();
    //reader.TransferRoots();
    //*shp = reader.OneShape();

    // document read
    IGESCAFControl_Reader reader2;
    reader2.ReadFile(fileName);
    reader2.Transfer(doc2);
  }
  else if (strcmp(fileFmt, "fmt_step") == 0)
  {
    // simple read
    //STEPControl_Reader reader;
    //reader.ReadFile(fileName);
    //reader.TransferRoots();
    //*shp = reader.OneShape();
    
    // document read
    STEPCAFControl_Reader reader2;
    reader2.ReadFile(fileName);
    reader2.Transfer(doc2);
  }

  // init app
  app->InitDocument(doc2);

  // get shapes
  Handle(XCAFDoc_ShapeTool) shapeTool = XCAFDoc_DocumentTool::ShapeTool(doc2->Main());
  TDF_LabelSequence labels;

  // get shape
  shapeTool->GetFreeShapes(labels);
  //printf("freeshapes length=%d\n", labels.Length());
  if (labels.Length() == 1)
  {
    TDF_Label label = labels.Value(1);
    *shp = shapeTool->GetShape(label);
  }
  else // build an assembly
  {
    TDF_Label assemblyLabel = shapeTool->NewShape();
    for (Standard_Integer i = 1; i <= labels.Length(); i++)
    {
      TDF_Label label = labels.Value(i);
      shapeTool->AddComponent(assemblyLabel, label, TopLoc_Location());
    }
#if OCC_VERSION_MAJOR < 7
    shapeTool->UpdateAssembly(assemblyLabel);
#else
    shapeTool->UpdateAssemblies();
#endif
    *shp = shapeTool->GetShape(assemblyLabel);
  }

  //XmlXCAFDrivers::DefineFormat(app); // register driver
  //PCDM_StoreStatus status = app->SaveAs(doc2, "toto.xml");
  //if (status != PCDM_SS_OK) printf("can not write document\n");

  //TDocStd_Document* doc = doc2.get();
  TDocStd_Document* doc = doc2.operator->();
  
  // Extract surfaces
  TopTools_IndexedMapOfShape* surfs = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_FACE, *surfs);
  
  // Extract edges
  TopTools_IndexedMapOfShape* edges = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_EDGE, *edges);
  printf("INFO: Nb edges=%d\n", edges->Extent());
  printf("INFO: Nb faces=%d\n", surfs->Extent());
  
  // Extract compounds and solids
  TopTools_IndexedMapOfShape* compounds = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_COMPOUND, *compounds);
  printf("INFO: Nb compounds=%d\n", compounds->Extent());
  TopTools_IndexedMapOfShape* solids = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_SOLID, *solids);
  printf("INFO: Nb solids=%d\n", solids->Extent());

  // copy the CAD file name and format
  E_Int l = strlen(fileName);
  char* fileNameC = new char [l+1];
  strcpy(fileNameC, fileName);
  l = strlen(fileFmt);
  char* fileFmtC = new char [l+1];
  strcpy(fileFmtC, fileFmt);

  // capsule 
  PyObject* hook;
  E_Int sizePacket = 6;
  void** packet = new void* [sizePacket];
  packet[0] = shp; // the top shape
  packet[1] = surfs; // the face map
  packet[2] = edges; // the edge map
  packet[3] = fileNameC; // CAD file name
  packet[4] = fileFmtC; // CAD file format
  packet[5] = doc; // OCAF document

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif

  return hook;
} 

// ============================================================================
// Retourne le nbre de faces dans le hook
// ============================================================================
PyObject* K_OCC::getNbFaces(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;  

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  return Py_BuildValue("l", surfaces.Extent());
}

// ============================================================================
// Retourne le nbre d'edges dans le hook
// ============================================================================
PyObject* K_OCC::getNbEdges(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;  

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  return Py_BuildValue("l", edges.Extent());
}

// ============================================================================
// Retourne le nom du fichier et le format
// ============================================================================
PyObject* K_OCC::getFileAndFormat(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;  

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  return Py_BuildValue("ss", packet[3], packet[4]);
}

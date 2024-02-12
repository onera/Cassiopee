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
// IGES/STEP
#include "IGESControl_Writer.hxx" 
#include "STEPControl_Writer.hxx"

//Data structure
#include "TColStd_HSequenceOfTransient.hxx"
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"

// ============================================================================
/* Write CAD file from OpenCascade hook */
// ============================================================================
PyObject* K_OCC::writeCAD(PyObject* self, PyObject* args)
{
  PyObject* hook;
  char* fileName; char* fileFmt;
  if (!PyArg_ParseTuple(args, "Oss", &hook, &fileName, &fileFmt)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopoDS_Shape* shp = (TopoDS_Shape*) packet[0];
  
  if (strcmp(fileFmt, "fmt_iges") == 0)
  {
    IGESControl_Writer writer;
    writer.AddShape(*shp);
    writer.Write(fileName);
  }
  else if (strcmp(fileFmt, "fmt_step") == 0)
  {
    // Read the file
    STEPControl_Writer writer;
    writer.Transfer(*shp, STEPControl_AsIs);
    writer.Write(fileName);
  }  

  Py_INCREF(Py_None);
  return Py_None;
} 
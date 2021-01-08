/*    
    Copyright 2013-2020 Onera.

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

//Data structure
#include "TColStd_HSequenceOfTransient.hxx"
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"

// ============================================================================
/* Convert CAD to OpenCascade hook */
// ============================================================================
PyObject* K_OCC::readCAD(PyObject* self, PyObject* args)
{
  char* fileName; char* fileFmt;
  if (!PyArg_ParseTuple(args, "ss", &fileName, &fileFmt)) return NULL;

  TopoDS_Shape* shp = new TopoDS_Shape();
  
  if (strcmp(fileFmt, "fmt_iges") == 0)
  {
    IGESControl_Reader reader;
    reader.ReadFile(fileName);
      
    // Transfer CAD faces (only) into a OCC list
    Handle(TColStd_HSequenceOfTransient) occ_list = reader.GiveList("iges-faces");
    
    //if (occ_list.IsNull()) return 1;

    //Standard_Integer nb_transfered_faces = reader.TransferList(occ_list);
    Standard_Integer nb_cad_faces = occ_list->Length();
    printf("IGES Faces: %d\n", nb_cad_faces);
    *shp = reader.OneShape();
  }
  else if (strcmp(fileFmt, "fmt_step") == 0)
  {
    // Read the file
    STEPControl_Reader reader;
    reader.ReadFile(fileName);
    *shp = reader.OneShape();
  }
  
  // Extract surfaces
  TopTools_IndexedMapOfShape* surfs = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_FACE, *surfs);
  
  // Extract edges
  TopTools_IndexedMapOfShape* edges = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*shp, TopAbs_EDGE, *edges);
  
  // capsule 
  PyObject* hook;
  E_Int sizePacket = 3;
  void** packet = new void* [sizePacket];
  packet[0] = shp;
  packet[1] = surfs;
  packet[2] = edges;
  
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif

  return hook;
} 



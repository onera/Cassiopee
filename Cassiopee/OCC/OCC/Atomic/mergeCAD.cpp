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
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Builder.hxx"

// ============================================================================
/* Merge two CAD hooks in a single hook 
   Caller must dealloc input hooks */
// ============================================================================
PyObject* K_OCC::mergeCAD(PyObject* self, PyObject* args)
{
  PyObject* hook1; PyObject* hook2;
  if (!PYPARSETUPLE_(args, OO_, &hook1, &hook2)) return NULL;

  void** packet1 = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet1 = (void**) PyCObject_AsVoidPtr(hook1);
#else
  packet1 = (void**) PyCapsule_GetPointer(hook1, NULL);
#endif
  //TopoDS_Shape* shp1 = (TopoDS_Shape*) packet1[0];
  TopTools_IndexedMapOfShape& surfaces1 = *(TopTools_IndexedMapOfShape*)packet1[1];

  void** packet2 = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet2 = (void**) PyCObject_AsVoidPtr(hook2);
#else
  packet2 = (void**) PyCapsule_GetPointer(hook2, NULL);
#endif
  //TopoDS_Shape* shp2 = (TopoDS_Shape*) packet2[0];
  TopTools_IndexedMapOfShape& surfaces2 = *(TopTools_IndexedMapOfShape*)packet2[1];

  // Rebuild a single compound
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
    
  for (E_Int i = 1; i <= surfaces1.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces1(i));
    builder.Add(compound, F);
  }

  for (E_Int i = 1; i <= surfaces2.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces2(i));
    builder.Add(compound, F);
  }
  
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);

  // capsule 
  PyObject* hook;
  E_Int sizePacket = 5;
  void** packet = new void* [sizePacket];
  packet[0] = newshp;

  // Extract surfaces
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf);
  packet[1] = sf;

  // Extract edges
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se);
  packet[2] = se;
  printf("INFO: after merge: Nb edges=%d\n", se->Extent());
  printf("INFO: after merge: Nb faces=%d\n", sf->Extent());
  
  // copy filenames
  char* fileName = (char*)packet2[3];
  E_Int l = strlen(fileName);
  char* fileNameC = new char [l+1];
  strcpy(fileNameC, fileName);
  packet[3] = fileNameC;
  char* fileFmt = (char*)packet2[4];
  l = strlen(fileFmt);
  char* fileFmtC = new char [l+1];
  strcpy(fileFmtC, fileFmt);
  packet[4] = fileFmtC;

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif

  return hook;
} 
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

#include "occ.h"
#include "TopoDS.hxx"
#include "TopoDS_Edge.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Builder.hxx"
#include "TDocStd_Document.hxx"

// ============================================================================
/* Merge two CAD hooks in a single hook 
   Caller must eventually dealloc input hooks */
// ============================================================================
PyObject* K_OCC::mergeCAD(PyObject* self, PyObject* args)
{
  PyObject* listHooks;
  if (!PYPARSETUPLE_(args, O_, &listHooks)) return NULL;

  // Rebuild a single compound
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);

  E_Int size = PyList_Size(listHooks);
  for (E_Int i = 0; i < size; i++)
  {
    PyObject* hook = PyList_GetItem(listHooks, i);
    GETPACKET;
    GETMAPSURFACES;
    GETMAPEDGES;
    
    for (E_Int i = 1; i <= surfaces.Extent(); i++)
    {
      TopoDS_Face F = TopoDS::Face(surfaces(i));
      builder.Add(compound, F);
    }
    for (E_Int i = 1; i <= edges.Extent(); i++)
    {
      TopoDS_Edge E = TopoDS::Edge(edges(i));
      builder.Add(compound, E);
    }
  }
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);

  // capsule
  CREATEHOOK;
  packet[0] = newshp;

  SETMAPEDGES;
  SETMAPSURFACES;
  printf("INFO: after merge: Nb edges=%d\n", se->Extent());
  printf("INFO: after merge: Nb faces=%d\n", sf->Extent());
  
  // copy filenames
  PyObject* hook2 = PyList_GetItem(listHooks, 0);
  void** packet2 = (void**) PyCapsule_GetPointer(hook2, NULL);  

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
  TDocStd_Document* doc = (TDocStd_Document*)packet2[5]; // todo: must merge document
  packet[5] = doc;

  return hook;
} 
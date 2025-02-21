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

#include <ShapeUpgrade_UnifySameDomain.hxx>
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "BRep_Builder.hxx"
#include "TopoDS.hxx"
#include "TopExp.hxx"
#include <TopoDS_Shape.hxx>
#include <BRepBuilderAPI_Sewing.hxx>

// ============================================================================
/* Merge a list of faces in a single face */
// ============================================================================
PyObject* K_OCC::mergeFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces;
  if (!PYPARSETUPLE_(args, OO_, &hook, &listFaces)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];

  E_Int nfaces = PyList_Size(listFaces);

  TopoDS_Shape* shp = NULL;
  if (nfaces == 0)
  {
    shp = (TopoDS_Shape*)packet[0];
  }
  else
  {
    // Build compounds made of faces
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
      
    for (E_Int no = 0; no < nfaces; no++)
    {
      PyObject* noFaceO = PyList_GetItem(listFaces, no);
      E_Int noFace = PyInt_AsLong(noFaceO);
      const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
      builder.Add(compound, F);
    }
    BRepBuilderAPI_Sewing sewer(1.e-6);
    sewer.Add(compound);
    sewer.Perform();
    shp = new TopoDS_Shape(sewer.SewedShape());
  }

  // Unify the faces
  ShapeUpgrade_UnifySameDomain unifier(*shp, true, true, true);
  unifier.Build();
  TopoDS_Shape unifiedShape = unifier.Shape();
  
  TopoDS_Shape* newshp = NULL;

  if (nfaces == 0)
  {
    newshp = new TopoDS_Shape(unifiedShape);
  }
  else if (nfaces > 0)
  {
    delete shp;

    // rebuild compound
    E_Int* tag = new E_Int [surfaces.Extent()];
    for (E_Int i = 0; i < surfaces.Extent(); i++) tag[i] = 1;
    for (E_Int i = 0; i < nfaces; i++)
    {
      PyObject* noFaceO = PyList_GetItem(listFaces, i);
      E_Int noFace = PyInt_AsLong(noFaceO);
      tag[noFace-1] = 0;
    }

    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
 
    for (E_Int i = 0; i < surfaces.Extent(); i++)
    {
      if (tag[i] == 1) builder.Add(compound, surfaces(i+1));
    }
    delete [] tag;

    TopTools_IndexedMapOfShape sf = TopTools_IndexedMapOfShape();
    TopExp::MapShapes(unifiedShape, TopAbs_FACE, sf);
    printf("number of faces in unified=%d\n", sf.Extent());
    for (E_Int i = 0; i < sf.Extent(); i++) builder.Add(compound, sf(i+1));

    BRepBuilderAPI_Sewing sewer(1.e-6);
    sewer.Add(compound);
    sewer.Perform();
    newshp = new TopoDS_Shape(sewer.SewedShape());
  }
  else
  {
    // DBX
    newshp = new TopoDS_Shape(*shp);
  }

  // export
  packet[0] = newshp;
  // Extract surfaces
  TopTools_IndexedMapOfShape* ptr = (TopTools_IndexedMapOfShape*)packet[1];
  delete ptr;
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf);
  packet[1] = sf;
  
  TopTools_IndexedMapOfShape* ptr2 = (TopTools_IndexedMapOfShape*)packet[2];
  delete ptr2;
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se);
  packet[2] = se;
  
  printf("INFO: after mergeFaces: Nb edges=%d\n", se->Extent());
  printf("INFO: after mergeFaces: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
  
}
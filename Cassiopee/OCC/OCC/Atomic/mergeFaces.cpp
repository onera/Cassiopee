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

  GETSHAPE;
  GETMAPSURFACES;

  TopoDS_Shape* shp = NULL;
  if (listFaces == Py_None)
  {
    shp = (TopoDS_Shape*)packet[0];
  }
  else
  {
    E_Int nfaces = PyList_Size(listFaces);
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

  if (listFaces == Py_None)
  {
    newshp = new TopoDS_Shape(unifiedShape);
  }
  else
  {
    delete shp;
    E_Int nfaces = PyList_Size(listFaces);
    
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

  delete shape;
  SETSHAPE(newshp);
  
  printf("INFO: after mergeFaces: Nb edges=%d\n", se->Extent());
  printf("INFO: after mergeFaces: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
  
}
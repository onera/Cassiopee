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
// extract faces
#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "ShapeAnalysis.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeUpgrade_FaceDivide.hxx"
#include "ShapeUpgrade_ShapeDivideArea.hxx"
#include "ShapeUpgrade_ShapeDivideClosed.hxx"
#include "ShapeUpgrade_ClosedFaceDivide.hxx"
#include "ShapeUpgrade_SplitSurfaceArea.hxx"
#include "TColGeom_SequenceOfSurface.hxx"
#include "ShapeExtend_CompositeSurface.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "BRep_Builder.hxx"
#include "ShapeUpgrade_ShapeDivideClosedEdges.hxx"
#include "BRepBuilderAPI_Sewing.hxx"

//=====================================================================
// extract faces
//=====================================================================
PyObject* K_OCC::extractFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces;
  if (!PYPARSETUPLE_(args, O_ O_, &hook, &listFaces)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPSURFACES;
  
  TopoDS_Shape* newshp = NULL;
  E_Int nfaces = PyList_Size(listFaces);
  // Build a compound
  BRep_Builder builder;
  TopoDS_Compound shc;
  builder.MakeCompound(shc);
    
  for (E_Int no = 0; no < nfaces; no++)
  {
    PyObject* noFaceO = PyList_GetItem(listFaces, no);
    E_Int noFace = PyInt_AsLong(noFaceO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
    builder.Add(shc, F);
  }

  BRepBuilderAPI_Sewing sewer(1.e-6);
  sewer.Add(shc);
  sewer.Perform();
  newshp = new TopoDS_Shape(sewer.SewedShape());

  // export
  delete shape;
  SETSHAPE(newshp);
  
  printf("INFO: after extractFaces: Nb edges=%d\n", se->Extent());
  printf("INFO: after extractFaces: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

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
// CAD split for parallel
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
// Sew faces removing extra edges
//=====================================================================
PyObject* K_OCC::sewing(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float tol;
  if (!PYPARSETUPLE_(args, O_ R_, &hook, &tol)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  //TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  //TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];

  // try on all shape
  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];
  
  
  printf("sew top shape\n");

  const Standard_Real tolerance = tol;
  BRepBuilderAPI_Sewing sewer(tolerance);
  sewer.Add(*shp);
  sewer.Perform();

  TopoDS_Shape* newshp = new TopoDS_Shape(sewer.SewedShape());

  packet[0] = newshp;
  // Extract surfaces
  TopTools_IndexedMapOfShape* ptr = (TopTools_IndexedMapOfShape*)packet[1];
  delete ptr;
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf);
  packet[1] = sf;

  // Extract edges
  TopTools_IndexedMapOfShape* ptr2 = (TopTools_IndexedMapOfShape*)packet[2];
  delete ptr2;
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se);
  packet[2] = se;
  printf("INFO: after sewing: Nb edges=%d\n", se->Extent());
  printf("INFO: after sewing: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

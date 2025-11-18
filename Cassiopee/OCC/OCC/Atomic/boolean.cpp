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
// boolean operations

#include "occ.h"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "Precision.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "BRepAlgoAPI_Fuse.hxx"
#include "ShapeUpgrade_UnifySameDomain.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_Common.hxx"

//=====================================================================
// boolean operations
// op=0 (fuse), 1 (cut), 2 (common)
// rev1, rev2: reverse solid
//=====================================================================
PyObject* K_OCC::boolean(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces1; PyObject* listFaces2;
  E_Int op, rev1, rev2;
  if (!PYPARSETUPLE_(args, OOO_ III_ , &hook, &listFaces1, &listFaces2, &op, &rev1, &rev2)) return NULL;

  GETSHAPE;
  GETMAPEDGES;
  GETMAPSURFACES;

  E_Int nfaces1 = PyList_Size(listFaces1);
  E_Int nfaces2 = PyList_Size(listFaces2);
  
  // Get compound1
  BRep_Builder builder1;
  TopoDS_Compound compound1;
  builder1.MakeCompound(compound1);
  for (E_Int i = 0; i < nfaces1; i++)
  {
    PyObject* noO = PyList_GetItem(listFaces1, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(no));
    builder1.Add(compound1, F);
  }

  BRepBuilderAPI_Sewing sewingTool1;
  sewingTool1.Add(compound1);
  sewingTool1.Perform();
  TopoDS_Shape sewedShape1 = sewingTool1.SewedShape();
  TopoDS_Shell shell1 = TopoDS::Shell(sewedShape1);
  TopoDS_Solid solid1 = BRepBuilderAPI_MakeSolid(shell1);
  if (rev1 == 1) solid1 = TopoDS::Solid(solid1.Reversed());
  
  // Get compound2
  BRep_Builder builder2;
  TopoDS_Compound compound2;
  builder2.MakeCompound(compound2);
  for (E_Int i = 0; i < nfaces2; i++)
  {
    PyObject* noO = PyList_GetItem(listFaces2, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Face& F = TopoDS::Face(surfaces(no));
    builder2.Add(compound2, F);
  }

  BRepBuilderAPI_Sewing sewingTool2;
  sewingTool2.Add(compound2);
  sewingTool2.Perform();
  TopoDS_Shape sewedShape2 = sewingTool2.SewedShape();
  TopoDS_Shell shell2 = TopoDS::Shell(sewedShape2);
  TopoDS_Solid solid2 = BRepBuilderAPI_MakeSolid(shell2);
  if (rev2 == 1) solid2 = TopoDS::Solid(solid2.Reversed());

  TopoDS_Shape result;
  if (op == 0) result = BRepAlgoAPI_Fuse(solid1, solid2);
  else if (op == 1) result = BRepAlgoAPI_Cut(solid1, solid2);
  else result = BRepAlgoAPI_Common(solid1, solid2);

  //ShapeUpgrade_UnifySameDomain unify(result);
  //unify.Build();
  //TopoDS_Shape unified = unify.Shape();

  TopoDS_Shape* newshp = new TopoDS_Shape(result);
  
  // Rebuild the hook
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after boolUnion: Nb edges=%d\n", se->Extent());
  printf("INFO: after boolUnion: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

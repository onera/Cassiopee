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
// boolean operations on wires (in nearly the same plane)

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
#include "ShapeUpgrade_UnifySameDomain.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepAlgoAPI_Common.hxx"

#include "XCAFDoc_DocumentTool.hxx"
#include "TDocStd_Document.hxx"
#include "XCAFDoc_ShapeTool.hxx"
#include "XCAFDoc_ShapeMapTool.hxx"
#include "TDataStd_Name.hxx"

#include "BOPAlgo_PaveFiller.hxx"
#include "BOPAlgo_Builder.hxx"
#include "BRepAlgoAPI_Fuse.hxx"
#include "BRepTools.hxx"
#include "StdFail_NotDone.hxx"

//=====================================================================
// boolean operations on wires
// op=0 (fuse), 1 (cut), 2 (common)
//=====================================================================
PyObject* K_OCC::booleanEdges(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges1; PyObject* listEdges2;
  E_Int op;
  if (!PYPARSETUPLE_(args, OOO_ I_ , &hook, &listEdges1, &listEdges2, &op)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPEDGES;

  // Get wire1
  E_Int nedges1 = PyList_Size(listEdges1);
  BRepBuilderAPI_MakeWire builder1;
  for (E_Int i = 0; i < nedges1; i++)
  {
    PyObject* noO = PyList_GetItem(listEdges1, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    builder1.Add(E);
  }
  TopoDS_Wire w1 = builder1.Wire();

  // Get wire2
  E_Int nedges2 = PyList_Size(listEdges2);
  BRepBuilderAPI_MakeWire builder2;
  for (E_Int i = 0; i < nedges2; i++)
  {
    PyObject* noO = PyList_GetItem(listEdges2, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    builder2.Add(E);
  }
  TopoDS_Wire w2 = builder2.Wire();

  // Compute intersections and split edges
  //BOPAlgo_PaveFiller paveFiller;
  //paveFiller.AddArgument(w1);
  //paveFiller.AddArgument(w2);
  //paveFiller.Perform();
  //if (paveFiller.HasErrors()) 
  //{
  //  PyErr_SetString(PyExc_TypeError, "boolean: intersection/split of wire fails.");  
  //  return NULL;
  //}

  // build faces (may fail if edges are not closed)
  TopoDS_Face f1 = BRepBuilderAPI_MakeFace(w1);
  TopoDS_Face f2 = BRepBuilderAPI_MakeFace(w2);
  
  // bool on faces
  TopoDS_Shape result;
  if (op == 0) result = BRepAlgoAPI_Fuse(f1, f2);
  else if (op == 1) result = BRepAlgoAPI_Cut(f1, f2);
  else result = BRepAlgoAPI_Common(f1, f2);

  // sew result
  BRepBuilderAPI_Sewing sewingTool;
  sewingTool.Add(result);
  sewingTool.Perform();
  TopoDS_Shape result2 = sewingTool.SewedShape();

  // try to unify faces in one face
  TopoDS_Shape result3;
  ShapeUpgrade_UnifySameDomain usd(result2, 1, 1, 1);
  usd.AllowInternalEdges(Standard_True);
  usd.SetLinearTolerance(1.e-6);
  usd.SetAngularTolerance(1.*K_CONST::E_PI/180.);
  try 
  {
    usd.Build();
    result3 = usd.Shape();
  } catch (StdFail_NotDone& e) 
  {
    result3 = result2;
  }
    
  // Build new shape from outer wire
  BRep_Builder builderc;
  TopoDS_Compound compound;
  builderc.MakeCompound(compound);
  for (TopExp_Explorer exp(result3, TopAbs_FACE); exp.More(); exp.Next()) 
  {
    TopoDS_Face f = TopoDS::Face(exp.Current());
    TopoDS_Wire outer = BRepTools::OuterWire(f);
    builderc.Add(compound, outer);
  }

  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after boolean: Nb edges=%d\n", se->Extent());
  printf("INFO: after boolean: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

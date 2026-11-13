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
// utilities

#include "occ.h"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "BRepLib.hxx"
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_Sewing.hxx"
#include "BRepCheck_Analyzer.hxx"
#include "TopExp_Explorer.hxx"

//=================================================================
// Build a compound from list of faces
//=================================================================
void K_OCC::buildfaceCompound(TopoDS_Shape& shape, PyObject* listFaces, TopoDS_Compound& compound)
{
  TopTools_IndexedMapOfShape sf;
  TopExp::MapShapes(shape, TopAbs_FACE, sf);

  if (listFaces != Py_None)
  {
    E_Int nfaces = PyList_Size(listFaces);
    BRep_Builder builder;
    builder.MakeCompound(compound);
    for (E_Int i = 0; i < nfaces; i++)
    {
      PyObject* noO = PyList_GetItem(listFaces, i);
      E_Int no = PyInt_AsLong(noO);
      const TopoDS_Face& F = TopoDS::Face(sf(no));
      builder.Add(compound, F);
    }
  }
  else 
  {
    BRep_Builder builder;
    builder.MakeCompound(compound);
    builder.Add(compound, shape);
  }
}

//=================================================================
// Build a compound from faces not in list
//=================================================================
void K_OCC::buildOutfaceCompound(TopoDS_Shape& shape, PyObject* listFaces, TopoDS_Compound& compound)
{
  if (listFaces == Py_None) return;
  TopTools_IndexedMapOfShape sf;
  TopExp::MapShapes(shape, TopAbs_FACE, sf);
  E_Int nfacesTot = sf.Extent();
  std::vector<E_Int> tag(nfacesTot);
  for (E_Int i = 0; i < nfacesTot; i++) tag[i] = 0;

  E_Int nfaces = PyList_Size(listFaces);
  for (E_Int i = 0; i < nfaces; i++)
  {
    PyObject* noO = PyList_GetItem(listFaces, i);
    E_Int no = PyInt_AsLong(noO);
    tag[no-1] = 1; 
  }
  
  BRep_Builder builder;
  builder.MakeCompound(compound);  
  for (E_Int i = 0; i < nfacesTot; i++)
  {
    if (tag[i] == 0)
    {
      const TopoDS_Face& F = TopoDS::Face(sf(i+1));
      builder.Add(compound, F);
    }
  }
}

//=================================================================
// sew a shape
//=================================================================
void K_OCC::sewShape(TopoDS_Shape& shape, E_Float tol)
{
  BRepBuilderAPI_Sewing sewingTool;
  sewingTool.Add(shape);
  sewingTool.Perform();
}

//=================================================================
// Get tag = 0 or 1 for given list of edges
//=================================================================
void K_OCC::getEdgeTag(TopoDS_Shape& shape, PyObject* listEdges, std::vector<E_Int>& tag)
{
  if (listEdges == Py_None) return;
  TopTools_IndexedMapOfShape se;
  TopExp::MapShapes(shape, TopAbs_EDGE, se);
   
  E_Int nedgesTot = se.Extent();
  tag.resize(nedgesTot);
  for (E_Int i = 0; i < nedgesTot; i++) tag[i] = 0;
  for (E_Int no = 0; no < PyList_Size(listEdges); no++)
  {
    PyObject* noEdgeO = PyList_GetItem(listEdges, no);
    E_Int noEdge = PyInt_AsLong(noEdgeO);
    tag[noEdge-1] = 1;
  }
}

//=================================================================
// Get tag = 0 or 1 for given list of faces
//=================================================================
void K_OCC::getFaceTag(TopoDS_Shape& shape, PyObject* listFaces, std::vector<E_Int>& tag)
{
  if (listFaces == Py_None) return;
  TopTools_IndexedMapOfShape sf;
  TopExp::MapShapes(shape, TopAbs_FACE, sf);
   
  E_Int nfacesTot = sf.Extent();
  tag.resize(nfacesTot);
  for (E_Int i = 0; i < nfacesTot; i++) tag[i] = 0;
  for (E_Int no = 0; no < PyList_Size(listFaces); no++)
  {
    PyObject* noFaceO = PyList_GetItem(listFaces, no);
    E_Int noFace = PyInt_AsLong(noFaceO);
    tag[noFace-1] = 1;
  }
}

//=================================================================
// check a face - is it valid, is the face wire valid?
// return 0 (valid), 1 (invalid)
//=================================================================
E_Int K_OCC::checkFace(TopoDS_Face& F)
{
  if (F.IsNull()) 
  {
    printf("Error: checkFace: face is null.\n");    
    return 1;
  }
  BRepCheck_Analyzer analyzer(F);
  if (analyzer.IsValid()) 
  {
    printf("Info: checkFace: face is valid.\n");    
    return 0;
  }

  E_Int ret = 0;

  // Try to find error
  TopExp_Explorer exp;
  // 1. Check Face-level status
  for (exp.Init(F, TopAbs_FACE); exp.More(); exp.Next()) 
  {
    const TopoDS_Face& f = TopoDS::Face(exp.Current());
    const Handle(BRepCheck_Result)& res = analyzer.Result(f);
    if (!res.IsNull()) 
    {
      for (const auto& stat : res->Status()) 
      {
        if (stat != 0) 
        {
          printf("Error: checkFace: face error code=%d.\n", (int)stat);
          ret = 1;
        }
      }
    }
  }

  // 2. Check Wire-level status (commonly catches boundary loops issues)
  for (exp.Init(F, TopAbs_WIRE); exp.More(); exp.Next()) 
  {
    const TopoDS_Wire& w = TopoDS::Wire(exp.Current());
    const Handle(BRepCheck_Result)& res = analyzer.Result(w);
    if (!res.IsNull())
    {
      for (const auto& stat : res->Status()) 
      {
        if (stat != 0) 
        {
          printf("Error: checkFace: wire error code=%d.\n", (int)stat);
          ret = 1;
        }
      }
    }
  }
  return ret;
}
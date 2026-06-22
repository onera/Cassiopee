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

// Build a compound from list of faces
void buildfaceCompound(TopoDS_Shape& shape, PyObject* listFaces, TopoDS_Compound& compound)
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

// Build a compound from faces not in list
void buildOutfaceCompound(TopoDS_Shape& shape, PyObject* listFaces, TopoDS_Compound& compound)
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

// sew a shape
void sewShape(TopoDS_Shape& shape, E_Float tol)
{
  BRepBuilderAPI_Sewing sewingTool;
  sewingTool.Add(shape);
  sewingTool.Perform();
}

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
// untrim faces

#include "occ.h"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS.hxx"
#include "Geom_Surface.hxx"
#include "BRep_Tool.hxx"
#include "BRepTools.hxx"
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"

//=====================================================================
// untrim faces
//=====================================================================
PyObject* K_OCC::untrimFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces;
  if (!PYPARSETUPLE_(args, OO_, &hook, &listFaces)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPSURFACES;

  E_Int nfaces = PyList_Size(listFaces);

  // rebuild compound
  E_Int nfacesTot = surfaces.Extent();
  std::vector<E_Int> tag(nfacesTot);
  for (E_Int i = 0; i < nfacesTot; i++) tag[i] = 0;
  for (E_Int i = 0; i < nfaces; i++)
  {
    PyObject* faceNoO = PyList_GetItem(listFaces, i);
    E_Int faceNo = PyLong_AsLong(faceNoO);
    tag[faceNo-1] = 1;
  }
  double u1, u2, v1, v2;
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  for (E_Int i = 0; i < nfacesTot; i++)
  {
    const TopoDS_Face& F = TopoDS::Face(surfaces(i+1));
    if (tag[i] == 0) builder.Add(compound, F);
    else
    {
      // get the bound of the input face
      BRepTools::UVBounds(F, u1, u2, v1, v2);
      E_Float extension = 10.;
      u1 -= extension;
      u2 += extension;
      v1 -= extension;
      v2 += extension;
      TopLoc_Location faceLocation;
      Handle(Geom_Surface) surf = BRep_Tool::Surface(F, faceLocation);
      //surf->Bounds(u1, u2, v1, v2);
      printf("%g %g %g %g\n", u1,u2,v1,v2);
      TopoDS_Face face = BRepBuilderAPI_MakeFace(surf, u1, u2, v1, v2, Precision::Confusion());
      //TopoDS_Face face = BRepBuilderAPI_MakeFace(surf, true, Precision::Confusion()).Face();
      //BRepBuilderAPI_MakeFace maker;
      //maker.Init(surf, true, Precision::Confusion());
      //TopoDS_Face face = maker.Face();
      face.Location(faceLocation);
      face.Orientation(F.Orientation());
      builder.Add(compound, face);
      //builder.Add(compound, BRepTools::OuterWire(face));
    }
  }

  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after untrim: Nb edges=%d\n", se->Extent());
  printf("INFO: after untrim: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;

}

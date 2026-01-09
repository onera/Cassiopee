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
// rotate CAD

#include "occ.h"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include <gp_Ax1.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include "BRep_Builder.hxx"
#include <BRepBuilderAPI_Sewing.hxx>

//=====================================================================
// Rotate the full shape or some faces
// from an axis and an angle in degree
//=====================================================================
PyObject* K_OCC::rotate(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float angle; E_Float xc, yc, zc; E_Float xaxis, yaxis, zaxis;
  PyObject* listFaces; 
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ R_ O_, &hook, &xc, &yc, &zc, 
    &xaxis, &yaxis, &zaxis, &angle, &listFaces)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;

  gp_Pnt center(xc, yc, zc);
  gp_Dir direction(xaxis, yaxis, zaxis);
  gp_Ax1 axis(center, direction);

  angle = angle * M_PI / 180.;
  gp_Trsf myTrsf;
  myTrsf.SetRotation(axis, angle); 

  TopoDS_Shape* newshp = new TopoDS_Shape();

  if (listFaces == Py_None) // on all shape
  {
    BRepBuilderAPI_Transform myTransform(*shape, myTrsf);
    TopoDS_Shape tShape = myTransform.Shape();
    *newshp = tShape;
  }
  else
  {
    // Build a compound
    BRep_Builder builder;
    TopoDS_Compound shc;
    builder.MakeCompound(shc);
    E_Int nf = surfaces.Extent();
    std::vector<E_Int> nos(nf);
    for (E_Int i = 0; i < nf; i++) nos[i] = -1;
    
    E_Int nfaces = PyList_Size(listFaces);
    for (E_Int no = 0; no < nfaces; no++)
    {
      PyObject* noFaceO = PyList_GetItem(listFaces, no);
      E_Int noFace = PyInt_AsLong(noFaceO);
      nos[noFace-1] = no;
      const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
      builder.Add(shc, F);
    }

    BRepBuilderAPI_Transform myTransform(shc, myTrsf);
    TopoDS_Shape tShape = myTransform.Shape();

    // Rebuild
    TopTools_IndexedMapOfShape surfaces2;
    TopExp::MapShapes(tShape, TopAbs_FACE, surfaces2);  

    BRep_Builder builder2;
    TopoDS_Compound shc2;
    builder2.MakeCompound(shc2);
    for (E_Int i = 0; i < nf; i++)
    {
      if (nos[i] == -1)
      {
        const TopoDS_Face& F = TopoDS::Face(surfaces(i+1));
        builder2.Add(shc2, F);
      }
      else
      {
        const TopoDS_Face& F = TopoDS::Face(surfaces2(nos[i]+1));
        builder2.Add(shc2, F);
      }
    }
    BRepBuilderAPI_Sewing sewer(1.e-6);
    sewer.Add(shc2);
    sewer.Perform();
    *newshp = sewer.SewedShape();
  }

  delete shape;
  packet[0] = newshp;

  SETSHAPE(newshp);

  Py_INCREF(Py_None);
  return Py_None;
}

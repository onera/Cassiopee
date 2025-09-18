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
// translate CAD

#include "occ.h"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "BRepBuilderAPI_Transform.hxx"
#include "BRep_Builder.hxx"
#include <BRepBuilderAPI_Sewing.hxx>

//=====================================================================
// Translate the full shape or some faces
// from vector
//=====================================================================
PyObject* K_OCC::translate(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float dx, dy, dz; PyObject* listFaces; 
  if (!PYPARSETUPLE_(args, O_ TRRR_ O_, &hook, &dx, &dy, &dz, &listFaces)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  
  gp_Trsf myTrsf;
  myTrsf.SetTranslation(gp_Vec(dx, dy, dz)); // Translate by (dx, dy, dz)

  TopoDS_Shape* newshp = new TopoDS_Shape();

  if (listFaces == Py_None) // on all shape
  {
    BRepBuilderAPI_Transform myTransform(*shape, myTrsf);
    TopoDS_Shape tShape = myTransform.Shape();
    *newshp = tShape;
  }
  else // on face list
  {
    E_Int nfaces = PyList_Size(listFaces);
    // Build a compound
    BRep_Builder builder;
    TopoDS_Compound shc;
    builder.MakeCompound(shc);
    E_Int nf = surfaces.Extent();
    std::vector<E_Int> nos(nf);
    for (E_Int i = 0; i < nf; i++) nos[i] = -1;
    
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

  // Rebuild the hook
  delete shape;
  SETSHAPE(newshp);

  Py_INCREF(Py_None);
  return Py_None;
}

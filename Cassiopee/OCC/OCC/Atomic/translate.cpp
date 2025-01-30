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

//=====================================================================
// Translate the full shape or some faces
// from vector
//=====================================================================
PyObject* K_OCC::translate(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float dx, dy, dz; PyObject* listFaces; 
  if (!PYPARSETUPLE_(args, O_ TRRR_ O_, &hook, &dx, &dy, &dz, &listFaces)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  
  gp_Trsf myTrsf;
  myTrsf.SetTranslation(gp_Vec(dx, dy, dz)); // Translate by (dx, dy, dz)

  E_Int nfaces = PyList_Size(listFaces);

  TopoDS_Shape* newshp = new TopoDS_Shape();

  if (nfaces == 0) // on all shape
  {
    BRepBuilderAPI_Transform myTransform(*shp, myTrsf);
    TopoDS_Shape tShape = myTransform.Shape();
    *newshp = tShape;
  }
  else // on face list
  {
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
    *newshp = shc2;
  }

  // Rebuild the hook
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

  Py_INCREF(Py_None);
  return Py_None;
}

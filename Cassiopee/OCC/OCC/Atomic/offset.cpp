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
// offset CAD

#include "occ.h"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "BRep_Builder.hxx"
#include "BRepOffsetAPI_MakeOffsetShape.hxx"

//=====================================================================
// offset the full shape or some faces
// from a constant distance
// if no face list: full shape offset (with intersection)
// if face list : simple by face offset
//=====================================================================
PyObject* K_OCC::offset(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float distance;
  PyObject* listFaces; 
  if (!PYPARSETUPLE_(args, O_ R_ O_, &hook, &distance, &listFaces)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  
  TopoDS_Shape* newshp = new TopoDS_Shape();

  if (listFaces == Py_None) // on all shape
  {    
    BRepOffsetAPI_MakeOffsetShape offsetMaker;
    //offsetMaker.PerformBySimple(*shape, distance);
    offsetMaker.PerformByJoin(
      *shape,
      distance,          // offset
      1e-3,             // tolerance
      BRepOffset_Skin,
      Standard_True,  // compute intersections
      Standard_False, // detect self-intersections
      GeomAbs_Arc     // join type
    );
    *newshp = offsetMaker.Shape();
  }
  else
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
    
    BRepOffsetAPI_MakeOffsetShape offsetMaker; 
    offsetMaker.PerformBySimple(shc, distance);
    TopoDS_Shape tShape = offsetMaker.Shape();

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
  delete shape;
  SETSHAPE(newshp);

  Py_INCREF(Py_None);
  return Py_None;
}

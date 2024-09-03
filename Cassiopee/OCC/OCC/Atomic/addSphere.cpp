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

#include "occ.h"
#include "TopoDS_Shape.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"

//=====================================================================
// Add a sphere to hook
//=====================================================================
PyObject* K_OCC::addSphere(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float xc, yc, zc, R;
  if (!PYPARSETUPLE_(args, O_ RRRR_, &hook, &xc, &yc, &zc, &R)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  /* previous shape */
  //TopoDS_Shape* psh = (TopoDS_Shape*)packet[0];
  //TopAbs_ShapeEnum ptype = curShape.ShapeType();
  // == TopAbs_COMPOUND 
  //for(TopoDS_Iterator anExp(curShape); anExp.More(); anExp.Next()){
  //const TopoDS_Shape &curShape1 = anExp.Value();

  /* new sphere */
  gp_Pnt center(xc, yc, zc);
  BRepPrimAPI_MakeSphere makerSphere(center, R);
  TopoDS_Shape sh = makerSphere.Shape();

  // Building the Resulting Compound
  /*
  TopoDS_Compound aRes;
  BRep_Builder aBuilder;
  aBuilder.MakeCompound(aRes);
  aBuilder.Add(aRes, myBody);
  aBuilder.Add(aRes, myThreading);
  */

  /* export */
  TopoDS_Shape* newshp = new TopoDS_Shape(sh);
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

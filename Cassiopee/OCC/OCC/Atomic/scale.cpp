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
// scale CAD

#include "occ.h"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "BRepBuilderAPI_Transform.hxx"

//=====================================================================
// scale the full shape or some faces
// from a constant scale factor
//=====================================================================
PyObject* K_OCC::scale(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float factor; E_Float x0, y0, z0;
  if (!PYPARSETUPLE_(args, O_ R_ TRRR_, &hook, &factor, 
    &x0, &y0, &z0)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];

  // idem scale
  gp_Trsf myTrsf;
  myTrsf.SetScale(gp_Pnt(x0, y0, z0), factor);
  
  BRepBuilderAPI_Transform myTransform(*shp, myTrsf);
  TopoDS_Shape tShape = myTransform.Shape();

  TopoDS_Shape* newshp = new TopoDS_Shape();
  *newshp = tShape;

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

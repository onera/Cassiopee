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

//=====================================================================
// Rotate the full shape or some faces
// from an axis and an angle (point?
//=====================================================================
PyObject* K_OCC::rotate(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float angle; E_Float xc, yc, zc; E_Float xaxis, yaxis, zaxis;
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ R_, &hook, &xc, &yc, &zc, &xaxis, &yaxis, &zaxis, &angle)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];

  gp_Pnt center(xc, yc, zc);
  gp_Dir direction(xaxis, yaxis, zaxis);
  gp_Ax1 axis(center, direction);

  //printf("angle=%g\n", angle);
  angle = angle * M_PI / 180.;
  gp_Trsf myTrsf;
  myTrsf.SetRotation(axis, angle); 

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
  //printf("INFO: after fix: Nb edges=%d\n", se->Extent());
  //printf("INFO: after fix: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

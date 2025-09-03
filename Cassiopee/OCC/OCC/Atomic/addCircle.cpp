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

#include "occ.h"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include <gp_Ax2.hxx>
#include <gp_Dir.hxx>
#include <GC_MakeCircle.hxx>
#include <Geom_Circle.hxx>

//=====================================================================
// Add a circle to CAD hook
//=====================================================================
PyObject* K_OCC::addCircle(PyObject* self, PyObject* args)
{
  PyObject* hook; 
  E_Float xc, yc, zc, ax, ay, az, R;
  E_Int makeFace;
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ R_ I_, &hook, &xc, &yc, &zc, 
    &ax, &ay, &az, &R, &makeFace)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  //TopoDS_Shape* shp = (TopoDS_Shape*) packet[0];
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];

  /* new circle */
  gp_Pnt pc(xc, yc, zc); // Center
  gp_Dir normal(ax, ay, az); // Normal vector
  gp_Ax2 axis(pc, normal);

  GC_MakeCircle circleMaker(axis, R);
  Handle(Geom_Circle) circle = circleMaker.Value();

  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circle);
  TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge);
  TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);

  // Rebuild a single compound
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
    
  for (E_Int i = 1; i <= surfaces.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    builder.Add(compound, F);
  }
  if (makeFace == 1) builder.Add(compound, face);
  else builder.Add(compound, wire);
  
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
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
  printf("INFO: after addCircle: Nb edges=%d\n", se->Extent());
  printf("INFO: after addCircle: Nb faces=%d\n", sf->Extent());
  
  Py_INCREF(Py_None);
  return Py_None;

}

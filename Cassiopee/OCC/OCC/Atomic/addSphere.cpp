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

#include "occ.h"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"
#include "BRep_Builder.hxx"

//=====================================================================
// Add a sphere to CAD hook
//=====================================================================
PyObject* K_OCC::addSphere(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Float xc, yc, zc, R;
  if (!PYPARSETUPLE_(args, O_ TRRR_ R_, &hook, &xc, &yc, &zc, &R)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  GETMAPEDGES;

  /* new sphere */
  gp_Pnt center(xc, yc, zc);
  BRepPrimAPI_MakeSphere makerSphere(center, R);
  TopoDS_Shape sphere = makerSphere.Shape();

  // Rebuild a single compound
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
    
  for (E_Int i = 1; i <= surfaces.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    builder.Add(compound, F);
  }
  for (E_Int i = 1; i <= edges.Extent(); i++)
  {
    TopoDS_Edge E = TopoDS::Edge(edges(i));
    builder.Add(compound, E);
  }
  // Add the sphere faces
  TopTools_IndexedMapOfShape* sfs = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(sphere, TopAbs_FACE, *sfs);
  TopTools_IndexedMapOfShape& surfaces2 = *(TopTools_IndexedMapOfShape*)sfs;
  for (E_Int i = 1; i <= surfaces2.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces2(i));
    builder.Add(compound, F);
  }
  delete sfs;

  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after addSphere: Nb edges=%d\n", se->Extent());
  printf("INFO: after addSphere: Nb faces=%d\n", sf->Extent());
  
  Py_INCREF(Py_None);
  return Py_None;
}

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
// revolve a list of edges

#include "occ.h"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax1.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>

//=====================================================================
// Revolve
//=====================================================================
PyObject* K_OCC::revolve(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges;
  E_Float cx, cy, cz;
  E_Float ax, ay, az; 
  E_Float angle; 
  if (!PYPARSETUPLE_(args, OO_ TRRR_ TRRR_ R_, &hook, &listEdges, 
    &cx, &cy, &cz, &ax, &ay, &az, &angle)) return NULL;

  GETSHAPE;

  //TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];

  gp_Pnt center(cx, cy, cz);
  gp_Ax1 axis(center, gp_Dir(ax, ay, az));
  
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);

  E_Int nedges = PyList_Size(listEdges);
  for (E_Int i = 0; i < nedges; i++)
  {
    PyObject* noO = PyList_GetItem(listEdges, i);
    E_Int no = PyInt_AsLong(noO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    builder.Add(compound, E);
  }

  TopoDS_Shape revolvedSurface = BRepPrimAPI_MakeRevol(compound, axis);

  // rebuild
  TopoDS_Shape* newshp = new TopoDS_Shape(revolvedSurface);
  
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after revolve: Nb edges=%d\n", se->Extent());
  printf("INFO: after revolve: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

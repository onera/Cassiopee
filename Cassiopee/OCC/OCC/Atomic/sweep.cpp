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
// sweep a list of edges by a list of edges

#include "occ.h"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "Precision.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"
#include "GeomFill_BSplineCurves.hxx"
#include "Geom_BSplineCurve.hxx"
#include "Geom_BSplineSurface.hxx"
#include "Geom_Surface.hxx"
#include "Geom_Curve.hxx"
#include "BRepOffsetAPI_MakePipe.hxx"

//=====================================================================
// Sweep
//=====================================================================
PyObject* K_OCC::sweep(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listProfiles; PyObject* listPaths; 
  if (!PYPARSETUPLE_(args, OOO_ , &hook, &listProfiles, &listPaths)) return NULL;

  GETSHAPE;
  GETMAPEDGES;

  E_Int nprofiles = PyList_Size(listProfiles);
  E_Int npaths = PyList_Size(listPaths);
  TopoDS_Shape* newshp = NULL;

  // Try merge listProfiles in a single wire
  TopoDS_Wire W1;
  for (E_Int i = 0; i < nprofiles; i++)
  {
    PyObject* noO = PyList_GetItem(listProfiles, i);
    E_Int no = PyInt_AsLong(noO);

    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    if (i == 0)
    {
      W1 = BRepBuilderAPI_MakeWire(E);
    }
    else
    {
      W1 = BRepBuilderAPI_MakeWire(W1, E);
    }
  }

  // Try merge listPaths in a single wire
  TopoDS_Wire W2;
  for (E_Int i = 0; i < npaths; i++)
  {
    PyObject* noO = PyList_GetItem(listPaths, i);
    E_Int no = PyInt_AsLong(noO);

    const TopoDS_Edge& E = TopoDS::Edge(edges(no));
    if (i == 0)
    {
      W2 = BRepBuilderAPI_MakeWire(E);
    }
    else
    {
      W2 = BRepBuilderAPI_MakeWire(W1, E);
    }
  }

  // Sweep
  BRepOffsetAPI_MakePipe pipe(W2, W1);
  newshp = new TopoDS_Shape(pipe.Shape());

  // Rebuild the hook
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after sweep: Nb edges=%d\n", se->Extent());
  printf("INFO: after sweep: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

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
#include "BRep_Tool.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx" 
#include "GCPnts_UniformDeflection.hxx"
#include "GCPnts_UniformAbscissa.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp_Explorer.hxx"

// ============================================================================
/* Identify loops in edges  */
// ============================================================================
PyObject* K_OCC::identifyLoopsInEdges(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;  

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  
  TopExp_Explorer expl;
  PyObject* pl = PyList_New(0);
  for (E_Int i=1; i <= surfaces.Extent(); i++)
  {
    PyObject* l = PyList_New(0);
    //const TopoDS_Face& F = TopoDS::Face(surfaces(i));
    for (expl.Init(surfaces(i), TopAbs_EDGE); expl.More(); expl.Next())
    {
      const TopoDS_Edge& E = TopoDS::Edge(expl.Current());
      E_Int id = edges.FindIndex(E); // index de l'edge dans edges
      PyObject* o = Py_BuildValue("i", id);
      PyList_Append(l, o); Py_DECREF(o);
    }
    PyList_Append(pl, l); Py_DECREF(l);
  }  
    
  return pl;
}
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
// CAD analyser
#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "TopoDS_Wire.hxx"
#include "BRepCheck_Wire.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "ShapeAnalysis.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx" 
#include "GCPnts_UniformDeflection.hxx"
#include "GCPnts_UniformAbscissa.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp_Explorer.hxx"
#include "Geom2d_Curve.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepGProp.hxx"
#include "GProp_GProps.hxx"

E_Float __getLength(const TopoDS_Edge& E);

//=====================================================================
// Return min / max / mean length of all edges
//=====================================================================
PyObject* K_OCC::analyseEdges(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;
    
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  E_Float emin = K_CONST::E_MAX_FLOAT;
  E_Float emax = -K_CONST::E_MAX_FLOAT;
  E_Float ltot = 0.;

  for (E_Int i=1; i <= edges.Extent(); i++)
  {
    const TopoDS_Edge& E = TopoDS::Edge(edges(i));
    E_Float l = __getLength(E);
    emax = std::max(emax, l);
    emin = std::min(emin, l);
    ltot += l;
  }
  E_Float emean = ltot / edges.Extent();

  // calcul hmax, hausd : 20 pts par edges moyen
  E_Float hmax = emean / 20.;
  // calcul du nbre de points sur la longueur totale
  E_Int Np = ltot / hmax;
  if (Np > 20000) hmax = ltot / 20000.;
  E_Float hmin = emin / 20.;
  E_Float hausd = hmax / 10.;
  printf("INFO: suggested hmin=%g hmax=%g hausd=%g\n", hmin, hmax, hausd);
  return Py_BuildValue("ddd", hmin, hmax, hausd);
}

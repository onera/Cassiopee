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

#include "TColStd_HSequenceOfTransient.hxx"
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Tool.hxx"
#include "GCPnts_AbscissaPoint.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Wire.hxx"
#include "GeomAdaptor_Curve.hxx"

// sum edge length
void edgeLength(TopTools_IndexedMapOfShape& edges, E_Int no, E_Float& L)
{
  E_Float length = 0.;
  Standard_Real first, last;
  const TopoDS_Edge& E = TopoDS::Edge(edges(no));
  Handle(Geom_Curve) curve = BRep_Tool::Curve(E, first, last);
  if (curve.IsNull()) return;
  GeomAdaptor_Curve adaptorCurve(curve, first, last);
  length = GCPnts_AbscissaPoint::Length(adaptorCurve, first, last);
  L += length;
}

// ============================================================================
/* Return the min/max of edge lengths of given face */
// ============================================================================
PyObject* K_OCC::getEdgeLength(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges;
  if (!PYPARSETUPLE_(args, O_ O_, &hook, &listEdges)) return NULL;

  GETPACKET;
  GETMAPEDGES;

  E_Float L = 0.0;
  if (listEdges == Py_None) // all edges of topshape
  {
    for (E_Int noEdge=1; noEdge <= edges.Extent(); noEdge++)
    {
      edgeLength(edges, noEdge, L);
    }
  }
  else
  {
    E_Int nedges = PyList_Size(listEdges);
    for (E_Int no = 0; no < nedges; no++)
    {
      PyObject* noEdgeO = PyList_GetItem(listEdges, no);
      E_Int noEdge = PyInt_AsLong(noEdgeO);
      edgeLength(edges, noEdge, L);
    }
  }
  return Py_BuildValue("d", L);
} 
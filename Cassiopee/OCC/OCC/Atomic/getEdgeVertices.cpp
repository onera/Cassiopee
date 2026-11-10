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
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Tool.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "TopoDS_Vertex.hxx"

// ============================================================================
/* Return the coordinates of extreme vertices of edge */
// ============================================================================
PyObject* K_OCC::getEdgeVertices(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges;
  if (!PYPARSETUPLE_(args, O_ O_, &hook, &listEdges)) return NULL;

  GETPACKET;
  GETMAPEDGES;

  PyObject* l = PyList_New(0);
  
  E_Int nedges = PyList_Size(listEdges);
  for (E_Int no = 0; no < nedges; no++)
  {
    PyObject* noEdgeO = PyList_GetItem(listEdges, no);
    E_Int noEdge = PyInt_AsLong(noEdgeO);
    const TopoDS_Edge& E = TopoDS::Edge(edges(noEdge));
    
    TopoDS_Vertex vFirst, vLast;
    
    // Extract the start and end vertices, respecting edge orientation
    TopExp::Vertices(E, vFirst, vLast);
    
    // Get the 3D coordinates (gp_Pnt) from the vertices
    gp_Pnt pFirst = BRep_Tool::Pnt(vFirst);
    gp_Pnt pLast = BRep_Tool::Pnt(vLast);
    BRepAdaptor_Curve adaptor(E);
    // Get the parameter bounds of the edge
    Standard_Real firstParam = adaptor.FirstParameter();
    Standard_Real lastParam  = adaptor.LastParameter();
    Standard_Real midParam = (firstParam + lastParam) * 0.5;
    gp_Pnt pMid = adaptor.Value(midParam);
    
    // Example usage:
    // pFirst.X(), pFirst.Y(), pFirst.Z();
    // pLast.X(), pLast.Y(), pLast.Z();
    PyObject* tuple = Py_BuildValue("[(ddd)(ddd)(ddd)]",
      pFirst.X(), pFirst.Y(), pFirst.Z(),
      pLast.X(), pLast.Y(), pLast.Z(),
      pMid.X(), pMid.Y(), pMid.Z());

    PyList_Append(l, tuple); Py_DECREF(tuple);
  }

  return l;
}

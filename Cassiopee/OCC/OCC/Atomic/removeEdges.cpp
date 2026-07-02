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
// Remove CAD edges and rebuild compound
#include "occ.h"

#include "TopoDS.hxx"
#include "TopoDS_Edge.hxx"
#include "BRep_Tool.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "BRep_Builder.hxx"

//=====================================================================
// Remove some faces and rebuild compound
// output edgeMap et faceMap 
//=====================================================================
PyObject* K_OCC::removeEdges(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges; 
  if (!PYPARSETUPLE_(args, OO_, &hook, &listEdges)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPEDGES;

  // get top shape
  E_Int nbEdges = edges.Extent();

  ShapeBuild_ReShape reshaper;
  for (E_Int no = 0; no < PyList_Size(listEdges); no++)
  {
    PyObject* noEdgeO = PyList_GetItem(listEdges, no);
    E_Int noEdge = PyInt_AsLong(noEdgeO);
    if (noEdge >= 1 && noEdge <= nbEdges)
    {
      const TopoDS_Edge& E = TopoDS::Edge(edges(noEdge));
      printf("Info: removing Edge " SF_D_ "\n", noEdge);
      reshaper.Remove(E);
    }
    else printf("Warning: removeEdges: invalid edge number.\n");
  }
  TopoDS_Shape shc = reshaper.Apply(*shape);

  // export
  TopoDS_Shape* newshp = new TopoDS_Shape(shc);

  SETSHAPE(newshp);

  printf("INFO: after removeEdges: Nb edges=%d\n", se->Extent());
  printf("INFO: after removeEdges: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

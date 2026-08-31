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

#include <ShapeUpgrade_UnifySameDomain.hxx>
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "BRep_Builder.hxx"
#include "TopoDS.hxx"
#include "TopExp.hxx"
#include "TopoDS_Shape.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepTools_ReShape.hxx"

#include "BRepTools_WireExplorer.hxx"
#include "GeomConvert_CompCurveToBSplineCurve.hxx"
#include "Geom_TrimmedCurve.hxx"
#include "Geom_BSplineCurve.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"

// fuse a wire in a single edge
TopoDS_Edge wireToSingleEdge(const TopoDS_Wire& wire)
{
  GeomConvert_CompCurveToBSplineCurve curveConverter;
    
  // Iterate through edges in topological order along the wire
  BRepTools_WireExplorer explorer(wire);
  for (; explorer.More(); explorer.Next())
  {
    const TopoDS_Edge& currentEdge = explorer.Current();
    Standard_Real first, last;
    Handle(Geom_Curve) curve = BRep_Tool::Curve(currentEdge, first, last);
        
    if (!curve.IsNull())
    {
      // Truncate curve to edge bounds if necessary
      Handle(Geom_TrimmedCurve) trimmedCurve = new Geom_TrimmedCurve(curve, first, last);
            
      // Add curve segment (false = do not reverse parameterization unless needed)
      curveConverter.Add(trimmedCurve, Precision::Confusion());
    }
  }
    
  // Get the combined B-spline curve
  Handle(Geom_BSplineCurve) mergedCurve = curveConverter.BSplineCurve();
    
  if (mergedCurve.IsNull()) 
  {
    return TopoDS_Edge(); // Return null edge on failure
  }
    
  // Build and return the new single edge
  return BRepBuilderAPI_MakeEdge(mergedCurve);
}

// ============================================================================
/* Merge a list of edges in a single edge */
// ============================================================================
PyObject* K_OCC::mergeEdges(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listEdges;
  if (!PYPARSETUPLE_(args, OO_, &hook, &listEdges)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPEDGES;
  GETMAPSURFACES;

  TopoDS_Shape* newshp = NULL;
  if (listEdges == Py_None)
  {
    // Try merging all edges
    const TopoDS_Edge& E = TopoDS::Edge(edges(1));
    TopoDS_Wire W = BRepBuilderAPI_MakeWire(E);
    for (E_Int i = 2; i <= edges.Extent(); i++)
    {
      const TopoDS_Edge& E = TopoDS::Edge(edges(i));
      W = BRepBuilderAPI_MakeWire(W, E);
    }
    // rebuild new shape
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    builder.Add(compound, W);
    for (E_Int i = 1; i <= surfaces.Extent(); i++)
    {
      TopoDS_Face F = TopoDS::Face(surfaces(i));
      builder.Add(compound, F);
    }
    newshp = new TopoDS_Shape(compound);
  }
  else
  {
    E_Int nedges = PyList_Size(listEdges);

    PyObject* noEdgeO = PyList_GetItem(listEdges, 0);
    E_Int noEdge = PyInt_AsLong(noEdgeO);
    printf("no edge=%d\n", noEdge); fflush(stdout);
    TopoDS_Edge E = TopoDS::Edge(edges(noEdge));
    TopoDS_Wire W = BRepBuilderAPI_MakeWire(E);
    
    for (E_Int no = 1; no < nedges; no++)
    {
      noEdgeO = PyList_GetItem(listEdges, no);
      noEdge = PyInt_AsLong(noEdgeO);
      printf("no edge=%d\n", noEdge); fflush(stdout);
      E = TopoDS::Edge(edges(noEdge));
      W = BRepBuilderAPI_MakeWire(W, E);
    }

    printf("Getting wire\n"); fflush(stdout);

    /*
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);

    E_Int* tag = new E_Int [surfaces.Extent()];
    for (E_Int i = 0; i < edges.Extent(); i++) tag[i] = 1;
    for (E_Int no = 0; no < nedges; no++)
    {
      noEdgeO = PyList_GetItem(listEdges, no);
      noEdge = PyInt_AsLong(noEdgeO);
      tag[noEdge-1] = 0;
    }
    
    for (E_Int i = 1; i <= edges.Extent(); i++)
    {
      if (tag[i-1] == 1)
      {
        TopoDS_Edge E = TopoDS::Edge(edges(i));
        builder.Add(compound, E);
      }
    }
    for (E_Int i = 1; i <= surfaces.Extent(); i++)
    {
      TopoDS_Face F = TopoDS::Face(surfaces(i));
      builder.Add(compound, F);
    }
    delete [] tag;
    */

    BRepTools_ReShape reshaper;
    noEdgeO = PyList_GetItem(listEdges, 0);
    noEdge = PyInt_AsLong(noEdgeO);
    E = TopoDS::Edge(edges(noEdge));
    TopoDS_Edge WE = wireToSingleEdge(W);
    reshaper.Replace(E, WE);
    printf("replaced edge 0 with merged wire\n"); fflush(stdout);

    for (E_Int no = 1; no < nedges; no++)
    {
      noEdgeO = PyList_GetItem(listEdges, no);
      noEdge = PyInt_AsLong(noEdgeO);
      E = TopoDS::Edge(edges(noEdge));
      reshaper.Remove(E);
    }
    printf("removed edge\n"); fflush(stdout);

    /*
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    for (E_Int i = 1; i <= surfaces.Extent(); i++)
    {
      TopoDS_Face F = TopoDS::Face(surfaces(i));
      //reshaper.Apply(F);
      builder.Add(compound, F);
    }*/

    TopoDS_Shape compound = reshaper.Apply(*shape);
    newshp = new TopoDS_Shape(compound);
  }

  delete shape;
  SETSHAPE(newshp);
  
  printf("INFO: after mergeEdges: Nb edges=%d\n", se->Extent());
  printf("INFO: after mergeEdges: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
  
}
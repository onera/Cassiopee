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
// Query edge
#include "occ.h"
#include "TopoDS_Edge.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GeomAbs_CurveType.hxx"
#include "Geom_BSplineCurve.hxx"
#include "Geom_BezierCurve.hxx"
#include "gp_Lin.hxx"
#include "gp_Circ.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "gp_Pnt.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopoDS.hxx"

// build data dict from E
PyObject* queryEdge1(const TopoDS_Edge& E) 
{
  // return dict
  PyObject* dict=NULL;

  // 1. Wrap the edge into an adaptor
  BRepAdaptor_Curve adaptor(E);
  GeomAbs_CurveType curveType = adaptor.GetType();

  // 2. Determine the type and extract data accordingly
  switch (curveType) 
  {
    case GeomAbs_Line:
    {
      printf("Info: queryEdge: edge is a line.\n");
      //gp_Lin line = adaptor.Line();
      Standard_Real first = adaptor.FirstParameter();
      Standard_Real last = adaptor.LastParameter();
      gp_Pnt pStart = adaptor.Value(first);
      gp_Pnt pEnd = adaptor.Value(last);
      dict = Py_BuildValue("{s:s, s:O, s:O}", 
                           "type", "Line", 
                           "P0", Py_BuildValue("(ddd)", pStart.X(), pStart.Y(), pStart.Z()),
                           "P1", Py_BuildValue("(ddd)", pEnd.X(), pEnd.Y(), pEnd.Z()));
      break;
    }
    case GeomAbs_Circle: 
    {
      printf("Info: queryEdge: edge is a circle.\n");
      gp_Circ circle = adaptor.Circle();
      gp_Pnt center = circle.Location();
      gp_Ax1 axis = circle.Axis();
      gp_Dir normalDir = axis.Direction();
      Standard_Real radius = circle.Radius();

      // check if circle is full
      Standard_Real firstParam = adaptor.FirstParameter();
      Standard_Real lastParam = adaptor.LastParameter();
      Standard_Real sweepAngle = lastParam - firstParam;

      bool isFullCircle = (std::abs(sweepAngle - (2.0 * M_PI)) < 1e-7);
      if (isFullCircle)
      {
        dict = Py_BuildValue("{s:s, s:O, s:O, s:d}", 
                            "type", "Circle", 
                            "center", Py_BuildValue("(ddd)", center.X(), center.Y(), center.Z()),
                            "axis", Py_BuildValue("(ddd)", normalDir.X(), normalDir.Y(), normalDir.Z()),
                            "radius", radius);
      }
      else
      {
        //gp_Pnt pStart = adaptor.Value(firstParam);
        //gp_Pnt pEnd = adaptor.Value(lastParam);
        dict = Py_BuildValue("{s:s, s:O, s:O, s:d}", 
                             "type", "Arc", 
                             "center", Py_BuildValue("(ddd)", center.X(), center.Y(), center.Z()),
                             "axis", Py_BuildValue("(ddd)", normalDir.X(), normalDir.Y(), normalDir.Z()),
                             "radius", radius);
      }
      break;
    }
    case GeomAbs_BSplineCurve: 
    {
      printf("Info: queryEdge: edge is a spline.\n");
      Handle(Geom_BSplineCurve) bspline = adaptor.BSpline();
      Standard_Integer nbPoles = bspline->NbPoles();
            
      TColgp_Array1OfPnt poles(1, nbPoles);
      bspline->Poles(poles); // Extracts the control points (poles)

      PyObject* array = K_NUMPY::buildNumpyArray(nbPoles, 3, 0);
      E_Float* ap; E_Int size;
      K_NUMPY::getFromNumpyArray(array, ap, size);
  
      for (Standard_Integer i = 1; i <= nbPoles; ++i) 
      {
        // Access each control point coordinate
        ap[3*(i-1)] = poles(i).X();
        ap[3*(i-1)+1] = poles(i).Y();
        ap[3*(i-1)+2] = poles(i).Z();
      }
      E_Int degree = bspline->Degree();
      dict = Py_BuildValue("{s:s, s:O, s:i}", 
                           "type", "Spline1",
                           "Points", array,
                           "degree", degree);
      Py_DECREF(array);
      break;
    }
    case GeomAbs_BezierCurve: 
    {
      printf("Info: queryEdge: edge is a bezier.\n");
      Handle(Geom_BezierCurve) bezier = adaptor.Bezier();
      Standard_Integer nbPoles = bezier->NbPoles();
            
      TColgp_Array1OfPnt poles(1, nbPoles);
      bezier->Poles(poles); // Extracts Bezier control points

      PyObject* array = K_NUMPY::buildNumpyArray(nbPoles, 3, 0);
      E_Float* ap; E_Int size;
      K_NUMPY::getFromNumpyArray(array, ap, size);
  
      for (Standard_Integer i = 1; i <= nbPoles; ++i) 
      {
        // Access each control point coordinate
        ap[3*(i-1)] = poles(i).X();
        ap[3*(i-1)+1] = poles(i).Y();
        ap[3*(i-1)+2] = poles(i).Z();
      }
      E_Int degree = bezier->Degree();
      dict = Py_BuildValue("{s:s, s:O, s:i}", 
                           "type", "Bezier",
                           "Points", array,
                           "degree", degree);
      Py_DECREF(array);
      break;
    }
    case GeomAbs_Ellipse: 
    {
      printf("Info: queryEdge: edge is an ellipse.\n");
      gp_Elips ellipse = adaptor.Ellipse();
      gp_Pnt center = ellipse.Location();
      gp_Ax1 axis = ellipse.Axis();
      gp_Dir normalDir = axis.Direction();
      Standard_Real majorR = ellipse.MajorRadius();
      Standard_Real minorR = ellipse.MinorRadius();
      dict = Py_BuildValue("{s:s, s:O, s:O, s:d, s:d}", 
                           "type", "Ellipse",
                           "center", Py_BuildValue("(ddd)", center.X(), center.Y(), center.Z()),
                           "axis", Py_BuildValue("(ddd)", normalDir.X(), normalDir.Y(), normalDir.Z()),
                           "R1", majorR,
                           "R2", minorR);
      break;
    }
    default:
    {
      printf("Info: queryEdge: is not implemented for %d type of curve.\n", curveType);
      break;
    }
  }
  return dict;
}

//=====================================================================
// Query data for an edge
// return a dict
//=====================================================================
PyObject* K_OCC::queryEdge(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int edgeNo;
  if (!PYPARSETUPLE_(args, O_ I_, &hook, &edgeNo)) return NULL;

  GETPACKET;
  GETMAPEDGES;

  const TopoDS_Edge& E = TopoDS::Edge(edges(edgeNo));
  PyObject* dict = queryEdge1(E);

  return dict;
}

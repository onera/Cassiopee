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
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "gp_Ax2.hxx"
#include "gp_Dir.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "Geom2d_Curve.hxx"
#include "Geom_BSplineCurve.hxx"
#include "GeomAPI_PointsToBSpline.hxx"
#include "GeomAPI_Interpolate.hxx"

//=====================================================================
// Add a super ellipse to CAD hook
//=====================================================================
PyObject* K_OCC::addSuperEllipse(PyObject* self, PyObject* args)
{
  PyObject* hook; 
  E_Float xc, yc, zc, a, b;
  E_Int n, samples;
  E_Int makeFace;
  if (!PYPARSETUPLE_(args, O_ TRRR_ RR_ II_ I_, &hook, 
        &xc, &yc, &zc, &a, &b, 
        &n, &samples, &makeFace)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  GETMAPEDGES;

  Handle(TColgp_HArray1OfPnt) pointsArray = new TColgp_HArray1OfPnt(1, samples);
  if (n == 0) n = 1;
  E_Float exp = 2.0 / E_Float(n);
  E_Float t, x, y, z, cx, sx;

  for (E_Int i = 1; i <= samples; i++) 
  {
    t = 2.0 * M_PI * E_Float(i - 1) / E_Float(samples);
    cx = cos(t); sx = sin(t);
    x = a * K_FUNC::E_sign(cx) * pow(K_FUNC::E_abs(cx), exp);
    x += xc;
    y = b * K_FUNC::E_sign(sx) * pow(K_FUNC::E_abs(sx), exp);
    y += yc;
    z = zc;
    pointsArray->SetValue(i, gp_Pnt(x, y, z));
  }

  // Interpolate a closed BSpline through sampled points
  const Standard_Real tol = 1.0e-7;
  GeomAPI_Interpolate interp(pointsArray, /*Periodic*/ Standard_True, /*Tolerance*/ tol);
  interp.Perform();
  Handle(Geom_BSplineCurve) spline = interp.Curve();
  
  //Handle(Geom_BSplineCurve) spline = GeomAPI_PointsToBSpline(
  //      pointsArray, 
  //      Geom_BSplineCurve::MaxDegree() - 6,
  //      Geom_BSplineCurve::MaxDegree(),
  //      GeomAbs_C2, 
  //      Precision::Confusion()).Curve();

  //spline->SetNotPeriodic();
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(spline);

  TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge);
  TopoDS_Face face = BRepBuilderAPI_MakeFace(wire);

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
  if (makeFace == 1) builder.Add(compound, face);
  else builder.Add(compound, wire);
  
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after addSuperEllipse: Nb edges=%d\n", se->Extent());
  printf("INFO: after addSuperEllipse: Nb faces=%d\n", sf->Extent());
  
  Py_INCREF(Py_None);
  return Py_None;
}

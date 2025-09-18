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

#include "occ.h"
#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "BRepPrimAPI_MakeSphere.hxx"
#include "BRep_Builder.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeWire.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "Geom_BSplineCurve.hxx"
#include "GeomAPI_PointsToBSpline.hxx"
#include "TColStd_HArray1OfBoolean.hxx"
#include "TColgp_Array1OfVec.hxx"
#include "GeomAPI_Interpolate.hxx"

//=====================================================================
// Add a spline to CAD hook
//=====================================================================
PyObject* K_OCC::addSpline(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* opc; E_Int method = 0;
  if (!PYPARSETUPLE_(args, OO_ I_, &hook, &opc, &method)) return NULL;

  GETSHAPE;
  GETMAPSURFACES;
  GETMAPEDGES;

  /* get control points (method0) or through points (method1) */
  K_FLD::FldArrayF* pc;
  E_Int ret = K_NUMPY::getFromNumpyArray(opc, pc);
  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addSpline: invalid array.");
    return NULL;
  }
  E_Float* x = pc->begin(1);
  E_Float* y = pc->begin(2);
  E_Float* z = pc->begin(3);
  
  // incoming numpy
  E_Int ncp = pc->getSize();
  //for (E_Int i = 0; i < ncp; i++) printf("%d : %g %g %g\n", i, x[i], y[i], z[i]);

  TopoDS_Edge edge;

  if (method == 0) // linear knots
  {
    // compute total length and copy control points
    std::vector<gp_Pnt>::const_iterator iter;
    gp_Pnt lastP(x[0], y[0], z[0]);
    E_Float totalLen = 0.;
    TColgp_Array1OfPnt cp(1, ncp);
    for (E_Int i = 1; i <= ncp; i++) 
    {
      gp_Pnt p(x[i-1],y[i-1],z[i-1]);
      E_Float segLen = p.Distance(lastP);
      totalLen += segLen;
      cp.SetValue(i, p);
      lastP = p;
    }

    // compute knots
    TColStd_Array1OfReal knots(1, ncp);
    TColStd_Array1OfInteger mults(1, ncp);
  
    E_Float lastKnot = 0.;
    lastP.SetCoord(x[0], y[0], z[0]);
    for (E_Int i = 1; i <= ncp; ++i) 
    {
      if (i == 1 || i == ncp) mults.SetValue(i, 2);
      else mults.SetValue(i, 1);

      E_Float knot = cp.Value(i).Distance(lastP)/totalLen + lastKnot;
      knots.SetValue(i, knot);

      lastKnot = knot;
      lastP = cp.Value(i);
    }
    Handle(Geom_BSplineCurve) spline = new Geom_BSplineCurve(cp, knots, mults, 1, false);
    edge = BRepBuilderAPI_MakeEdge(spline);
  }
  else if (method == 1) // approximation through points
  {
    TColgp_Array1OfPnt pointsArray(1, static_cast<Standard_Integer>(ncp));
    for (E_Int i = 1; i <= ncp; i++) 
    {
      gp_Pnt p(x[i-1],y[i-1],z[i-1]);
      pointsArray.SetValue(i, p);
    }

    Handle(Geom_BSplineCurve) spline = GeomAPI_PointsToBSpline(
        pointsArray, 
        Geom_BSplineCurve::MaxDegree() - 6,
        Geom_BSplineCurve::MaxDegree(),
        GeomAbs_C2, 
        Precision::Confusion()).Curve();

    // This one works around a bug in OpenCascade if a curve is closed and
    // periodic. After calling this method, the curve is still closed but
    // no longer periodic, which leads to errors when creating the 3d-lofts
    // from the curves.
    spline->SetNotPeriodic();
    edge = BRepBuilderAPI_MakeEdge(spline);

    // Get control points
    /*
    TColgp_Array1OfPnt poles = spline->Poles();
    for (Standard_Integer i = poles.Lower(); i <= poles.Upper(); ++i) 
    {
      gp_Pnt pt = poles.Value(i);
      std::cout << "Control Point " << i << ": " << pt.X() << ", " << pt.Y() << ", " << pt.Z() << std::endl;
    }*/

    // Get knots
    /*
    TColStd_Array1OfReal knots = spline->Knots();
    for (Standard_Integer i = knots.Lower(); i <= knots.Upper(); ++i) 
    {
      std::cout << "Knot " << i << ": " << knots.Value(i) << std::endl;
    }*/

    // Get multiplicities
    /*
    TColStd_Array1OfInteger mults = spline->Multiplicities();
    for (Standard_Integer i = mults.Lower(); i <= mults.Upper(); ++i) 
    {
      std::cout << "Multiplicity " << i << ": " << mults.Value(i) << std::endl;
    }*/
  }
  else if (method == 3)
  {
    // a terminer: modele pour imposer les points et les tangentes
    Handle(TColgp_HArray1OfPnt) points = new TColgp_HArray1OfPnt(1, 3);
    points->SetValue(1, gp_Pnt(0, 0, 0));      // Fixed point
    points->SetValue(2, gp_Pnt(5, 5, 0));      // Tangency constraint
    points->SetValue(3, gp_Pnt(10, 0, 0));     // Free point

    TColgp_Array1OfVec tangents(1, 3);
    tangents.SetValue(1, gp_Vec(1, 0, 0));     // Tangent at fixed point
    tangents.SetValue(2, gp_Vec(0, 1, 0));     // Tangent at middle point
    tangents.SetValue(3, gp_Vec(0, 0, 0));     // No constraint

    Handle(TColStd_HArray1OfBoolean) tangentFlags = new TColStd_HArray1OfBoolean(1, 3);
    tangentFlags->SetValue(1, Standard_True);  // Enforce tangent
    tangentFlags->SetValue(2, Standard_True);  // Enforce tangent
    tangentFlags->SetValue(3, Standard_False); // No constraint

    GeomAPI_Interpolate interpolator(points, Standard_False, 1e-6);
    interpolator.Load(tangents, tangentFlags, Standard_True);
    interpolator.Perform();
    Handle(Geom_BSplineCurve) constrainedSpline = interpolator.Curve();
  }

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
  builder.Add(compound, edge);

  // export
  delete shape;
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
  SETSHAPE(newshp);

  printf("INFO: after addSpline: Nb edges=%d\n", se->Extent());
  printf("INFO: after addSpline: Nb faces=%d\n", sf->Extent());
  
  RELEASESHAREDN(opc, pc);

  Py_INCREF(Py_None);
  return Py_None;
}

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
#include "Geom_BSplineCurve.hxx"
#include "GeomAPI_PointsToBSpline.hxx"
#include "TColStd_HArray1OfBoolean.hxx"
#include "TColgp_Array1OfVec.hxx"
#include "GeomAPI_Interpolate.hxx"

//=====================================================================
// Add a spline to CAD hook
// method=0; chord length parametrization from control points
// method=1; interpolation of given points
// method=2; uniform parametrization
// degree: degree of spline
//=====================================================================
PyObject* K_OCC::addSpline(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* opc; E_Int method = 0; E_Int degree = 3;
  if (!PYPARSETUPLE_(args, OO_ II_, &hook, &opc, &method, &degree)) return NULL;

  printf("INFO: addSpline: method=%d, degree=%d\n", method, degree);

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

  if (method == 0) // chord length parametrization
  {
    // compute total length of control points
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

    if (ncp == 2) degree = 1;
    else if (ncp == 3) degree = 2;

    E_Int nk = ncp - degree +1;

    // compute knots
    TColStd_Array1OfReal knots(1, nk);
    TColStd_Array1OfInteger mults(1, nk);
  
    E_Float lastKnot = 0.;
    lastP.SetCoord(x[0], y[0], z[0]);
    for (E_Int i = 1; i <= nk; i++)
    {
      if (i == 1 || i == nk) mults.SetValue(i, degree+1);
      else mults.SetValue(i, 1);

      E_Float knot = cp.Value(i).Distance(lastP)/totalLen + lastKnot;
      knots.SetValue(i, knot);

      printf(SF_D_ ": knot=%f, mult=%d\n", i, knots(i), mults(i));

      lastKnot = knot;
      lastP = cp.Value(i);
    }
    Handle(Geom_BSplineCurve) spline = new Geom_BSplineCurve(cp, knots, mults, degree, false);
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

    // Get control points back
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
  else if (method == 2) // quasi uniform
  {
    TColgp_Array1OfPnt cp(1, ncp);
    for (E_Int i = 1; i <= ncp; i++)
    {
      gp_Pnt p(x[i-1],y[i-1],z[i-1]);
      cp.SetValue(i, p);
    }

    // compute knots
    if (ncp == 2) degree = 1;
    else if (ncp == 3) degree = 2;

    E_Int nk = ncp - degree +1;

    TColStd_Array1OfReal knots(1, nk);
    TColStd_Array1OfInteger mults(1, nk);

    for (E_Int i = 1; i <= nk; i++)
    {
      knots(i) = (i-1)*1./(nk-1);
      mults(i) = 1; 
    }
    mults(1) = degree+1;
    mults(nk) = degree+1;
    
    //for (E_Int i = 1; i <= nk; i++)
    //  printf("%d: knot=%g, mult=%d\n", i, knots(i), mults(i));

    Handle(Geom_BSplineCurve) spline =
      new Geom_BSplineCurve(cp, knots, mults, degree);
    edge = BRepBuilderAPI_MakeEdge(spline);
  }
  else
  {
    RELEASESHAREDN(opc, pc);
    PyErr_SetString(PyExc_TypeError,
                    "addSpline: invalid method.");
    return NULL;
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

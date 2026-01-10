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
#include "Geom_BSplineSurface.hxx"
#include "GeomAPI_PointsToBSpline.hxx"
#include "TColStd_HArray1OfBoolean.hxx"
#include "TColgp_Array1OfVec.hxx"
#include "TColgp_Array2OfPnt.hxx"
#include "GeomAPI_Interpolate.hxx"

//=====================================================================
// Add a spline surface to CAD hook
// method=0; chord length parametrization from control points
// method=1; interpolation of given points
// method=2; uniform parametrization
// degree: degree of spline
//=====================================================================
PyObject* K_OCC::addSplineSurface(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* opc; E_Int method = 0; E_Int degree = 3;
  if (!PYPARSETUPLE_(args, OO_ II_, &hook, &opc, &method, &degree)) return NULL;

  printf("INFO: addSplineSurface: method=%d, degree=%d\n", method, degree);

  GETSHAPE;
  GETMAPSURFACES;
  GETMAPEDGES;

  /* get control points (method0) or through points (method1) */
  E_Int im, jm, km;
  FldArrayF* pc; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray3(opc, varString, pc, im, jm, km, cn, eltType);
  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addSplineSurface: invalid array.");
    return NULL;
  }
  if (ret == 2)
  {
    RELEASESHAREDU(opc, pc, cn);
    PyErr_SetString(PyExc_TypeError,
                    "addSplineSurface: array must be structured.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDS(opc, pc);
    PyErr_SetString(PyExc_ValueError,
                    "addSplineSurface: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float* x = pc->begin(posx);
  E_Float* y = pc->begin(posy);
  E_Float* z = pc->begin(posz);
  
  TopoDS_Face face;

  if (method == 2) // quasi uniform
  {
    E_Int ind;
    TColgp_Array2OfPnt cp(1, im, 1, jm);
    for (E_Int j = 1; j <= jm; j++)
    {
      for (E_Int i = 1; i <= im; i++)
      {
        ind = (i-1) + (j-1)*im;
        gp_Pnt p(x[ind],y[ind],z[ind]);
        cp.SetValue(i, j, p);
      }
    }

    // compute knots
    E_Int udegree = degree;
    E_Int vdegree = degree;
    if (im == 2) udegree = 1;
    else if (im == 3) udegree = 2;
    if (jm == 2) vdegree = 1;
    else if (jm == 3) vdegree = 2;

    E_Int nk = im - udegree +1;
    TColStd_Array1OfReal uknots(1, nk);
    TColStd_Array1OfInteger umults(1, nk);
    for (E_Int i = 1; i <= nk; i++)
    {
      uknots(i) = (i-1)*1./(nk-1);
      umults(i) = 1; 
    }
    umults(1) = udegree+1;
    umults(nk) = udegree+1;
    
    nk = jm - vdegree +1;
    TColStd_Array1OfReal vknots(1, nk);
    TColStd_Array1OfInteger vmults(1, nk);
    for (E_Int i = 1; i <= nk; i++)
    {
      vknots(i) = (i-1)*1./(nk-1);
      vmults(i) = 1; 
    }
    vmults(1) = vdegree+1;
    vmults(nk) = vdegree+1;

    Handle(Geom_BSplineSurface) spline =
      new Geom_BSplineSurface(cp, uknots, vknots, umults, vmults, udegree, vdegree);
    face = BRepBuilderAPI_MakeFace(spline, Precision::Confusion());
  }
  else
  {
    RELEASESHAREDN(opc, pc);
    PyErr_SetString(PyExc_TypeError,
                    "addSplineSurface: invalid method.");
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
  builder.Add(compound, face);
  for (E_Int i = 1; i <= edges.Extent(); i++)
  {
    TopoDS_Edge E = TopoDS::Edge(edges(i));
    builder.Add(compound, E);
  }
  
  // export
  delete shape;
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);

  SETSHAPE(newshp);

  printf("INFO: after addSplineSurface: Nb edges=%d\n", se->Extent());
  printf("INFO: after addSplineSurface: Nb faces=%d\n", sf->Extent());
  
  RELEASESHAREDS(opc, pc);

  Py_INCREF(Py_None);
  return Py_None;
}

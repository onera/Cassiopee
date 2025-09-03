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

//=====================================================================
// Add a spline to CAD hook
//=====================================================================
PyObject* K_OCC::addSpline(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* opc;
  if (!PYPARSETUPLE_(args, OO_, &hook, &opc)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  //TopoDS_Shape* shp = (TopoDS_Shape*) packet[0];
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];

  /* generate knots */
  K_FLD::FldArrayF* pc;
  E_Int ret = K_NUMPY::getFromNumpyArray(opc, pc);
  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addSpline: invalid array.");
    return NULL;
  }
  E_Float* x = pc->begin(1);
  E_Float* y = pc->begin(1);
  E_Float* z = pc->begin(1);
  
  // incoming code
  E_Int ncp = pc->getSize();
  std::vector<gp_Pnt>::const_iterator iter;

  // compute total length and copy control points
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
  
  E_Float lastKnot = 0;
  lastP.SetCoord(x[0], y[0], z[0]);
  for (E_Int i = 1; i <= ncp; ++i) 
  {
    if (i == 1 || i == ncp) 
    {
      mults.SetValue(i, 2);
    }
    else 
    {
      mults.SetValue(i, 1);
    }

    E_Float knot = cp.Value(i).Distance(lastP)/totalLen + lastKnot;
    knots.SetValue(i, knot);

    lastKnot = knot;
    lastP = cp.Value(i);
  }

  Handle(Geom_BSplineCurve) spline = new Geom_BSplineCurve(cp, knots, mults, 1, false);

  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(spline);
  
  // Rebuild a single compound
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
    
  for (E_Int i = 1; i <= surfaces.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    builder.Add(compound, F);
  }
  builder.Add(compound, edge);

  // export
  TopoDS_Shape* newshp = new TopoDS_Shape(compound);
    
  // Rebuild the hook
  packet[0] = newshp;
  // Extract surfaces
  TopTools_IndexedMapOfShape* ptr = (TopTools_IndexedMapOfShape*)packet[1];
  delete ptr;
  TopTools_IndexedMapOfShape* sf = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_FACE, *sf);
  packet[1] = sf;

  // Extract edges
  TopTools_IndexedMapOfShape* ptr2 = (TopTools_IndexedMapOfShape*)packet[2];
  delete ptr2;
  TopTools_IndexedMapOfShape* se = new TopTools_IndexedMapOfShape();
  TopExp::MapShapes(*newshp, TopAbs_EDGE, *se);
  packet[2] = se;
  printf("INFO: after addSpline: Nb edges=%d\n", se->Extent());
  printf("INFO: after addSpline: Nb faces=%d\n", sf->Extent());
  
  RELEASESHAREDN(opc, pc);

  Py_INCREF(Py_None);
  return Py_None;

}

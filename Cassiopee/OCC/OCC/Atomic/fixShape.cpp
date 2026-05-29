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
// fixing shape

#include "occ.h"
#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "ShapeUpgrade_UnifySameDomain.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"
#include "BRepLib.hxx"
#include "ShapeFix_ShapeTolerance.hxx"

#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Solid.hxx>

//=====================================================================
// Fix the full shape
// Unify edges
//=====================================================================
PyObject* K_OCC::fixShape(PyObject* self, PyObject* args)
{
  PyObject* hook;
  E_Int unify; E_Int fixShape; E_Int fixWires;
  E_Float tol;
  if (!PYPARSETUPLE_(args, O_ III_ R_, &hook, &fixShape, &fixWires, &unify, &tol)) return NULL;

  GETPACKET;
  GETSHAPE;
  
  TopoDS_Shape* newshp;

  if (fixShape == 0 && fixWires == 0 && unify == 0)
  {
    Py_INCREF(Py_None);
    return Py_None;
  }

  // seems to have no effect
  const ShapeFix_ShapeTolerance& tolFix = ShapeFix_ShapeTolerance();
  tolFix.SetTolerance(*shape, tol, TopAbs_SHAPE);
  BRepLib::UpdateTolerances(*shape, tol);

  if (fixWires == 1)
  {
    // Fix wireframes
    Handle(ShapeFix_Wireframe) sfwf = new ShapeFix_Wireframe(*shape);
    sfwf->SetPrecision(tol*0.001);
    sfwf->SetMaxTolerance(tol);
    sfwf->FixSmallEdges();
    sfwf->FixWireGaps();
    newshp = new TopoDS_Shape();
    *newshp = sfwf->Shape();
    delete shape;
    shape = newshp;
  }

  if (fixShape == 1)
  {
    // Fix shape
    Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape; 
    sfs->Init(*shape); 
    sfs->SetPrecision(1.e-10);
    sfs->SetMaxTolerance(tol*10);
    sfs->SetMinTolerance(tol*0.1);
    sfs->Perform(); 
    newshp = new TopoDS_Shape();
    *newshp = sfs->Shape();
    delete shape;
    shape = newshp;
  }

  // unify edges / faces
  if (unify == 1)
  {
    ShapeUpgrade_UnifySameDomain usd(*shape, Standard_True, Standard_True, Standard_True); // UnifyFaces mode on, UnifyEdges mode on, ConcatBSplines mode on.
    usd.Build();
    newshp = new TopoDS_Shape();
    *newshp = usd.Shape();
    delete shape;
  }

  SETSHAPE(newshp);

  printf("INFO: after fixShape: Nb edges=%d\n", se->Extent());
  printf("INFO: after fixShape: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

  /*
  code to be tested to improve robustness
  // 2) Fix edges
  {
    TopExp_Explorer exp(*shape, TopAbs_EDGE);
    for (; exp.More(); exp.Next())
    {
      const TopoDS_Edge& E = TopoDS::Edge(exp.Current());
      ShapeFix_Edge edgeFix;
      edgeFix.Init(E);
      edgeFix.Perform();          // fix 3D/2D curves, SameParameter, etc.
      // edgeFix.FixReversed2d();  // optional
    }
  }

  // 3) Fix wires
  {
    TopExp_Explorer exp(*shape, TopAbs_WIRE);
    for (; exp.More(); exp.Next())
    {
      const TopoDS_Wire& W = TopoDS::Wire(exp.Current());
      ShapeFix_Wire wireFix;
      wireFix.Load(W);
      wireFix.FixReorder();
      wireFix.FixConnected();
      wireFix.FixEdgeCurves();
      wireFix.FixDegenerated();
      // You can retrieve the fixed wire with wireFix.Wire()
    }
  }

  // 4) Fix faces
  {
    TopExp_Explorer exp(*shape, TopAbs_FACE);
    for (; exp.More(); exp.Next())
    {
      const TopoDS_Face& F = TopoDS::Face(exp.Current());
      ShapeFix_Face faceFix(F);
      faceFix.Perform();
      // faceFix.FixOrientation(); // optional
    }
  }

  // 5) Fix shells and solids
  {
    TopExp_Explorer expShell(*shape, TopAbs_SHELL);
    for (; expShell.More(); expShell.Next())
    {
      const TopoDS_Shell& Sh = TopoDS::Shell(expShell.Current());
      ShapeFix_Shell shellFix(Sh);
      shellFix.Perform();
    }

    TopExp_Explorer expSolid(*shape, TopAbs_SOLID);
    for (; expSolid.More(); expSolid.Next())
    {
      const TopoDS_Solid& So = TopoDS::Solid(expSolid.Current());
      ShapeFix_Solid solidFix(So);
      solidFix.Perform();
    }
  }

  // 7) Global shape fix (wrap-up)
  {
    ShapeFix_Shape sfix(*shape);
    sfix.Perform();
    shape = sfix.Shape();
  }

  // 8) Unify same domain (merge edges/faces)
  {
    ShapeUpgrade_UnifySameDomain unify(shape,
                                        Standard_True,  // unify edges
                                        Standard_True,  // unify faces
                                        Standard_True); // concat BSplines, etc.
    unify.Build();
    shape = unify.Shape();
  }*/

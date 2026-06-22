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
#include "BRep_Builder.hxx"
#include "ShapeFix_ShapeTolerance.hxx"
#include <StdFail_NotDone.hxx>

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
  PyObject* hook; PyObject* listFaces;
  E_Int fixShape, fixWires;
  E_Int unifyEdges, unifyFaces; 
  E_Float tol; // use in fixShapes
  E_Float linearDeflection; // used in unify
  E_Float angularDeflection; // used in unify, in degrees
  if (!PYPARSETUPLE_(args, O_ O_ IIII_ RRR_, 
    &hook, &listFaces,
    &fixShape, &fixWires, 
    &unifyEdges, &unifyFaces, 
    &tol, &linearDeflection, &angularDeflection)) return NULL;

  GETPACKET;
  GETSHAPE;
  GETMAPSURFACES;
  
  TopoDS_Shape* newshp;

  if (fixShape == 0 && fixWires == 0 && unifyEdges == 0 && unifyFaces == 0)
  {
    Py_INCREF(Py_None);
    return Py_None;
  }

  if (listFaces != Py_None)
  {
    printf("Warning: fixShape: not implemented for face list.");
    Py_INCREF(Py_None);
    return Py_None;
  }

  // Build compound from face list
  TopoDS_Compound compound; // face compound
  TopoDS_Compound compound2; // remaining faces of shape
  if (listFaces != Py_None)
  {
    E_Int nfaces = PyList_Size(listFaces);
    BRep_Builder builder;
    builder.MakeCompound(compound);
    for (E_Int i = 0; i < nfaces; i++)
    {
      PyObject* noO = PyList_GetItem(listFaces, i);
      E_Int no = PyInt_AsLong(noO);
      const TopoDS_Face& F = TopoDS::Face(surfaces(no));
      builder.Add(compound, F);
    }
  }
  else 
  {
    BRep_Builder builder;
    builder.MakeCompound(compound);
    builder.Add(compound, *shape);
  }

  // seems to have no effect
  const ShapeFix_ShapeTolerance& tolFix = ShapeFix_ShapeTolerance();
  tolFix.SetTolerance(compound, tol, TopAbs_SHAPE);
  BRepLib::UpdateTolerances(compound, tol);

  if (fixWires == 1)
  {
    // Fix wireframes
    Handle(ShapeFix_Wireframe) sfwf = new ShapeFix_Wireframe(compound);
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
    sfs->Init(compound); 
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
  if (unifyEdges == 1 || unifyFaces == 1)
  {
    E_Bool concatSplines = false;
    if (unifyEdges == 1) concatSplines = true;    
    ShapeUpgrade_UnifySameDomain usd(compound, unifyFaces, unifyEdges, concatSplines);
    usd.AllowInternalEdges(Standard_True);
    usd.SetLinearTolerance(linearDeflection);
    usd.SetAngularTolerance(angularDeflection*K_CONST::E_PI/180.);
    try 
    {
      usd.Build();
    } catch (StdFail_NotDone& e) 
    { 
      PyErr_SetString(PyExc_TypeError, "fixShape: fail to unify (notDone).");  
      Py_INCREF(Py_None);
      return Py_None;
    }
    newshp = new TopoDS_Shape();
    *newshp = usd.Shape();
    delete shape;
    shape = newshp;
  }

  SETSHAPE(shape);

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

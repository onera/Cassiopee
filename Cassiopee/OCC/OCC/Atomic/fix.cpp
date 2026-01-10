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
// CAD fix

#include "occ.h"
#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "ShapeUpgrade_UnifySameDomain.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "TopoDS.hxx"

//=====================================================================
// Fix the full shape
// Unify edges
//=====================================================================
PyObject* K_OCC::fixShape(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];

  TopoDS_Shape* newshp = new TopoDS_Shape();

  // Build a shape fixer
  /*
  Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape; 
  sfs->Init(*shp); 
  sfs->SetPrecision(1.e-10); 
  sfs->SetMaxTolerance(1.e-6); 
  sfs->SetMinTolerance(1.e-8);
  sfs->Perform(); 
  *newshp = sfs->Shape(); 
  */

  // Fix wireframes
  /*
  Handle(ShapeFix_Wireframe) SFWF = new ShapeFix_Wireframe(*shp);
  SFWF->SetPrecision(1.e-10);
  SFWF->SetMaxTolerance(1.e-6);
  SFWF->FixSmallEdges();
  SFWF->FixWireGaps();
  *newshp = SFWF->Shape(); 
  */

  // Suppress two faces and rebuild shape
  /*
  ShapeBuild_ReShape builder;
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  E_Int nbFaces = surfaces.Extent();
  const TopoDS_Face& F1 = TopoDS::Face(surfaces(1));
  const TopoDS_Face& F2 = TopoDS::Face(surfaces(2));
  builder.Remove(F1);
  builder.Remove(F2);
  *newshp = builder.Apply(*shp, TopAbs_SHAPE);
  */

  // unify edges / faces
  ShapeUpgrade_UnifySameDomain USD(*shp, Standard_True, Standard_True, Standard_False); // UnifyFaces mode on, UnifyEdges mode on, ConcatBSplines mode off.
  USD.Build();
  *newshp = USD.Shape();

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
  printf("INFO: after fix: Nb edges=%d\n", se->Extent());
  printf("INFO: after fix: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

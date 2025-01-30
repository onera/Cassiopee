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
// CAD sewing
#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "ShapeAnalysis.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "ShapeUpgrade_FaceDivide.hxx"
#include "ShapeUpgrade_ShapeDivideArea.hxx"
#include "ShapeUpgrade_ShapeDivideClosed.hxx"
#include "ShapeUpgrade_ClosedFaceDivide.hxx"
#include "ShapeUpgrade_SplitSurfaceArea.hxx"
#include "TColGeom_SequenceOfSurface.hxx"
#include "ShapeExtend_CompositeSurface.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "BRep_Builder.hxx"
#include "ShapeUpgrade_ShapeDivideClosedEdges.hxx"
#include "BRepBuilderAPI_Sewing.hxx"

//=====================================================================
// Sew faces removing extra edges
//=====================================================================
PyObject* K_OCC::sewing(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* listFaces; E_Float tol;
  if (!PYPARSETUPLE_(args, OO_ R_, &hook, &listFaces, &tol)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  //TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];

  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];
  const Standard_Real tolerance = tol;
  BRepBuilderAPI_Sewing sewer(tolerance);
  
  TopoDS_Shape* newshp = NULL;
  E_Int nfaces = PyList_Size(listFaces);
  //nfaces = 0; // force car le code par subfaces semble ne pas marcher
  if (nfaces == 0)
  {
    // top shape
    TopoDS_Shape shc;
    printf("Info: sewing top shape.\n");
    sewer.Add(*shp);
    sewer.Perform();
    shc = sewer.SewedShape();
    newshp = new TopoDS_Shape(shc);
  }
  else
  {
    // Build remaining faces list
    std::list<E_Int> pl;
    E_Int nf = surfaces.Extent();
    printf("Info: sewing %d / %d faces.\n", nfaces, nf);
    for (E_Int i = 1; i <= nf; i++) pl.push_back(i);

    for (E_Int no = 0; no < PyList_Size(listFaces); no++)
    {
      PyObject* noFaceO = PyList_GetItem(listFaces, no);
      E_Int noFace = PyInt_AsLong(noFaceO);
      auto it = std::find(pl.begin(), pl.end(), noFace);
      if (it != pl.end()) pl.erase(it);
      const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
      sewer.Add(F);
    } 
    sewer.Perform();
    TopoDS_Shape shs = sewer.SewedShape();

    /*
    ShapeBuild_ReShape reshaper;
    TopTools_IndexedMapOfShape faces;
    TopExp::MapShapes(shs, TopAbs_FACE, faces);

    for (E_Int i = 0; i < faces.Extent(); i++) 
    {
      PyObject* noFaceO = PyList_GetItem(listFaces, i);
      E_Int noFace = PyInt_AsLong(noFaceO);
      const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));
      const TopoDS_Face& F2 = TopoDS::Face(faces(i+1));
      printf("replace face %d with %d\n", noFace, i+1);
      reshaper.Replace(F, F2);
    }
    shc = reshaper.Apply(*shp);
    */

    BRep_Builder builder;
    TopoDS_Compound shc;
    builder.MakeCompound(shc);
    for (auto& i : pl)
    {
      TopoDS_Face F = TopoDS::Face(surfaces(i));
      builder.Add(shc, F);
    }

    TopExp_Explorer expl1(shs, TopAbs_FACE);
    while (expl1.More())
    {
      TopoDS_Shape shape = expl1.Current();
      TopoDS_Face face = TopoDS::Face(shape);
      builder.Add(shc, face);
      expl1.Next();
    }
    newshp = new TopoDS_Shape(shc);
  }

  // export
  delete shp;
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
  printf("INFO: after sewing: Nb edges=%d\n", se->Extent());
  printf("INFO: after sewing: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}

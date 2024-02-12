/*    
    Copyright 2013-2024 Onera.

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
// CAD split for parallel
#include "occ.h"

#include "TopoDS.hxx"
#include "BRep_Tool.hxx"
#include "TopoDS_Wire.hxx"
#include "BRepCheck_Wire.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "ShapeAnalysis.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GCPnts_AbscissaPoint.hxx" 
#include "GCPnts_UniformDeflection.hxx"
#include "GCPnts_UniformAbscissa.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "Geom2d_Curve.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"
#include "BRepGProp.hxx"
#include "ShapeUpgrade_FaceDivide.hxx"
#include "ShapeUpgrade_ShapeDivideArea.hxx"
#include "ShapeUpgrade_ShapeDivideClosed.hxx"
#include "ShapeUpgrade_ClosedFaceDivide.hxx"
#include "ShapeUpgrade_SplitSurfaceArea.hxx"
#include "TColGeom_SequenceOfSurface.hxx"
#include "ShapeExtend_CompositeSurface.hxx"

#include "GProp_GProps.hxx"

//=====================================================================
// Split every faces
//=====================================================================
PyObject* K_OCC::splitFaces(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PYPARSETUPLE_(args, O_, &hook)) return NULL;

    void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];


  // try on all shape
  TopoDS_Shape* shp = (TopoDS_Shape*)packet[0];
  
  printf("splitting top shape\n");

  // Closed Shape Divide
  //ShapeUpgrade_ShapeDivideClosed splitter(*shp);
  //splitter.SetNbSplitPoints(20); // cree num+1 faces
  //splitter.Perform(Standard_False);

  // Shape Divide
  //ShapeUpgrade_ShapeDivide splitter(*shp);
  //splitter.SetNbSplitPoints(20); // cree num+1 faces
  //splitter.Perform(Standard_False);

  // Shape Divide by Area (not ok in OCC 7.6)
  //ShapeUpgrade_ShapeDivideArea splitter(*shp);
  //splitter.MaxArea();
  //splitter.Perform(Standard_False);
  

  TopExp_Explorer expl;
  E_Int nbFaces = surfaces.Extent();
  printf("number of faces=%d\n", nbFaces);

  for (E_Int noFace = 1; noFace <= nbFaces; noFace++)
  {
    const TopoDS_Face& F = TopoDS::Face(surfaces(noFace));

    // Closed Face Divide
    //ShapeUpgrade_ClosedFaceDivide splitter(F);
    //splitter.SetNbSplitPoints(5); // cree num+1 faces
    //splitter.SplitSurface();
    //Standard_Integer status = splitter.Status(ShapeExtend_DONE2);
    //printf("status=%d\n", status);
    //TopoDS_Shape result = splitter.Result();
    //TopTools_IndexedMapOfShape* surfs = new TopTools_IndexedMapOfShape();
    //TopExp::MapShapes(result, TopAbs_FACE, *surfs);
    //TopTools_IndexedMapOfShape& surfaces2 = *surfs;
    //printf("splitted in %d faces.\n", surfaces2.Extent());

    //ShapeUpgrade_FaceDivide splitter(F);
    //splitter.SetNbSplitPoints(2); // split in 2
    //splitter.Perform(Standard_True);
    //TopoDS_Shape result2 = splitter.Result();
        
    // Split a face with a wire (snippet)
    //TopoDS_Shape box = BRepPrimAPI_MakeBox(10, 10, 10).Shape();
    //gp_Circ circle(gp_Ax2(gp_Pnt(5, 5, 5), gp_Dir(0, 0, 1)), 3);
    //TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(circle).Edge();
    //TopoDS_Wire wire = BRepBuilderAPI_MakeWire(edge).Wire();
    //BRepFeat_SplitShape splitter(box);
    //splitter.Add(wire, box.Face(6));
    //splitter.Build();
    //TopoDS_Shape result = splitter.Shape();

    // Split a face in UV (snippet)
    // Create a surface to split
    Handle(Geom_Surface) aSurf = BRep_Tool::Surface(F);
    // Create a splitter with default constructor
    ShapeUpgrade_SplitSurfaceArea aSplitter;
    // Set the number of parts to split the surface into
    aSplitter.NbParts() = 4;
    // Set the splitting mode to true (split into squares)
    //aSplitter.SetSplittingIntoSquares(Standard_True);
    // Initialize the splitter with the surface
    aSplitter.Init(aSurf);
    // Perform the splitting
    aSplitter.Perform();
    // Get the resulting surfaces as a sequence
    Handle(ShapeExtend_CompositeSurface) aResSurfs = aSplitter.ResSurfaces();
    // Comment reshaper?
    // ShapeBuild_ReShape use replace, remove


  }

  // repush and update shape
  

  Py_INCREF(Py_None);
  return Py_None;
}

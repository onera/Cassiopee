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
// split face from two vertices
#include "occ.h"

#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Vertex.hxx"
#include "TopoDS_Edge.hxx"
#include "BRep_Tool.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "Geom_Surface.hxx"
#include "Geom2d_Curve.hxx"
#include "GCE2d_MakeSegment.hxx"
#include "gp_Pnt2d.hxx"
#include "TopLoc_Location.hxx"
#include "TopoDS_Shape.hxx"
#include "BRepFeat_SplitShape.hxx"
#include "Standard_Transient.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRep_Builder.hxx"
#include "TopExp.hxx"
#include "TopExp_Explorer.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "GeomAPI_ProjectPointOnSurf.hxx"
#include "BRepAlgoAPI_Splitter.hxx"
#include "BRepCheck_Analyzer.hxx"
#include "gp_Pln.hxx"
#include "BRepBuilderAPI_MakeFace.hxx"

//=========================================================================
// Creates an edge on a face's surface connecting two vertices via a 
// straight parametric (UV) line.
//=========================================================================
TopoDS_Edge createParametricEdgeBetweenVertices(const TopoDS_Face& face,
                                                //const TopoDS_Vertex& v1,
                                                //const TopoDS_Vertex& v2,
                                                E_Float u1, E_Float v1,
                                                E_Float u2, E_Float v2)
{
  // 1. Get the UV coordinates of the two vertices on the given face
  //Standard_Real pu1, pv1, pu2, pv2;
  gp_Pnt2d p2d1(u1, v1);
  gp_Pnt2d p2d2(u2, v2);

  //printf("getting parameters\n");
  //gp_Pnt2d p2d1 = BRep_Tool::Parameters(v1, face);
  //printf("uv1 = %g %g\n", p2d1.X(), p2d1.Y());
  //gp_Pnt2d p2d2 = BRep_Tool::Parameters(v2, face);
  //printf("get UV OK\n");

  // 2. Create a 2D line segment in the parameter space (UV space)
  GCE2d_MakeSegment segBuilder(p2d1, p2d2);
  if (!segBuilder.IsDone()) 
  {
    return TopoDS_Edge(); // Return null edge if invalid (e.g., identical points)
  }
  Handle(Geom2d_TrimmedCurve) curve2d = segBuilder.Value();
  printf("trimmed curve ok\n");

  // 3. Retrieve the underlying surface and location of the face
  TopLoc_Location loc;
  Handle(Geom_Surface) surface = BRep_Tool::Surface(face, loc);

  // 4. Build the 3D edge from the 2D curve and the surface
  // This automatically projects the UV line onto the 3D curved surface geometry.
  BRepBuilderAPI_MakeEdge edgeBuilder(curve2d, surface);
  printf("edge on surface ok\n");

  if (!edgeBuilder.IsDone()) 
  {
    return TopoDS_Edge();
  }
  printf("before return\n");

  return edgeBuilder.Edge();
}

//==========================================
// Splits a face with edge
//==========================================
TopoDS_Shape splitFaceByEdge(const TopoDS_Face& face, 
                             const TopoDS_Edge& cutEdge)
{
  //BRepFeat_SplitShape splitter(face);
  //splitter.Add(cutEdge, face);

  TopTools_ListOfShape faces;
  faces.Append(face);
  TopTools_ListOfShape cuttings;
  cuttings.Append(cutEdge);

  BRepAlgoAPI_Splitter splitter;  
  splitter.SetArguments(faces);
  splitter.SetTools(cuttings);
  //splitter.SetRunParallel(Standard_True);
  //splitter.SetFuzzyValue(1e-5);
  splitter.Build();
  
  if (splitter.IsDone()) 
  {
    return splitter.Shape();
  }

  return face;
}

//=====================================================================
// Split face by two points
//=====================================================================
PyObject* K_OCC::splitFace(PyObject* self, PyObject* args)
{
  PyObject* hook; E_Int no; E_Float x1, y1, z1; E_Float x2, y2, z2;
  if (!PYPARSETUPLE_(args, O_ I_ TRRR_ TRRR_, &hook, &no,
                    &x1, &y1, &z1, &x2, &y2, &z2)) return NULL;

  /*{
    // 1. Create a bounded 100x100 planar face on XY plane
    gp_Pln plane(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1));
    TopoDS_Face face = BRepBuilderAPI_MakeFace(plane, 0.0, 100.0, 0.0, 100.0).Face();

    // 2. Create a finite splitting edge lying on that plane
    TopoDS_Edge cutEdge = BRepBuilderAPI_MakeEdge(gp_Pnt(0, 0, 0), gp_Pnt(100, 100, 0)).Edge();

    // 3. Split
    TopTools_ListOfShape args, tools;
    args.Append(face);
    tools.Append(cutEdge);

    BRepAlgoAPI_Splitter splitter;
    splitter.SetArguments(args);
    splitter.SetTools(tools);
    splitter.Build();

    if (!splitter.IsDone()) {
        std::cerr << "Splitter failed!" << std::endl;
        return 1;
    }

    // 4. Count resulting faces (should print 2)
    int faceCount = 0;
    for (TopExp_Explorer exp(splitter.Shape(), TopAbs_FACE); exp.More(); exp.Next()) {
        faceCount++;
    }
    std::cout << "Resulting faces: " << faceCount << std::endl;
    return 0;
  }*/

  GETPACKET;
  GETSHAPE;
  GETMAPSURFACES;

  TopoDS_Face face = TopoDS::Face(surfaces(no)); // existing face

  gp_Pnt p1(x1, y1, z1);
  gp_Pnt p2(x2, y2, z2);
  printf("P1 %g %g %g\n",x1,y1,z1);
  printf("P2 %g %g %g\n",x2,y2,z2);
  
  TopoDS_Vertex vertex1 = BRepBuilderAPI_MakeVertex(p1);
  TopoDS_Vertex vertex2 = BRepBuilderAPI_MakeVertex(p2);

  // project point on face
  BRepExtrema_DistShapeShape extrema(vertex1, face);
  if (extrema.IsDone() && extrema.NbSolution() > 0)
  {
    p1 = extrema.PointOnShape2(1);
    vertex1 = BRepBuilderAPI_MakeVertex(p1);
    printf("projected1: %g %g %g\n", p1.X(),p1.Y(),p1.Z()); 
  }
  BRepExtrema_DistShapeShape extrema2(vertex2, face);
  if (extrema2.IsDone() && extrema2.NbSolution() > 0) 
  {
    p2 = extrema2.PointOnShape2(1);
    vertex2 = BRepBuilderAPI_MakeVertex(p2);
    printf("projected2: %g %g %g\n", p2.X(),p2.Y(),p2.Z()); 
  }
  printf("projection of points on face.\n");

  Standard_Real u1, v1, u2, v2;
  Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
  GeomAPI_ProjectPointOnSurf projector1(p1, surface);
  projector1.LowerDistanceParameters(u1, v1);
  GeomAPI_ProjectPointOnSurf projector2(p2, surface);
  projector2.LowerDistanceParameters(u2, v2);
  printf("params %g %g / %g %g\n", u1, v1, u2, v2);

  // create cutting edge
  //TopoDS_Edge cutEdge = createParametricEdgeBetweenVertices(face, u1, v1, u2, v2);
  // DEBUG
  TopoDS_Edge cutEdge = BRepBuilderAPI_MakeEdge(p1, p2).Edge();

  // split face with cutting edge
  TopoDS_Shape sh = splitFaceByEdge(face, cutEdge);
  // DEBUG
  //TopoDS_Shape sh;

  printf("split done.\n");

  // rebuild compound
  BRep_Builder builder;
  TopoDS_Compound compound;
  builder.MakeCompound(compound);
  for (E_Int i = 1; i <= surfaces.Extent(); i++)
  {
    TopoDS_Face F = TopoDS::Face(surfaces(i));
    if (i != no) builder.Add(compound, F);
    else
    {
      for (TopExp_Explorer explorer(sh, TopAbs_FACE); explorer.More(); explorer.Next()) 
      {
        TopoDS_Face F2 = TopoDS::Face(explorer.Current());
        builder.Add(compound, F2);
      }
    }
  }
  builder.Add(compound, cutEdge);

  TopoDS_Shape* newshp = new TopoDS_Shape(compound);

  delete shape;
  SETSHAPE(newshp);

  printf("INFO: after splitFace: Nb edges=%d\n", se->Extent());
  printf("INFO: after splitFace: Nb faces=%d\n", sf->Extent());

  Py_INCREF(Py_None);
  return Py_None;
}
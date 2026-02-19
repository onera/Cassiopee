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
#include "OCCSurface.h"
#include "TopoDS.hxx"
#include "TopoDS_Edge.hxx"
#include "TopoDS_Face.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "Poly_Triangulation.hxx"
#include "BRep_Tool.hxx"
#include "Standard_Version.hxx"

//===========================================================
// occmesh: mesh with TRI on a CAD shape
// use occ mesher based only on hausd
// exports faces and edges
//===========================================================
PyObject* K_OCC::occmesh(PyObject* self, PyObject* args)
{
  PyObject* hook;
  E_Float hausd;
  E_Float angularDeflection; // in degrees
  if (!PYPARSETUPLE_(args, O_ R_ R_, &hook, &hausd, &angularDeflection)) return NULL;  
  
  GETSHAPE;
  
  // Triangulate shape
  E_Float angle = angularDeflection * M_PI / 180.; // en radians
  Standard_Boolean relative = Standard_False;
  if (hausd < 0) { relative = Standard_True; hausd = -hausd; }
  BRepMesh_IncrementalMesh Mesh(*shape, hausd, relative, angle, Standard_True);
  Mesh.Perform(); 

  //===================
  // Recupere les faces
  //===================
  PyObject* dfaces = PyList_New(0);
  for (TopExp_Explorer exp(*shape, TopAbs_FACE); exp.More(); exp.Next())
  {
    TopoDS_Face face = TopoDS::Face(exp.Current());

    // get the triangulation
    TopLoc_Location loc;
    Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
    if (tri.IsNull()) continue;
#if OCC_VERSION_MAJOR >= 7 && OCC_VERSION_MINOR >= 6
    Handle(TColgp_HArray1OfPnt) nodes = tri->MapNodeArray();
#else
    const TColgp_Array1OfPnt* nodes = &(tri->Nodes());
#endif    

    // Get dims
    E_Int nbNodes = 0;
    for (Standard_Integer iCount = nodes->Lower(); iCount <= nodes->Upper(); iCount++)
    {
      nbNodes++;
    }
    E_Int nbTris = tri->NbTriangles();

    // buildArray
    PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", nbNodes, nbTris, "TRI", false, 1);
    FldArrayF* f; FldArrayI* cn;
    K_ARRAY::getFromArray3(o, f, cn);
    E_Float* fx = f->begin(1);
    E_Float* fy = f->begin(2);
    E_Float* fz = f->begin(3);

    // create a transformation from the location
    gp_Trsf xloc = loc;

    // get the nodes
    Standard_Real x,y,z;
    for (Standard_Integer i = nodes->Lower(); i <= nodes->Upper(); i++)
    {
      (*nodes)(i).Coord(x,y,z);
      //tri->Node(iCount);
      xloc.Transforms(x,y,z);
      fx[i-1] = x;
      fy[i-1] = y;
      fz[i-1] = z;
      //if (i > nbNodes) printf("danger i > nbNodes %d %d\n", i, nbNodes);
      //printf("point %d: %f %f %f\n",i,x,y,z);
    }
    
    // copy the TRIs
    Standard_Integer i1, i2, i3;
#if OCC_VERSION_MAJOR < 7
    const Poly_Array1OfTriangle& tris = tri->Triangles();
#endif
    for (Standard_Integer i = 1; i <= tri->NbTriangles(); i++) 
    {
      // get the node indexes for this triangle
#if OCC_VERSION_MAJOR < 7
      Poly_Triangle tril = tris(i);
#else
      Poly_Triangle tril = tri->Triangle(i);
#endif
      tril.Get(i1, i2, i3);
      (*cn)(i-1,1) = i1;
      (*cn)(i-1,2) = i2;
      (*cn)(i-1,3) = i3;
      //if (i > nbTris) printf("danger nbTris %d %d\n", i, nbNodes);
      //if (i1 > nbNodes) printf("danger1 %d %d\n",i1,nbNodes);
      //if (i2 > nbNodes) printf("danger2 %d %d\n",i2,nbNodes);
      //if (i3 > nbNodes) printf("danger3 %d %d\n",i3,nbNodes);
      //printf("TRI %d: %d %d %d\n", i, i1,i2,i3);
    }
    PyList_Append(dfaces, o); Py_DECREF(o);
    RELEASESHAREDU(o, f, cn);
  }

  //====================
  // recupere les edges
  //====================
  PyObject* dedges = PyList_New(0);

  TopExp_Explorer expEdge;
  for (expEdge.Init(*shape, TopAbs_EDGE); expEdge.More(); expEdge.Next()) 
  {
    TopoDS_Edge edge = TopoDS::Edge(expEdge.Current());
    PyObject* o = NULL; FldArrayF* f;

    // 1) Essayer d'obtenir le polygone 3D directement
    TopLoc_Location loc;
    Handle(Poly_Polygon3D) poly3d = BRep_Tool::Polygon3D(edge, loc);
    if (!poly3d.IsNull())
    {
      const TColgp_Array1OfPnt& nodes = poly3d->Nodes();
      std::cout << "Edge discretized points (3D):\n";
      E_Int nbNodes = nodes.Upper()-nodes.Lower()+1;
      o = K_ARRAY::buildArray3(3, "x,y,z", nbNodes, 1, 1, 1);
      K_ARRAY::getFromArray3(o, f);
      E_Float* fx = f->begin(1);
      E_Float* fy = f->begin(2);
      E_Float* fz = f->begin(3);

      for (E_Int i = nodes.Lower(); i <= nodes.Upper(); ++i) 
      {
        gp_Pnt p = nodes(i);
        fx[i-1] = p.X();
        fy[i-1] = p.Y();
        fz[i-1] = p.Z();
        //std::cout << "  " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
      }
    }
    else
    {
      // 2) Sinon, on essaie via PolygonOnTriangulation
      Handle(Poly_PolygonOnTriangulation) polyOnTri;
      Handle(Poly_Triangulation) triangulation;
      BRep_Tool::PolygonOnTriangulation(edge, polyOnTri, triangulation, loc);

      if (!polyOnTri.IsNull() && !triangulation.IsNull())
      {
        const TColStd_Array1OfInteger& indices = polyOnTri->Nodes();

#if OCC_VERSION_MAJOR >= 7 && OCC_VERSION_MINOR >= 6
        Handle(TColgp_HArray1OfPnt) triNodesPt = triangulation->MapNodeArray();
        TColgp_Array1OfPnt& triNodes = *triNodesPt;
#else
        const TColgp_Array1OfPnt& triNodes = triangulation->Nodes();
#endif
        gp_Trsf xloc = loc.Transformation();
        //std::cout << "Edge discretized points (via triangulation):\n";
        E_Int nbNodes = indices.Upper()-indices.Lower()+1;
        o = K_ARRAY::buildArray3(3, "x,y,z", nbNodes, 1, 1, 1);
        K_ARRAY::getFromArray3(o, f);
        E_Float* fx = f->begin(1);
        E_Float* fy = f->begin(2);
        E_Float* fz = f->begin(3);

        for (E_Int i = indices.Lower(); i <= indices.Upper(); ++i) 
        {
          gp_Pnt p = triNodes(indices(i));
          p.Transform(xloc);
          fx[i-1] = p.X();
          fy[i-1] = p.Y();
          fz[i-1] = p.Z();
          //std::cout << "  " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
        }
      }
    }
    if (o != NULL)
    {
      PyList_Append(dedges, o); Py_DECREF(o);
      RELEASESHAREDS(o, f);
    }
  }

  // export both
  PyObject* tpl = PyList_New(0);
  PyList_Append(tpl, dedges); Py_DECREF(dedges);
  PyList_Append(tpl, dfaces); Py_DECREF(dfaces);
  
  return tpl;
}

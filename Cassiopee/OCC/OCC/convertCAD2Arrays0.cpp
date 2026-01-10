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
// convertCAD2Arrays (algo=0)

#define OCEVERSION 0

#include "occ.h"
#include "kcore.h"
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "STEPControl_Reader.hxx"
#include "IGESControl_Reader.hxx"
#include "StlAPI_Writer.hxx"
#include "VrmlAPI_Writer.hxx"
#include "IGESControl_Controller.hxx"
#include "STEPControl_Writer.hxx"

#include "TopoDS_Shape.hxx"
#include "Standard_ConstructionError.hxx"

#include "BRepPrimAPI_MakeCylinder.hxx"
#include "BRepPrimAPI_MakeBox.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepGProp.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "GProp_GProps.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "Poly_Triangulation.hxx"

#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"

#include "Standard_Version.hxx"

// ============================================================================
/* Essai avec utilisation de BrepMesh d'open cascade */
// ============================================================================
PyObject* K_OCC::convertCAD2Arrays0(PyObject* self, PyObject* args)
{ 
  // inFileFormat: fmt_iges ou fmt_step
  // outFileFormat: None ou fmt_stl ou fmt_vrml
  // deflection: deviation pour le maillage
  char* inFileName; char* inFileFormat; 
  char* outFileName; char* outFileFormat;
  E_Float deflection;
  if (!PYPARSETUPLE_(args, SSSS_ R_,
                    &inFileName, &inFileFormat, 
                    &outFileName, &outFileFormat,
                    &deflection)) return NULL;

  // Read CAD
  TopoDS_Shape shape;
  if (strcmp(inFileFormat, "fmt_step") == 0)
  {
    STEPControl_Reader reader;
    reader.ReadFile(inFileName);
    reader.TransferRoots();
    shape = reader.OneShape();
    printf("done reading %s.\n", inFileName);
  }
  if (strcmp(inFileFormat, "fmt_iges") == 0)
  {
    IGESControl_Controller::Init();
    IGESControl_Reader reader;
    reader.ReadFile(inFileName);
    reader.SetReadVisible(Standard_True);
    reader.PrintCheckLoad(Standard_True, IFSelect_GeneralInfo);
    reader.ClearShapes();
    reader.TransferRoots();
    shape = reader.OneShape();
    printf("done reading %s.\n", inFileName); 
  }

  /*
  E_Float precision = deflection*1.e-2;
  E_Float maxTol = deflection*100;
  E_Float minTol = deflection*10;
  printf("%f %f %f\n", precision, maxTol, minTol);
  */

  // FIX cad
  /*
  Handle(ShapeFix_Shape) aFixShape = new ShapeFix_Shape(shape);
  aFixShape->SetPrecision(precision);
  aFixShape->SetMaxTolerance(maxTol);
  aFixShape->SetMinTolerance(minTol);
  aFixShape->Perform();
  shape = aFixShape->Shape();
  */

  // Fix Cad wire
  /*
  Handle(ShapeFix_Wireframe) aFixWire = new ShapeFix_Wireframe(shape);
  aFixWire->SetPrecision(precision);
  aFixWire->SetMaxTolerance(maxTol);
  //aFixWire->DropSmallEdgesMode() = true;
  aFixWire->FixSmallEdges();
  aFixWire->FixWireGaps();
  shape = aFixWire->Shape();
  */
 
  // Triangulate
  E_Float angularDeflection = 0.5; // en degres
  Standard_Boolean relative = Standard_False;
  if (deflection < 0) { relative = Standard_True; deflection = -deflection; }
  BRepMesh_IncrementalMesh Mesh(shape, deflection, relative, angularDeflection, Standard_True);
  
  Mesh.Perform(); 

  PyObject* out = PyList_New(0);
  
  E_Int nbNodes = 0;
  E_Int nbTris = 0;
  
  // Dimensionnement
  for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
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
    for (Standard_Integer iCount = nodes->Lower(); iCount <= nodes->Upper(); iCount++)
    {
      nbNodes++;  
    }
    nbTris += tri->NbTriangles();
  }
  
  printf("INFO: total number of nodes: " SF_D_ "\n", nbNodes);
  printf("INFO: total number of triangles: " SF_D_ "\n", nbTris);
  
  // buildArray
  PyObject* o = K_ARRAY::buildArray3(3, "x,y,z", nbNodes, nbTris, "TRI", false, 1);
  FldArrayF* f; FldArrayI* c;
  K_ARRAY::getFromArray3(o, f, c);
  E_Float* fx = f->begin(1);
  E_Float* fy = f->begin(2);
  E_Float* fz = f->begin(3);
  E_Int* c1 = c->begin(1);
  E_Int* c2 = c->begin(2);
  E_Int* c3 = c->begin(3);
  E_Int stride = c->getStride(); 

  E_Int cNodes = 0; E_Int cTris = 0;
  for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
  {
    TopoDS_Face face = TopoDS::Face(exp.Current());

    // get the triangulation
    TopLoc_Location loc;
    Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);
    if (tri.IsNull()) continue;

    // create a transformation from the location
    gp_Trsf xloc = loc;

#if OCC_VERSION_MAJOR >= 7 && OCC_VERSION_MINOR >= 6
    Handle(TColgp_HArray1OfPnt) nodes = tri->MapNodeArray();
#else
    const TColgp_Array1OfPnt* nodes = &(tri->Nodes());
#endif    
    // get the nodes
    for (Standard_Integer iCount = nodes->Lower(); iCount <= nodes->Upper(); iCount++)
    {
      Standard_Real x,y,z;
      (*nodes)(iCount).Coord(x,y,z);
      //tri->Node(iCount);
      xloc.Transforms(x,y,z);
      fx[cNodes+iCount-1] = x;
      fy[cNodes+iCount-1] = y;
      fz[cNodes+iCount-1] = z;
      //if (iCount > nbNodes) printf("danger iCount>nbNodes %d %d\n", iCount, nbNodes);
      //printf("point %d: %f %f %f\n",iCount,x,y,z);
    }
    
    // copy the polygons
    Standard_Integer i1, i2, i3;
#if OCC_VERSION_MAJOR < 7
    const Poly_Array1OfTriangle &tris = tri->Triangles();
#endif
    for (Standard_Integer iCount = 1; iCount <= tri->NbTriangles(); iCount++) 
    {
      // get the node indexes for this triangle
#if OCC_VERSION_MAJOR < 7
      Poly_Triangle tril = tris(iCount);
#else
      Poly_Triangle tril = tri->Triangle(iCount);
#endif
      tril.Get(i1, i2, i3);
      c1[stride*(cTris+iCount-1)] = i1+cNodes;
      c2[stride*(cTris+iCount-1)] = i2+cNodes;
      c3[stride*(cTris+iCount-1)] = i3+cNodes;
      //if (iCount > nbTris) printf("danger nbTris %d %d\n", iCount, nbNodes);
      //if (i1 > nbNodes) printf("danger1 %d %d\n",i1,nbNodes);
      //if (i2 > nbNodes) printf("danger2 %d %d\n",i2,nbNodes);
      //if (i3 > nbNodes) printf("danger3 %d %d\n",i3,nbNodes);
      //printf("TRI %d: %d %d %d\n", iCount, i1,i2,i3);
    }
    cNodes += nodes->Upper()-nodes->Lower()+1;
    cTris += tri->NbTriangles();
  }
  
  PyList_Append(out, o); Py_DECREF(o);
  RELEASESHAREDU(o, f, c);

  // Pour exporter en stl il faut deja creer la triangulation dans la topo
  if (strcmp(outFileFormat, "fmt_stl") == 0)
  {
    // push binary
    StlAPI_Writer writer;
    writer.Write(shape, outFileName);
    printf("done writing %s.\n", outFileName);
  }
  
  // export en VRML
  if (strcmp(outFileFormat, "fmt_vrml") == 0)
  {
    VrmlAPI_Writer writer;
    writer.Write(shape, outFileName);
    printf("done writing %s.\n", outFileName);
  }
  
  return out;
}

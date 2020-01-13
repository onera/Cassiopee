/*    
    Copyright 2013-2020 Onera.

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

#include "STEPControl_Reader.hxx"
#include "IGESControl_Reader.hxx"
#include "StlAPI_Writer.hxx"
#include "VrmlAPI_Writer.hxx"
#include "IGESControl_Controller.hxx"
#include "STEPControl_Writer.hxx"

#include "TopoDS_Shape.hxx"
#include "Standard_ConstructionError.hxx"
#include "occ.h"
#include "kcore.h"
#include <stdio.h>

#include <iostream>
#include <iomanip>
#include "BRepPrimAPI_MakeCylinder.hxx"
#include "BRepPrimAPI_MakeBox.hxx"
#include "BRepAlgoAPI_Cut.hxx"
#include "BRepGProp.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "GProp_GProps.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "Poly_Triangulation.hxx"

// ============================================================================
/* Essai avec utilisation de BrepMesh */
// ============================================================================
PyObject* K_OCC::convertCAD2Arrays(PyObject* self, PyObject* args)
{ 
  // inFileFormat : fmt_iges ou fmt_step
  // outFileFormat : None ou fmt_stl ou fmt_vrml
  // deflection : deviation pour le maillage
  char* inFileName; char* inFileFormat; 
  char* outFileName; char* outFileFormat;
  E_Float deflection;
  if (!PYPARSETUPLE(args, "ssssd", "ssssd",
                    "ssssf", "ssssf", &inFileName, &inFileFormat, 
                    &outFileName, &outFileFormat,
                    &deflection)) return NULL;

  //Create a simple box with a size 100x100x50, centered around the origin
  /*
  gp_Pnt lowerLeftCornerOfBox(-50.0,-50.0,0.0);
  BRepPrimAPI_MakeBox boxMaker(lowerLeftCornerOfBox,100,100,50);
  TopoDS_Shape box = boxMaker.Shape();
  */
    
  //Create a cylinder with a radius 25.0 and height 50.0, centered at the origin 
  /*
  BRepPrimAPI_MakeCylinder cylinderMaker(25.0,50.0);
  TopoDS_Shape cylinder = cylinderMaker.Shape();
  */
    
  //Cut the cylinder out from the box
  /*
  BRepAlgoAPI_Cut cutMaker(box,cylinder);
  TopoDS_Shape boxWithHole = cutMaker.Shape();
  */
    
  //Write the resulting shape to a file
  /*
  STEPControl_Writer writer;
  writer.Transfer(boxWithHole, STEPControl_AsIs);
  writer.Write("boxWithHole.stp");
  */
    
  //We compute some volumetric properties of the resulting shape
  /*
  GProp_GProps volumeProperties;
  BRepGProp::VolumeProperties(boxWithHole,volumeProperties);
  */
    
  //Compute the volume of the model
  //std::cout << std::setprecision(14) << "Volume of the model is: " << volumeProperties.Mass() << std::endl;
  
  //Compute the center of mass
  //std::cout << "Center of mass is: " << volumeProperties.CentreOfMass().X() << " " << volumeProperties.CentreOfMass().Y() << " " << volumeProperties.CentreOfMass().Z() << std::endl;

  /*
  const Standard_Real aLinearDeflection=0.1;
  if (BrepTools::Triangulation(boxWithHole, aLinearDeflection) == false)
  {
    BrepMesh_IncrementalMesh aMesh(boxWithHole, 0.1);
    Standard_Boolean bDone = aMesh.isDone();
    if (aMesh.isModified()) boxWithHole = aMesh.Shape();
  }
  */
  
  TopoDS_Shape shape;
  if (strcmp(inFileFormat, "fmt_step") == 0)
  {
    STEPControl_Reader reader;
    reader.ReadFile(inFileName);
    reader.TransferRoots();
    shape = reader.OneShape();
    printf("done reading %s.\n",inFileName);
  }
  if (strcmp(inFileFormat, "fmt_iges") == 0)
  {
    IGESControl_Controller::Init();
    IGESControl_Reader reader;
    IFSelect_ReturnStatus stat = reader.ReadFile(inFileName);
    reader.SetReadVisible(Standard_True);
    reader.PrintCheckLoad(Standard_True, IFSelect_GeneralInfo);
    reader.ClearShapes();
    reader.TransferRoots();
    shape = reader.OneShape();
    printf("done reading %s.\n", inFileName); 
  }
  // Triangulate
  BRepMesh_IncrementalMesh Mesh(shape, deflection);
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

    const TColgp_Array1OfPnt &nodes = tri->Nodes();
    for (Standard_Integer iCount = nodes.Lower(); iCount <= nodes.Upper(); iCount++)
    {
      nbNodes++;  
    }

    const Poly_Array1OfTriangle &tris = tri->Triangles();
    for (Standard_Integer iCount = tris.Lower(); iCount <= tris.Upper(); iCount++)
    {
      nbTris++;
    }
  }
  
  printf("Nbre total de noeuds: %d\n", nbNodes);
  printf("Nbre total de tris: %d\n", nbTris);
  
  // buildArray
  PyObject* o = K_ARRAY::buildArray2(3, "x,y,z", nbNodes, nbTris, -1, "TRI", false, 0, 0, 0, 1);
  FldArrayF* f; FldArrayI* c;
  K_ARRAY::getFromArray2(o, f, c);
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

    const TColgp_Array1OfPnt &nodes = tri->Nodes();
    // get the nodes
    for (Standard_Integer iCount = nodes.Lower(); iCount <= nodes.Upper(); iCount++)
    {
      Standard_Real x,y,z;
      nodes(iCount).Coord(x,y,z);
      xloc.Transforms(x,y,z);
      fx[cNodes+iCount-1] = x;
      fy[cNodes+iCount-1] = y;
      fz[cNodes+iCount-1] = z;
      //if (iCount > nbNodes) printf("danger iCount>nbNodes %d %d\n", iCount, nbNodes);
      //printf("point %d: %f %f %f\n",iCount,x,y,z);
    }
    
    // copy the polygons
    Standard_Integer i1, i2, i3;
    const Poly_Array1OfTriangle &tris = tri->Triangles();
    for (Standard_Integer iCount = tris.Lower(); iCount <= tris.Upper(); iCount++)
    {
      // get the node indexes for this triangle
      Poly_Triangle tri = tris(iCount);
      tri.Get(i1, i2, i3);
      c1[stride*(cTris+iCount-1)] = i1+cNodes;
      c2[stride*(cTris+iCount-1)] = i2+cNodes;
      c3[stride*(cTris+iCount-1)] = i3+cNodes;
      //if (iCount > nbTris) printf("danger nbTris %d %d\n", iCount, nbNodes);
      //if (i1 > nbNodes) printf("danger1 %d %d\n",i1,nbNodes);
      //if (i2 > nbNodes) printf("danger2 %d %d\n",i2,nbNodes);
      //if (i3 > nbNodes) printf("danger3 %d %d\n",i3,nbNodes);
      //printf("TRI %d: %d %d %d\n", iCount, i1,i2,i3);
    }
    cNodes += nodes.Upper()-nodes.Lower()+1;
    cTris += tris.Upper()-tris.Lower()+1;
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

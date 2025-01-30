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
# include "generator.h"
using namespace std;
using namespace K_FLD;

#include <mystdlib.h>
#include <meshing.hpp>
#include <nginterface.h>
#include <stlgeom.hpp>
using namespace netgen;

namespace netgen
{
  MeshingParameters mparam;
  int id = 0, ntasks = 1;
  //void Ng_PrintDest(const char* s) { printf("%s.", s); }
  //void MyError(const char* ch) { printf("Error: %s.", ch); }
  void Ng_PrintDest(const char* s) { }
  void MyError(const char* ch) { }
  double GetTime () { return 0; }
  void Render() {}
}

namespace nglib {

// Data type for NETGEN mesh
typedef void* Ng_Mesh;
typedef void* Ng_Geometry_2D;
typedef void* Ng_STL_Geometry;

  Array<STLReadTriangle> readtrias; //only before initstlgeometry
  Array<Point<3> > readedges; //only before init stlgeometry

// *** Special Enum types used within Netgen ***********
// Currently implemented surface element types
enum Ng_Surface_Element_Type 
{ NG_TRIG=1, NG_QUAD=2, NG_TRIG6=3, NG_QUAD6=4, NG_QUAD8=5 };

// Currently implemented volume element types
enum Ng_Volume_Element_Type 
{ NG_TET=1, NG_PYRAMID=2, NG_PRISM=3, NG_TET10=4 };

// initialize, deconstruct Netgen library:
void Ng_Init()
{
  mycout = &cout;
  myerr = &cerr;
  //netgen::testout->SetOutStream (new ofstream ("test.out"));
  testout = new ofstream ("test.out");
}

// Clean-up functions before ending usage of nglib
void Ng_Exit()
{
  ;
}
 
// Create a new netgen mesh object
Ng_Mesh* Ng_NewMesh()
{
   Mesh* mesh = new Mesh;  
   mesh->AddFaceDescriptor (FaceDescriptor (1, 1, 0, 1));
   return (Ng_Mesh*)(void*)mesh;
}

// Delete an existing netgen mesh object
void Ng_DeleteMesh (Ng_Mesh* mesh)
{
  if (mesh != NULL)
  {
    // Delete the Mesh structures
    ((Mesh*)mesh)->DeleteMesh();
    
    // Now delete the Mesh class itself
    delete (Mesh*)mesh;
    
    // Set the Ng_Mesh pointer to NULL
    mesh = NULL;
  }
}

// Manually add a point to an existing mesh object
void Ng_AddPoint (Ng_Mesh* mesh, double* x)
{
  Mesh* m = (Mesh*)mesh;
  m->AddPoint (Point3d (x[0], x[1], x[2]));
}

// Manually add a surface element of a given type to an existing mesh object
void Ng_AddSurfaceElement (Ng_Mesh* mesh, Ng_Surface_Element_Type et,
                           int* pi)
{
  Mesh * m = (Mesh*)mesh;
  Element2d el (3);
  el.SetIndex (1);
  el.PNum(1) = pi[0];
  el.PNum(2) = pi[1];
  el.PNum(3) = pi[2];
  m->AddSurfaceElement (el);
}

// Manually add a volume element of a given type to an existing mesh object
void Ng_AddVolumeElement (Ng_Mesh* mesh, Ng_Volume_Element_Type et,
                          int* pi)
{
  Mesh * m = (Mesh*)mesh;
  Element el (4);
  el.SetIndex (1);
  el.PNum(1) = pi[0];
  el.PNum(2) = pi[1];
  el.PNum(3) = pi[2];
  el.PNum(4) = pi[3];
  m->AddVolumeElement (el);
}

// Obtain the number of points in the mesh
int Ng_GetNP (Ng_Mesh * mesh)
{
  return ((Mesh*)mesh) -> GetNP();
}

// Obtain the number of surface elements in the mesh
int Ng_GetNSE (Ng_Mesh * mesh)
{
  return ((Mesh*)mesh) -> GetNSE();
}

// Obtain the number of volume elements in the mesh
int Ng_GetNE (Ng_Mesh * mesh)
{
  return ((Mesh*)mesh) -> GetNE();
}

//  Return point coordinates of a given point index in the mesh
void Ng_GetPoint (Ng_Mesh * mesh, int num, double * x)
{
  const Point3d & p = ((Mesh*)mesh)->Point(num);
  x[0] = p.X();
  x[1] = p.Y();
  x[2] = p.Z();
}

// Return the surface element at a given index "pi"
Ng_Surface_Element_Type Ng_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi)
{
  const Element2d & el = ((Mesh*)mesh)->SurfaceElement(num);
  for (int i = 1; i <= el.GetNP(); i++)
    pi[i-1] = el.PNum(i);
  Ng_Surface_Element_Type et;
  switch (el.GetNP())
  {
    case 3: et = NG_TRIG; break;
    case 4: et = NG_QUAD; break;
    case 6: 
      switch (el.GetNV())
      {
        case 3: et = NG_TRIG6; break;
        case 4: et = NG_QUAD6; break;
        default:
          et = NG_TRIG6; break;
      }
      break;
    case 8: et = NG_QUAD8; break;
    default:
      et = NG_TRIG; break; // for the compiler
  }
  return et;
}

// Return the volume element at a given index "pi"
Ng_Volume_Element_Type
Ng_GetVolumeElement(Ng_Mesh* mesh, int num, int *pi)
{
  const Element & el = ((Mesh*)mesh)->VolumeElement(num);
  for (int i = 1; i <= el.GetNP(); i++)
    pi[i-1] = el.PNum(i);
  Ng_Volume_Element_Type et;
  switch (el.GetNP())
  {
    case 4: et = NG_TET; break;
    case 5: et = NG_PYRAMID; break;
    case 6: et = NG_PRISM; break;
    case 10: et = NG_TET10; break;
    default:
      et = NG_TET; break; // for the compiler
  }
  return et;
}

// Set a global limit on the maximum mesh size allowed
void Ng_RestrictMeshSizeGlobal(Ng_Mesh* mesh, double h)
{
  ((Mesh*)mesh) -> SetGlobalH(h);
}

// Set a local limit on the maximum mesh size allowed around the given point
void Ng_RestrictMeshSizePoint(Ng_Mesh * mesh, double* p, double h)
{
  ((Mesh*)mesh) -> RestrictLocalH (Point3d (p[0], p[1], p[2]), h);
}

// Set a local limit on the maximum mesh size allowed within a given box region
void Ng_RestrictMeshSizeBox(Ng_Mesh* mesh, double* pmin, double* pmax, double h)
{
  for (double x = pmin[0]; x < pmax[0]; x += h)
    for (double y = pmin[1]; y < pmax[1]; y += h)
      for (double z = pmin[2]; z < pmax[2]; z += h)
        ((Mesh*)mesh) -> RestrictLocalH (Point3d (x, y, z), h);
}

// Generates volume mesh from an existing surface mesh
E_Int Ng_GenerateVolumeMesh(Ng_Mesh* mesh)
{
  Mesh* m = (Mesh*)mesh;
  m->CalcLocalH(mparam.grading);

  E_Int ret = MeshVolume(mparam, *m);
  RemoveIllegalElements(*m);
  OptimizeVolume(mparam, *m);
  return ret;
}

void Ng_AddBoundarySeg_2D(Ng_Mesh * mesh, int pi1, int pi2)
{
  Mesh* m = (Mesh*)mesh;  
  Segment seg;
  seg[0] = pi1;
  seg[1] = pi2;
  m->AddSegment(seg);
}

int Ng_GetNP_2D(Ng_Mesh* mesh)
{
  Mesh * m = (Mesh*)mesh;
  return m->GetNP();
}

int Ng_GetNE_2D(Ng_Mesh* mesh)
{
  Mesh * m = (Mesh*)mesh;
  return m->GetNSE();
}

int Ng_GetNSeg_2D(Ng_Mesh* mesh)
{
  Mesh * m = (Mesh*)mesh;
  return m->GetNSeg();
}

Ng_Surface_Element_Type
Ng_GetElement_2D(Ng_Mesh* mesh, int num, int* pi, int* matnum)
{
  const Element2d & el = ((Mesh*)mesh)->SurfaceElement(num);
  for (int i = 1; i <= el.GetNP(); i++)
    pi[i-1] = el.PNum(i);

  Ng_Surface_Element_Type et;
  switch (el.GetNP())
  {
    case 3: et = NG_TRIG; break;
    case 4: et = NG_QUAD; break;
    case 6: 
      switch (el.GetNV())
      {
        case 3: et = NG_TRIG6; break;
        case 4: et = NG_QUAD6; break;
        default:
          et = NG_TRIG6; break;
      }
      break;
    case 8: et = NG_QUAD8; break;
    default:
      et = NG_TRIG; break; // for the compiler
  }
  
  if (matnum) *matnum = el.GetIndex();

  return et;
}

// generate new STL Geometry
Ng_STL_Geometry* Ng_STL_NewGeometry ()
{
  return (Ng_STL_Geometry*)(void*)new STLGeometry;
} 

// fills STL Geometry
// positive orientation
// normal vector may be null-pointer
void Ng_STL_AddTriangle (Ng_STL_Geometry * geom, 
                         double * p1, double * p2, double * p3, 
                         double * nv)
{
  Point<3> apts[3];
  apts[0] = Point<3>(p1[0],p1[1],p1[2]);
  apts[1] = Point<3>(p2[0],p2[1],p2[2]);
  apts[2] = Point<3>(p3[0],p3[1],p3[2]);
  
  Vec<3> n;
  if (!nv)
    n = Cross (apts[0]-apts[1], apts[0]-apts[2]);
  else
    n = Vec<3>(nv[0],nv[1],nv[2]);
  
  readtrias.Append(STLReadTriangle(apts,n));
}

// add (optional) edges:
void Ng_STL_AddEdge (Ng_STL_Geometry * geom, 
                     double * p1, double * p2)
{
  readedges.Append(Point3d(p1[0],p1[1],p1[2]));
  readedges.Append(Point3d(p2[0],p2[1],p2[2]));
}

// loads geometry from STL file
Ng_STL_Geometry * Ng_STL_LoadGeometry (const char * filename, int binary)
{
  int i;
  STLGeometry geom;
  STLGeometry* geo;
  ifstream ist(filename);

  if (binary) geo = geom.LoadBinary(ist);
  else geo = geom.Load(ist);

  readtrias.SetSize(0);
  readedges.SetSize(0);

  Point3d p;
  Vec3d normal;
  double p1[3]; double p2[3]; double p3[3];
  double n[3];
  
  Ng_STL_Geometry * geo2 = Ng_STL_NewGeometry();
  
  for (i = 1; i <= geo->GetNT(); i++)
  {
    const STLTriangle& t = geo->GetTriangle(i);
    p = geo->GetPoint(t.PNum(1));
    p1[0] = p.X(); p1[1] = p.Y(); p1[2] = p.Z(); 
    p = geo->GetPoint(t.PNum(2));
    p2[0] = p.X(); p2[1] = p.Y(); p2[2] = p.Z(); 
    p = geo->GetPoint(t.PNum(3));
    p3[0] = p.X(); p3[1] = p.Y(); p3[2] = p.Z();
    normal = t.Normal();
    n[0] = normal.X(); n[1] = normal.Y(); n[2] = normal.Z();
    
    Ng_STL_AddTriangle(geo2, p1, p2, p3, n);
  }
  
  return geo2;
}


// after adding triangles (and edges) initialize
E_Int Ng_STL_InitSTLGeometry (Ng_STL_Geometry * geom)
{
  STLGeometry* geo = (STLGeometry*)geom;
  geo->InitSTLGeometry(readtrias);
  readtrias.SetSize(0);
  
  if (readedges.Size() != 0)
  {
    /*
      for (int i = 1; i <= readedges.Size(); i+=2)
      {
      cout << "e(" << readedges.Get(i) << "," << readedges.Get(i+1) << ")" << endl;
      }
    */
    geo->AddEdges(readedges);
  }
  
  if (geo->GetStatus() == STLTopology::STL_GOOD || geo->GetStatus() == STLTopology::STL_WARNING) return 0;
  return 1;
}

// automatically generates edges:
E_Int Ng_STL_MakeEdges(Ng_STL_Geometry * geom,
                       Ng_Mesh* mesh)
{
  STLGeometry* stlgeometry = (STLGeometry*)geom;
  Mesh* me = (Mesh*)mesh;
  
  me -> SetGlobalH(mparam.maxh);
  me -> SetLocalH(stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
                  stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
                  0.3);
  
  me -> LoadLocalMeshSize (mparam.meshsizefilename);
  /*
    if (mparam.meshsizefilename)
    {
    ifstream infile (mparam.meshsizefilename);
    if (!infile.good()) return NG_FILE_NOT_FOUND;
    me -> LoadLocalMeshSize (infile);
    }
  */
  
  STLMeshing (*stlgeometry, *me);
  stlgeometry->edgesfound = 1;
  stlgeometry->surfacemeshed = 0;
  stlgeometry->surfaceoptimized = 0;
  stlgeometry->volumemeshed = 0;
  
  return 0;
}

// generates mesh, empty mesh be already created.
E_Int Ng_STL_GenerateSurfaceMesh (Ng_STL_Geometry* geom,
                                  Ng_Mesh* mesh)
{
  STLGeometry* stlgeometry = (STLGeometry*)geom;
  Mesh* me = (Mesh*)mesh;

  /*
    me -> SetGlobalH (mparam.maxh);
    me -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
    stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
    0.3);
  */
  /*
    STLMeshing (*stlgeometry, *me);
    
    stlgeometry->edgesfound = 1;
    stlgeometry->surfacemeshed = 0;
    stlgeometry->surfaceoptimized = 0;
    stlgeometry->volumemeshed = 0;
  */  
  int retval = STLSurfaceMeshing (*stlgeometry, *me);
  if (retval == MESHING3_OK)
  {
    (*mycout) << "Success !!!!" << endl;
    stlgeometry->surfacemeshed = 1;
    stlgeometry->surfaceoptimized = 0;
    stlgeometry->volumemeshed = 0;
  } 
  else if (retval == MESHING3_OUTERSTEPSEXCEEDED)
  {
    (*mycout) << "ERROR: Give up because of too many trials. Meshing aborted!" << endl;
  }
  else if (retval == MESHING3_TERMINATE)
  {
    (*mycout) << "Meshing Stopped!" << endl;
  }
  else
  {
    (*mycout) << "ERROR: Surface meshing not successful. Meshing aborted!" << endl;
  }

  STLSurfaceOptimization (*stlgeometry, *me, mparam);

  return 0;
}

} // namespace nglib

//=========================================================================
/* Generation de maillage tetra a partir d'un maillage surfacique (netgen) 
   Le maillage surfacique n'est pas modifie */
//=========================================================================
PyObject* K_GENERATOR::netgen1(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float grading, maxh;
#if defined E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "Odd", &array, &maxh, &grading)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Off", &array, &maxh, &grading)) return NULL;
#endif
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res <= 0) return NULL;
  if (res == 1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "netgen: input must be TRI.");
    return NULL;
  }

  if (K_STRING::cmp(eltType, "TRI") != 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "netgen: input must be TRI.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "netgen: coordinates not found in array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  
  // Reglages des parametres de netgen

  // Normal
  mparam.uselocalh = 1;
  mparam.maxh = maxh;
  mparam.grading = grading;

  Ng_Mesh* mesh; 
  double point[3]; int trig[3]; int tet[4];
  nglib::Ng_Init();
  mesh = nglib::Ng_NewMesh();

  E_Int np = f->getSize();
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  /* Ajoute les points */
  for (E_Int i = 0; i < np; i++)
  {
    point[0] = x[i]; point[1] = y[i]; point[2] = z[i];
    //printf("%f %f %f\n", point[0], point[1], point[2]);
    nglib::Ng_AddPoint(mesh, point);
  }

  /* Ajoute les connectivites */
  E_Int ne = cn->getSize();
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  for (E_Int i = 0; i < ne; i++) 
  {
    trig[0] = cn1[i]; trig[1] = cn2[i]; trig[2] = cn3[i];
    //printf("%d %d %d\n", trig[0], trig[1], trig[2]);
    nglib::Ng_AddSurfaceElement(mesh, nglib::NG_TRIG, trig);
  }

  // Maillage volumique
  printf("start meshing.\n");
  E_Int ret = nglib::Ng_GenerateVolumeMesh(mesh);
  if (ret != 0)
  {
    printf("meshing failed (%d).\n", ret);
    nglib::Ng_DeleteMesh(mesh);
    RELEASESHAREDU(array, f, cn);
    if (ret == 1)
    {
      PyErr_SetString(PyExc_ValueError,
                      "netgen: mesh has still open quads (FAILED).");
    }
    else if (ret == 3)
    {
      PyErr_SetString(PyExc_ValueError,
                      "netgen: algorithm fails (FAILED).");
    }
    else if (ret == 5)
    {
      PyErr_SetString(PyExc_ValueError,
                      "netgen: input surface is invalid, overlap?, hole? (FAILED).");
    }
    else
    {
      PyErr_SetString(PyExc_ValueError,
                      "netgen: mesh generation fails.");
    }
    return NULL;
  } else printf("meshing done.\n");

  PyObject* out = NULL;
  // TETRA mesh output
  np = nglib::Ng_GetNP(mesh);
  printf("Generate %d points.\n",np);
  ne = nglib::Ng_GetNE(mesh);
  printf("Generate %d elements.\n",ne);
  
  /* Build output array */
  out = K_ARRAY::buildArray(3, "x,y,z", np, ne, 4, NULL, false, 0);
  E_Float* fp = K_ARRAY::getFieldPtr(out);
  E_Float* fx = fp; E_Float* fy = fp+np; E_Float* fz = fp+2*np;
  E_Int* cp = K_ARRAY::getConnectPtr(out);
  E_Int* cp1 = cp; E_Int* cp2 = cp+ne; 
  E_Int* cp3 = cp+2*ne; E_Int* cp4 = cp+3*ne; 

  for (E_Int i = 0; i < np; i++)
  {
    nglib::Ng_GetPoint (mesh, i+1, point);
    fx[i] = point[0]; fy[i] = point[1]; fz[i] = point[2];
    //printf("%f %f %f\n", fx[i], fy[i], fz[i]);
  }
  
  for (E_Int i = 0; i < ne; i++)
  {
    nglib::Ng_GetVolumeElement (mesh, i+1, tet);
    cp1[i] = tet[0]; cp2[i] = tet[1];
    cp3[i] = tet[2]; cp4[i] = tet[3];
    //printf("%d %d %d %d\n", cp1[i], cp2[i], cp3[i], cp4[i]);
  }
 
  nglib::Ng_DeleteMesh(mesh);

  RELEASESHAREDU(array, f, cn);

  return out;
}

//=========================================================================
/* Generation de maillage tetra a partir d'un maillage surfacique (netgen) 
   avec remaillage des surfaces */
//=========================================================================
PyObject* K_GENERATOR::netgen2(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float grading, maxh;
#if defined E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "Odd", &array, &maxh, &grading)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Off", &array, &maxh, &grading)) return NULL;
#endif
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res <= 0) return NULL;
  if (res == 1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "netgen: input must be TRI.");
    return NULL;
  }

  if (K_STRING::cmp(eltType, "TRI") != 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "netgen: input must be TRI.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "netgen: coordinates not found in array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  
  // Reglages des parametres de netgen
  mparam.uselocalh = 1;
  mparam.maxh = maxh;
  mparam.grading = grading;
  mparam.secondorder = 0;

  nglib::Ng_Mesh* mesh; 
  nglib::Ng_STL_Geometry* stl_geom;
  E_Int ret;

  nglib::Ng_Init();
  mesh = nglib::Ng_NewMesh();

  // Dump and read for now (in memory will be better!)
  stl_geom = nglib::Ng_STL_LoadGeometry("netgen.0120.stl", 0);

  ret = nglib::Ng_STL_InitSTLGeometry(stl_geom);
  if (ret != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "netgen: fails in surface analysis.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  printf("Start Edge Meshing....\n");
  ret = nglib::Ng_STL_MakeEdges(stl_geom, mesh);
  if (ret != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "netgen: fails in edge remeshing.");
    RELEASESHAREDU(array, f, cn); return NULL;
   }

   printf("Start Surface Meshing....\n");
   ret = nglib::Ng_STL_GenerateSurfaceMesh(stl_geom, mesh);
   if (ret != 0)
   {
     PyErr_SetString(PyExc_TypeError,
                     "netgen: fails in surface remeshing.");
     RELEASESHAREDU(array, f, cn); return NULL;
   }
   
   printf("Start Volume Meshing....\n");
   ret = nglib::Ng_GenerateVolumeMesh (mesh);
   if (ret != 0)
   {
     PyErr_SetString(PyExc_TypeError,
                     "netgen: fails in volume meshing.");
     RELEASESHAREDU(array, f, cn); return NULL;
   }
   
   printf("Meshing successfully completed....\n");

   E_Int np, ne;
   double point[3]; int tet[4];

   PyObject* out = NULL;
   // TETRA mesh output
   np = nglib::Ng_GetNP(mesh);
   printf("Generate %d points.\n",np);
   ne = nglib::Ng_GetNE(mesh);
   printf("Generate %d elements.\n",ne);
  
  /* Build output array */
  out = K_ARRAY::buildArray(3, "x,y,z", np, ne, 4, NULL, false, 0);
  E_Float* fp = K_ARRAY::getFieldPtr(out);
  E_Float* fx = fp; E_Float* fy = fp+np; E_Float* fz = fp+2*np;
  E_Int* cp = K_ARRAY::getConnectPtr(out);
  E_Int* cp1 = cp; E_Int* cp2 = cp+ne; 
  E_Int* cp3 = cp+2*ne; E_Int* cp4 = cp+3*ne; 

  for (E_Int i = 0; i < np; i++)
  {
    nglib::Ng_GetPoint (mesh, i+1, point);
    fx[i] = point[0]; fy[i] = point[1]; fz[i] = point[2];
    //printf("%f %f %f\n", fx[i], fy[i], fz[i]);
  }
  
  for (E_Int i = 0; i < ne; i++)
  {
    nglib::Ng_GetVolumeElement (mesh, i+1, tet);
    cp1[i] = tet[0]; cp2[i] = tet[1];
    cp3[i] = tet[2]; cp4[i] = tet[3];
    //printf("%d %d %d %d\n", cp1[i], cp2[i], cp3[i], cp4[i]);
  }
 
  nglib::Ng_DeleteMesh(mesh);

  RELEASESHAREDU(array, f, cn);

  return out;
}

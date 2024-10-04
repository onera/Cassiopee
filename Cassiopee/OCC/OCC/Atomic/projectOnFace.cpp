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

#include "occ.h"
#include "TopoDS_Face.hxx"
#include "GeomAPI_ProjectPointOnSurf.hxx"
#include "BRep_Tool.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "StdFail_NotDone.hxx"

#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepExtrema_DistShapeShape.hxx"

#include "Nuga/include/BbTree.h"
#include "Bnd_Box.hxx"
#include "BRepBndLib.hxx"

// Project coords on CAD face
void projectOnFace__(E_Int npts, E_Float* px, E_Float* py, E_Float* pz, const TopoDS_Face& F)
{
  Handle(Geom_Surface) face = BRep_Tool::Surface(F);

#pragma omp parallel
  {
    gp_Pnt Point;

#pragma omp for
    for (E_Int i=0; i < npts; i++)
    {
      Point.SetCoord(px[i], py[i], pz[i]);
      try
      {  
        GeomAPI_ProjectPointOnSurf o(Point, face, Extrema_ExtAlgo_Tree);
        gp_Pnt Pj = o.NearestPoint();
        printf("Info: projection %f %f %f -> %f %f %f\n",px[i],py[i],pz[i],Pj.X(),Pj.Y(),Pj.Z());
        px[i] = Pj.X(); py[i] = Pj.Y(); pz[i] = Pj.Z();
      }
      catch (StdFail_NotDone& e) { printf("Info: projectOnFace: fail for point %g %g %g.\n", px[i],py[i],pz[i]); }
    }
  }
}

// ============================================================================
/* Project array in place 
  IN: hook: CAD tree hook
  IN: array: array to project
  IN: faceList: list of no of faces (starting 1)
*/
// ============================================================================
// 0: proj ortho, 1: par EXTREMA
#define PROJMETHOD 1
PyObject* K_OCC::projectOnFaces(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* array; PyObject* faceList;
  if (!PYPARSETUPLE_(args, OOO_, &hook, &array, &faceList)) return NULL;  

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  // array a projeter
  FldArrayF* fi; E_Int ni, nj, nk;
  char* varString; FldArrayI* c; char* eltType;
  E_Int ret = K_ARRAY::getFromArray2(array, varString, fi, ni, nj, nk, c, eltType);
  if (ret != 1 && ret != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectOnFaces: invalid array.");
    return NULL;
  }

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  //TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];

  // liste des no des faces sur lesquelles on projete
  FldArrayI faces;
  if (faceList == Py_None)
  { 
    E_Int nfaces = surfaces.Extent(); 
    faces.malloc(nfaces);
    for (E_Int i = 0; i < nfaces; i++) faces[i] = i+1;
  }
  else K_ARRAY::getFromList(faceList, faces);

  E_Float* px = fi->begin(1); // fix
  E_Float* py = fi->begin(2);
  E_Float* pz = fi->begin(3);
  E_Int npts = fi->getSize();

  TopExp_Explorer expl;
  E_Int nfaces = faces.getSize();

  E_Float* ptx = new E_Float [npts];
  E_Float* pty = new E_Float [npts];
  E_Float* ptz = new E_Float [npts];
  E_Float* pox = new E_Float [npts];
  E_Float* poy = new E_Float [npts];
  E_Float* poz = new E_Float [npts];
  E_Float* dist = new E_Float [npts];

  //printf("nfaces=%d, npts=%d\n", nfaces, npts); fflush(stdout);
  //npts = 500;
  

#pragma omp parallel
  {
    gp_Pnt Point;
    gp_Pnt Pj;
    E_Float dx,dy,dz,d;

#pragma omp for
    for (E_Int i = 0; i < npts; i++) pox[i] = px[i];
#pragma omp for
    for (E_Int i = 0; i < npts; i++) poy[i] = py[i];
#pragma omp for
    for (E_Int i = 0; i < npts; i++) poz[i] = pz[i];
#pragma omp for
    for (E_Int i = 0; i < npts; i++) dist[i] = K_CONST::E_MAX_FLOAT;

    for (E_Int j = 0; j < nfaces; j++)
    {
        const TopoDS_Face& F = TopoDS::Face(surfaces(faces[j]));
        Handle(Geom_Surface) face = BRep_Tool::Surface(F);
        
#if PROJMETHOD == 0
        GeomAPI_ProjectPointOnSurf o;
#else
        BRepExtrema_DistShapeShape tool;
#endif

#pragma omp for
      for (E_Int i=0; i < npts; i++)
      {      
        Point.SetCoord(pox[i], poy[i], poz[i]);

#if PROJMETHOD == 1    
        BRepBuilderAPI_MakeVertex v(Point);
        const TopoDS_Vertex& Vertex = v.Vertex();
#endif
        try
        {  
#if PROJMETHOD == 0
          //GeomAPI_ProjectPointOnSurf o(Point, face, Extrema_ExtAlgo_Tree);
          GeomAPI_ProjectPointOnSurf o(Point, face, Extrema_ExtAlgo_Grad);
          
          //o.Init(Point, face, Extrema_ExtAlgo_Tree);
          o.Init(Point, face, Extrema_ExtAlgo_Grad); 
          Pj = o.NearestPoint();
          
          Pj = o.Point(1);
          //printf("projection %f %f %f -> %f %f %f\n",px[i],py[i],pz[i],Pj.X(),Pj.Y(),Pj.Z());
          ptx[i] = Pj.X(); pty[i] = Pj.Y(); ptz[i] = Pj.Z();
#else
          tool.LoadS1(F); tool.LoadS2(Vertex);
          tool.Perform();
          //E_Float dist = tool.Value();
          const gp_Pnt& Po = tool.PointOnShape1(1);
          ptx[i] = Po.X(); pty[i] = Po.Y(); ptz[i] = Po.Z();
#endif
        }
        catch (StdFail_NotDone& e)
        { 
          //printf("FAIL for point %g %g %g\n", px[i],py[i],pz[i]); 
          ptx[i] = K_CONST::E_MAX_FLOAT;
          pty[i] = K_CONST::E_MAX_FLOAT;
          ptz[i] = K_CONST::E_MAX_FLOAT;
        }

        dx = ptx[i]-pox[i];
        dy = pty[i]-poy[i];
        dz = ptz[i]-poz[i];
        d = dx*dx+dy*dy+dz*dz;
        if (d < dist[i])
        { dist[i] = d; px[i] = ptx[i]; py[i] = pty[i]; pz[i] = ptz[i]; }
        //printf("projection %f %f %f -> %f %f %f\n",pox[i],poy[i],poz[i],px[i],py[i],pz[i]);
      }
    }
  }

  delete [] pox; delete [] poy; delete [] poz;
  delete [] ptx; delete [] pty; delete [] ptz;
  delete [] dist;
  RELEASESHAREDB(ret, array, fi, c);
  Py_DECREF(Py_None);
  return Py_None;
}

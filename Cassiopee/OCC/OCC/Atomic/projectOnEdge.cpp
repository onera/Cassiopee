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
#include <TopoDS_Face.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <BRep_Tool.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <StdFail_NotDone.hxx>
#include <BRepAdaptor_Curve.hxx>

// ============================================================================
/* Project array in place 
  IN: hook: CAD tree hook
  IN: array: array to project
  IN: EdgeList: list of no of edges (starting 1)
*/
// ============================================================================
PyObject* K_OCC::projectOnEdges(PyObject* self, PyObject* args)
{
  PyObject* hook; PyObject* array; PyObject* edgeList;
  if (!PYPARSETUPLE_(args, OOO_, &hook, &array, &edgeList)) return NULL;  

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  // array a projeter
  FldArrayI* c; FldArrayF* fi;
  E_Int ni, nj, nk;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray2(array, varString, fi, ni, nj, nk, c, eltType);
  if (ret != 1 && ret != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectOnEdges: invalid array.");
    return NULL;
  }

  //TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopTools_IndexedMapOfShape& cadEdges = *(TopTools_IndexedMapOfShape*)packet[2];

  // liste des no des edges sur lesquelles on projete
  FldArrayI edges;
  if (edgeList == Py_None)
  { 
    E_Int nedges = cadEdges.Extent(); 
    edges.malloc(nedges);
    for (E_Int i = 0; i < nedges; i++) edges[i] = i+1;
  }
  else K_ARRAY::getFromList(edgeList, edges);

  E_Float* px = fi->begin(1); // fix
  E_Float* py = fi->begin(2);
  E_Float* pz = fi->begin(3);
  E_Int npts = fi->getSize();

  TopExp_Explorer expl;
  E_Int nedges = edges.getSize();

  E_Float* ptx = new E_Float [npts];
  E_Float* pty = new E_Float [npts];
  E_Float* ptz = new E_Float [npts];
  E_Float* pox = new E_Float [npts];
  E_Float* poy = new E_Float [npts];
  E_Float* poz = new E_Float [npts];
  E_Float* dist = new E_Float [npts];

#pragma omp parallel
  {
    gp_Pnt Point;
    E_Float dx,dy,dz,d;

#pragma omp for
    for (E_Int i = 0; i < npts; i++) pox[i] = px[i];
#pragma omp for
    for (E_Int i = 0; i < npts; i++) poy[i] = py[i];
#pragma omp for
    for (E_Int i = 0; i < npts; i++) poz[i] = pz[i];
#pragma omp for
    for (E_Int i = 0; i < npts; i++) dist[i] = K_CONST::E_MAX_FLOAT;

#pragma omp for
    for (E_Int i=0; i < npts; i++)
    {
      for (E_Int j = 0; j < nedges; j++)
      {
        const TopoDS_Edge& E = TopoDS::Edge(cadEdges(edges[j]));
        BRepAdaptor_Curve C0(E);
        Standard_Real aFirst=C0.FirstParameter(), aEnd=C0.LastParameter();
        Handle(Geom_Curve) aCurve = BRep_Tool::Curve(E, aFirst, aEnd);

        Point.SetCoord(pox[i], poy[i], poz[i]);
        try
        {
          GeomAPI_ProjectPointOnCurve o(Point, aCurve);
          gp_Pnt Pj = o.NearestPoint();
          //printf("projection %f %f %f -> %f %f %f\n",px[i],py[i],pz[i],Pj.X(),Pj.Y(),Pj.Z());
          ptx[i] = Pj.X(); pty[i] = Pj.Y(); ptz[i] = Pj.Z();
        }
        catch (StdFail_NotDone& e) 
        { 
          //printf("FAIL for point %g %g %g\n", px[i],py[i],pz[i]); 
          ptx[i] = K_CONST::E_MAX_FLOAT;
          pty[i] = K_CONST::E_MAX_FLOAT;
          ptz[i] = K_CONST::E_MAX_FLOAT;
        }
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

  delete [] pox; delete [] poy; delete [] poz;
  delete [] ptx; delete [] pty; delete [] ptz;
  delete [] dist;
  RELEASESHAREDB(ret, array, fi, c);
  Py_DECREF(Py_None);
  return Py_None;
}

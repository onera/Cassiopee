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
#include "TopoDS_Edge.hxx"
#include "Geom_Surface.hxx"
#include "BRep_Tool.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "BRepAdaptor_Curve.hxx"

// evalue une courbe u en x,y,z a partir de l'edge 
void evalEdge__(E_Int npts, E_Float* u, const TopoDS_Edge& E,
                E_Float* x, E_Float* y, E_Float* z)
{  
#pragma omp parallel
{
  BRepAdaptor_Curve C0(E);
  gp_Pnt Pt;
#pragma omp for
  for (E_Int i = 0; i < npts; i++)
  {
    C0.D0(u[i], Pt);
    x[i] = Pt.X(); y[i] = Pt.Y(); z[i] = Pt.Z();
  }
}
}

// evalEdge
// IN: arrayU: arrayU corresponding to the edge edgeNo of CAD
PyObject* K_OCC::evalEdge(PyObject* self, PyObject* args)
{
  PyObject* hook;
  PyObject* arrayU;
  E_Int edgeNo; // No de l'edge
  if (!PYPARSETUPLE_(args, OO_ I_, &hook, &arrayU, &edgeNo)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& edges = *(TopTools_IndexedMapOfShape*)packet[2];
  TopExp_Explorer expl;

  const TopoDS_Edge& E = TopoDS::Edge(edges(edgeNo));
  FldArrayF* fi; E_Int ni, nj, nk;
  char* varString; FldArrayI* ci; char* eltType;
  E_Int ret = K_ARRAY::getFromArray3(arrayU, varString, fi, ni, nj, nk, ci, eltType);
  E_Float* pu = fi->begin(1);
  PyObject* o;
  if (ret == 1) o = K_ARRAY::buildArray3(3, "x,y,z", ni, nj, nk, 1);
  else o = K_ARRAY::buildArray3(3, "x,y,z", fi->getSize(), ci->getSize(), eltType);
  FldArrayF* fo; FldArrayI* co;
  if (ret == 1) K_ARRAY::getFromArray3(o, fo);
  else
  {
    K_ARRAY::getFromArray3(o, fo, co);
    E_Int* pci = ci->begin(); E_Int* pco = co->begin();
    for (E_Int i = 0; i < ci->getSize()*ci->getNfld(); i++) pco[i]  = pci[i];
  }
  E_Float* px = fo->begin(1);
  E_Float* py = fo->begin(2);
  E_Float* pz = fo->begin(3);
  evalEdge__(fo->getSize(), pu, E, px, py, pz);
  RELEASESHAREDB(ret, arrayU, fi, ci);
  if (ret == 1) { RELEASESHAREDS(o, fo); }
  else { RELEASESHAREDU(o, fo, co); }
  return o;
}
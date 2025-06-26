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

#include "occ.h"
#include <gp_Pnt.hxx>
#include "TopoDS_Face.hxx"
#include "Geom_Surface.hxx"
#include "BRep_Tool.hxx"
#include "TopExp_Explorer.hxx"
#include "TopoDS.hxx"
#include "TopTools_IndexedMapOfShape.hxx"

// evalue une maillage uv en x,y,z a partir de la face 
void evalFace__(E_Int npts, E_Float* u, E_Float* v, const TopoDS_Face& F,
                E_Float* x, E_Float* y, E_Float* z)
{
#pragma omp parallel
  {
    Handle(Geom_Surface) face = BRep_Tool::Surface(F);
    gp_Pnt Pt;
    #pragma omp for 
    for (E_Int i = 0; i < npts; i++)
    {
      face->D0(u[i], v[i], Pt);
      x[i] = Pt.X(); y[i] = Pt.Y(); z[i] = Pt.Z();
    }
  }
}

// evalFace
// IN: arrayUV: arrayUV corresponding to the face faceNo of CAD
PyObject* K_OCC::evalFace(PyObject* self, PyObject* args)
{
  PyObject* hook;
  PyObject* arrayUV;
  E_Int faceNo; // No de la face 
  if (!PYPARSETUPLE_(args, OO_ I_, &hook, &arrayUV, &faceNo)) return NULL;

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopExp_Explorer expl;

  const TopoDS_Face& F = TopoDS::Face(surfaces(faceNo));
  FldArrayF* fi; E_Int ni, nj, nk;
  char* varString; FldArrayI* ci; char* eltType;
  E_Int ret = K_ARRAY::getFromArray2(arrayUV, varString, fi, ni, nj, nk, ci, eltType);
  E_Float* pu = fi->begin(1);
  E_Float* pv = fi->begin(2);
  PyObject* o = NULL;
  E_Int nfld = 3;
  if (nfld == 3)
  {
    if (ret == 1) o = K_ARRAY::buildArray2(nfld, "x,y,z", ni, nj, nk, 1);
    else o = K_ARRAY::buildArray2(nfld, "x,y,z,u,v", fi->getSize(), ci->getSize(), -1, eltType);
  }
  else if (nfld == 5)
  {
    if (ret == 1) o = K_ARRAY::buildArray2(nfld, "x,y,z", ni, nj, nk, 1);
    else o = K_ARRAY::buildArray2(nfld, "x,y,z,u,v", fi->getSize(), ci->getSize(), -1, eltType);
  }
  FldArrayF* fo; FldArrayI* co;
  if (ret == 1) K_ARRAY::getFromArray2(o, fo);
  else 
  {
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int* pci = ci->begin(); E_Int* pco = co->begin();
    for (E_Int i = 0; i < ci->getSize()*ci->getNfld(); i++) pco[i]  = pci[i];
  }
  E_Float* px = fo->begin(1);
  E_Float* py = fo->begin(2);
  E_Float* pz = fo->begin(3);
  evalFace__(fo->getSize(), pu, pv, F, px, py, pz);
  if (nfld == 5)
  {
    E_Float* pu2 = fo->begin(4);
    E_Float* pv2 = fo->begin(5);
    for (E_Int i = 0; i < fo->getSize(); i++) { pu2[i] = pu[i]; pv2[i] = pv[i]; } // keep uv
  }

  RELEASESHAREDB(ret, arrayUV, fi, ci);
  if (ret == 1) { RELEASESHAREDS(o, fo); }
  else { RELEASESHAREDU(o, fo, co); }
  return o;
}
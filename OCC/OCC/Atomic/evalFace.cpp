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

#include "occ.h"
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
    Handle(Geom_Surface) face = BRep_Tool::Surface(F);
#pragma omp parallel
{
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
// IN: faceUVArrays: list of arrays correspondinf to the face of CAD
PyObject* K_OCC::evalFace(PyObject* self, PyObject* args)
{
  PyObject* hook;
  PyObject* faceUVArrays;
  if (!PYPARSETUPLEF(args, "O", "O", &hook, &faceUVArrays)) return NULL;  

  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif

  TopTools_IndexedMapOfShape& surfaces = *(TopTools_IndexedMapOfShape*)packet[1];
  TopExp_Explorer expl;
  PyObject* out = PyList_New(0);
  for (E_Int i=1; i <= surfaces.Extent(); i++)
  {
    const TopoDS_Face& F = TopoDS::Face(surfaces(i));
    PyObject* array = PyList_GetItem(faceUVArrays, i-1);
    FldArrayF* fi; E_Int ni, nj, nk;
    char* varString; FldArrayI* c; char* eltType;
    E_Int ret = K_ARRAY::getFromArray2(array, varString, fi, ni, nj, nk, c, eltType);
    E_Float* pu = fi->begin(1);
    E_Float* pv = fi->begin(2);
    PyObject* o = K_ARRAY::buildArray2(3, "x,y,z", ni, nj, nk, 1);
    FldArrayF* fo; K_ARRAY::getFromArray2(o, fo);
    E_Float* px = fo->begin(1);
    E_Float* py = fo->begin(2);
    E_Float* pz = fo->begin(3);
    evalFace__(ni*nj, pu, pv, F, px, py, pz);
    RELEASESHAREDB(ret, array, fi, c);
    PyList_Append(out, o); Py_DECREF(o);
  }

  return out;

}
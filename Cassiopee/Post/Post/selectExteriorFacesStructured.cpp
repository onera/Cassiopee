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

// selectExteriorFacesStructured

# include <stdio.h>
# include <string.h>

# include "post.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Selectionne les facettes exterieures d'un array structure et les retourne 
   sous forme d'une liste d' arrays structures. */
// ============================================================================
PyObject* K_POST::selectExteriorFacesStructured(PyObject* self, PyObject* args)
{
  /* Check */ 
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorFacesStructured: invalid array.");
    return NULL;
  }
  if (res != 1) 
  {
    if (res == 2) { RELEASESHAREDU(array, f, cn); }
    PyErr_SetString(PyExc_TypeError,
                    "exteriorFacesStructured: array must be structured.");
    return NULL;
  }

  /* Exterior faces */
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();

  PyObject* l = PyList_New(0);
  PyObject* tpl;

  // 1D array
  if ((ni == 1 && nj == 1) || (nj == 1 && nk == 1) || (ni == 1 && nk == 1))
  {
    tpl = K_ARRAY::buildArray(nfld, varString, 1, 1, 1);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    for (E_Int eq = 0; eq < nfld; eq++) fnp[eq] = (*f)(0, eq+1);
    PyList_Append(l, tpl); Py_DECREF(tpl);

    tpl = K_ARRAY::buildArray(nfld, varString, 1, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    for (E_Int eq = 0; eq < nfld; eq++) fnp[eq] = (*f)(npts-1, eq+1);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  else if (ni == 1 || nj == 1 || nk == 1) // 2d
  {
    selectExteriorFacesStruct2D(ni, nj, nk, *f, varString, l);
  }
  else // 3d
  {
    selectExteriorFacesStruct3D(ni, nj, nk, *f, varString, l);
  }
  RELEASESHAREDS(array, f);
  return l;
}
//=============================================================================
void K_POST::selectExteriorFacesStruct3D(E_Int ni, E_Int nj, E_Int nk,
                                         FldArrayF& f, char* varString,
                                         PyObject* l)
{
  E_Int nfld = f.getNfld();
  E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
  E_Int ninj = ni*nj;
  E_Int shift;
  E_Float* fnp; PyObject* tpl;

  // i = 1
  tpl = K_ARRAY::buildArray(nfld, varString, nj, nk, 1);
  fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f2(nj*nk, nfld, fnp, true);

#pragma omp parallel default(shared)
  {
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f2.begin(nv);
      for (E_Int k = 0; k < nk; k++)
#pragma omp for
        for (E_Int j = 0; j < nj; j++)
          fv2[j+k*nj] = fv[j*ni+k*ninj];
    }
  }
  PyList_Append(l, tpl); Py_DECREF(tpl);

  // i = imax
  tpl = K_ARRAY::buildArray(nfld, varString, nj, nk, 1);
  fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f3(nj*nk, nfld, fnp, true);
  shift = ni1;
#pragma omp parallel default(shared)
  {
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f3.begin(nv);
      for (E_Int k = 0; k < nk; k++)
#pragma omp for
        for (E_Int j = 0; j < nj; j++)
          fv2[j+k*nj] = fv[shift+j*ni+k*ninj];
    }
  }
  PyList_Append(l, tpl); Py_DECREF(tpl);

  // j = 1
  tpl = K_ARRAY::buildArray(nfld, varString, ni, nk, 1);
  fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f4(ni*nk, nfld, fnp, true);
#pragma omp parallel default(shared)
  {
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f4.begin(nv);
      for (E_Int k = 0; k < nk; k++)
#pragma omp for
        for (E_Int i = 0; i < ni; i++)
          fv2[i+k*ni] = fv[i+k*ninj];
    }
  }
  PyList_Append(l, tpl); Py_DECREF(tpl);

  // j = jmax
  tpl = K_ARRAY::buildArray(nfld, varString, ni, nk, 1);
  fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f5(ni*nk, nfld, fnp, true);
  shift = nj1*ni;
#pragma omp parallel default(shared)
  {
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f5.begin(nv);
      for (E_Int k = 0; k < nk; k++)
#pragma omp for
        for (E_Int i = 0; i < ni; i++)
          fv2[i+k*ni] = fv[i+shift+k*ninj];
    }
  }
  PyList_Append(l, tpl); Py_DECREF(tpl);

  // k = 1
  tpl = K_ARRAY::buildArray(nfld, varString, ni, nj, 1);
  fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f6(ni*nj, nfld, fnp, true);
#pragma omp parallel default(shared)
  {
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f6.begin(nv);
      for (E_Int j = 0; j < nj; j++)
#pragma omp for
        for (E_Int i = 0; i < ni; i++)
          fv2[i+j*ni] = fv[i+j*ni];
    }
  }
  PyList_Append(l, tpl); Py_DECREF(tpl);

  // k = kmax
  tpl = K_ARRAY::buildArray(nfld, varString, ni, nj, 1);
  fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f7(ni*nj, nfld, fnp, true);
  shift = nk1*ninj;
#pragma omp parallel default(shared)
  {
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f7.begin(nv);
      for (E_Int j = 0; j < nj; j++)
#pragma omp for
        for (E_Int i = 0; i < ni; i++)
          fv2[i+j*ni] = fv[i+j*ni+shift];
    }
  }
  PyList_Append(l, tpl); Py_DECREF(tpl);
}
//=============================================================================
void K_POST::selectExteriorFacesStruct2D(E_Int ni, E_Int nj, E_Int nk,
                                         FldArrayF& f, char* varString, 
                                         PyObject* l)
{
  E_Int nfld = f.getNfld();
  E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
  E_Int ninj = ni*nj;
  E_Int shift;
  E_Float* fnp; PyObject* tpl;

  if (nk == 1)
  {
    // i = 1
    tpl = K_ARRAY::buildArray(nfld, varString, nj, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f2(nj, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f2.begin(nv);
      for (E_Int j = 0; j < nj; j++) fv2[j] = fv[j*ni];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
    
    // i = imax
    tpl = K_ARRAY::buildArray(nfld, varString, nj, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f3(nj, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f3.begin(nv);
      for (E_Int j = 0; j < nj; j++) fv2[j] = fv[ni1+j*ni];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
    
    // j = 1
    tpl = K_ARRAY::buildArray(nfld, varString, ni, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f4(ni, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f4.begin(nv);
      for (E_Int i = 0; i < ni; i++) fv2[i] = fv[i];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
    
    // j = jmax
    tpl = K_ARRAY::buildArray(nfld, varString, ni, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f5(ni, nfld, fnp, true);
    shift = nj1*ni;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f5.begin(nv);
      for (E_Int i = 0; i < ni; i++) fv2[i] = fv[i+shift];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  else if (nj == 1)
  {
    // i = 1
    tpl = K_ARRAY::buildArray(nfld, varString, nk, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f2(nk, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f2.begin(nv);
      for (E_Int k = 0; k < nk; k++) fv2[k] = fv[k*ninj];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
    
    // i = imax
    tpl = K_ARRAY::buildArray(nfld, varString, nk, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f3(nk, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f3.begin(nv);
      for (E_Int k = 0; k < nk; k++) fv2[k] = fv[ni1+k*ninj];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
    
    // k = 1
    tpl = K_ARRAY::buildArray(nfld, varString, ni, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f4(ni, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f4.begin(nv);
      for (E_Int i = 0; i < ni; i++) fv2[i] = fv[i];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
    
    // k = kmax
    tpl = K_ARRAY::buildArray(nfld, varString, ni, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f5(ni, nfld, fnp, true);
    shift = nk1*ninj;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f5.begin(nv);
      for (E_Int i = 0; i < ni; i++) fv2[i] = fv[i+shift];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  else //ni = 1
  {
    // j = 1
    tpl = K_ARRAY::buildArray(nfld, varString, nk, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f2(nk, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f2.begin(nv);
      for (E_Int k = 0; k < nk; k++) fv2[k] = fv[k*ninj];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
    
    // j= jmax
    tpl = K_ARRAY::buildArray(nfld, varString, nk, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f3(nk, nfld, fnp, true);
    shift = nj1*ni;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f3.begin(nv);
      for (E_Int k = 0; k < nk; k++) fv2[k] = fv[k*ninj+shift];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);

    // k = 1
    tpl = K_ARRAY::buildArray(nfld, varString, nj, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f4(nj, nfld, fnp, true);
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f4.begin(nv);
      for (E_Int j = 0; j < nj; j++) fv2[j] = fv[j*ni];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);

    // k = kmax
    tpl = K_ARRAY::buildArray(nfld, varString, nj, 1, 1);
    fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f5(nj, nfld, fnp, true);
    shift = nk1*ninj;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fv2 = f5.begin(nv);
      for (E_Int j = 0; j < nj; j++) fv2[j] = fv[j*ni+shift];
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
}
//=============== Post/selectExteriorFacesStructured.cpp ==================

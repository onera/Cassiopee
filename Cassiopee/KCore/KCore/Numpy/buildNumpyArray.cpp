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
# include "Numpy/Numpy.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Build a numpy array from a FldArray (Float)
   IN: field: Fld field
   IN: fortran=0 (C ordering), fortran=1 (Fortran ordering)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_NUMPY::buildNumpyArray(FldArrayF& field, E_Int fortran)
{
  npy_intp dim[2];
  PyArrayObject* a;
  IMPORTNUMPY;

  if (field.getNfld() == 1)
  {
    dim[0] = field.getSize();
    a = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, fortran);
    dim[1] = 1;
  }
  else
  {
    if (fortran == 0)
    {
      dim[1] = field.getSize(); dim[0] = field.getNfld();
    }
    else
    {
      dim[0] = field.getSize(); dim[1] = field.getNfld();
    }
    a = (PyArrayObject*)PyArray_EMPTY(2, dim, NPY_DOUBLE, fortran);
  }

  E_Float* ptrf = field.begin();
  E_Float* ptra = (E_Float*)PyArray_DATA(a);
  E_Int s = dim[0]*dim[1];
#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < s; i++) ptra[i] = ptrf[i];
  return (PyObject*)a;
}

//=============================================================================
/* Build a numpy array from a E_Float* + size
   IN: f: pointer (size, nfld)
   IN: size, nfld: sizes of f
   IN: fortran=0 (C ordering), fortran=1 (Fortran ordering)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_NUMPY::buildNumpyArray(E_Float* f, E_Int size, E_Int nfld,
                                   E_Int fortran)
{
  npy_intp dim[2];
  PyArrayObject* a;
  IMPORTNUMPY;

  if (nfld == 1)
  {
    dim[0] = size;
    a = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, fortran);
    dim[1] = 1;
  }
  else
  {
    if (fortran == 0) { dim[1] = size; dim[0] = nfld; }
    else { dim[0] = size; dim[1] = nfld; }
    a = (PyArrayObject*)PyArray_EMPTY(2, dim, NPY_DOUBLE, fortran);
  }

  E_Float* ptrf = f;
  E_Float* ptra = (E_Float*)PyArray_DATA(a);
  E_Int s = dim[0]*dim[1];
#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < s; i++) ptra[i] = ptrf[i];
  return (PyObject*)a;
}

//=============================================================================
/* Build a numpy array from a FldArray (Int)
   IN: field: Fld field
   IN: fortran=0 (C ordering), fortran=1 (Fortran ordering)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_NUMPY::buildNumpyArray(FldArrayI& field, E_Int fortran)
{
  npy_intp dim[2];
  PyArrayObject* a;
  IMPORTNUMPY;

  if (field.getNfld() == 1)
  {
    dim[0] = field.getSize();
    a = (PyArrayObject*)PyArray_EMPTY(1, dim, E_NPY_INT, fortran);
    dim[1] = 1;
  }
  else
  {
    if (fortran == 0) { dim[1] = field.getSize(); dim[0] = field.getNfld(); }
    else { dim[0] = field.getSize(); dim[1] = field.getNfld(); }
    a = (PyArrayObject*)PyArray_EMPTY(2, dim, E_NPY_INT, fortran);
  }
  E_Int* ptrf = field.begin();
  E_Int* ptra = (E_Int*)PyArray_DATA(a);
  E_Int s = dim[0]*dim[1];
#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < s; i++) ptra[i] = ptrf[i];
  return (PyObject*)a;
}

//=============================================================================
/* Build an empty numpy array
   IN: size, nfld: size and number of fields
   IN: type: type=0 (E_Float), type=1 (E_Int)
   IN: fortran=0 (C ordering), fortran=1 (Fortran ordering)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_NUMPY::buildNumpyArray(E_Int size, E_Int nfld, E_Int type, 
                                   E_Int fortran)
{
  npy_intp dim[2];
  PyArrayObject* a;
  IMPORTNUMPY;
  E_Int ndim = 2;
  if (nfld == 1) { dim[0] = size; ndim = 1; }
  else
  {
    if (fortran == 0) { dim[1] = size; dim[0] = nfld; }
    else { dim[0] = size; dim[1] = nfld; }
    ndim = 2;
  }
  if (type == 0) // E_Float
    a = (PyArrayObject*)PyArray_EMPTY(ndim, dim, NPY_DOUBLE, fortran);
  else
    a = (PyArrayObject*)PyArray_EMPTY(ndim, dim, E_NPY_INT, fortran);
  return (PyObject*)a;
}

//=============================================================================
/* Build a numpy array from a E_Int* + size
   IN: f: pointer (size, nfld)
   IN: size, nfld: sizes of f
   IN: fortran=0 (C ordering), fortran=1 (Fortran ordering)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_NUMPY::buildNumpyArray(E_Int* f, E_Int size, E_Int nfld,
                                   E_Int fortran)
{
  npy_intp dim[2];
  PyArrayObject* a;
  IMPORTNUMPY;

  if (nfld == 1)
  {
    dim[0] = size;
    a = (PyArrayObject*)PyArray_EMPTY(1, dim, E_NPY_INT, fortran);
    dim[1] = 1;
  }
  else
  {
    if (fortran == 0) { dim[1] = size; dim[0] = nfld; }
    else { dim[0] = size; dim[1] = nfld; }
    a = (PyArrayObject*)PyArray_EMPTY(2, dim, E_NPY_INT, fortran);
  }

  E_Int* ptrf = f;
  E_Int* ptra = (E_Int*)PyArray_DATA(a);
  E_Int s = dim[0]*dim[1];
#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < s; i++) ptra[i] = ptrf[i];
  return (PyObject*)a;
}

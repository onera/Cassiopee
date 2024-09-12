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

#include "kcore.h"
#include <stdlib.h>

#ifdef NPY_1_7_API_VERSION
#define FLAG(arr)                                                   \
  PyArray_ENABLEFLAGS((PyArrayObject*)arr, NPY_ARRAY_OWNDATA);      \
  PyArray_ENABLEFLAGS((PyArrayObject*)arr, NPY_ARRAY_F_CONTIGUOUS); \
  PyArray_CLEARFLAGS((PyArrayObject*)arr, NPY_ARRAY_C_CONTIGUOUS);
#define FLAG2(arr)                              \
  PyArray_CLEARFLAGS((PyArrayObject*)arr, NPY_ARRAY_OWNDATA); \
  PyArray_ENABLEFLAGS((PyArrayObject*)arr, NPY_ARRAY_F_CONTIGUOUS);     \
  PyArray_CLEARFLAGS((PyArrayObject*)arr, NPY_ARRAY_C_CONTIGUOUS);
#else
#define FLAG(arr)                                               \
  ((PyArrayObject*)arr)->flags |= NPY_OWNDATA;                  \
  ((PyArrayObject*)arr)->flags |= NPY_F_CONTIGUOUS;             \
  ((PyArrayObject*)arr)->flags ^= NPY_C_CONTIGUOUS;
#define FLAG2(arr)                                              \
  ((PyArrayObject*)arr)->flags |= NPY_F_CONTIGUOUS;             \
  ((PyArrayObject*)arr)->flags ^= NPY_C_CONTIGUOUS;             \
  ((PyArrayObject*)arr)->flags ^= NPY_OWNDATA;                  \
  ((PyArrayObject*)arr)->flags ^= NPY_OWNDATA;
#endif

//=============================================================================
// Threaded memcpy
//==============================================================================
void K_KCORE::memcpy__(E_Int* a, E_Int* b, E_Int s)
{
  if (s < 100) memcpy(a, b, s*sizeof(E_Int));
  else
  {
#pragma omp parallel for default (shared)
    for (E_Int i = 0; i < s; i++) a[i] = b[i];
  }
}
void K_KCORE::memcpy__(E_Float* a, E_Float* b, E_Int s)
{
  if (s < 100) memcpy(a, b, s*sizeof(E_Float));
  else
  {
#pragma omp parallel for default (shared)
    for (E_Int i = 0; i < s; i++) a[i] = b[i];
  }
}

//=============================================================================
// Equivalent de la fonction empty de numpy mais avec des tableaux 
// de double fortrans alignes sur align bits
// IN: shape: un tuple (120,3)
// IN: align: nbre octets pour l'alignement
//=============================================================================
PyObject* K_KCORE::empty(PyObject* self, PyObject* args)
{
  PyObject* shape; E_Int align;
  if (!PYPARSETUPLE_(args, O_ I_, &shape, &align)) return NULL;
  IMPORTNUMPY;
  npy_intp dims[3];
  E_Int nd; E_Int size;
  if (PyTuple_Check(shape)) 
  {
    nd = PyTuple_Size(shape);
    size = 1;
    for (E_Int i = 0; i < nd; i++)
    {
      PyObject* s = PyTuple_GetItem(shape, i);
      E_Int sl = PyLong_AsLong(s);
      size *= sl;
      dims[i] = sl;
    }
  }
  else
  {
    nd = 1;
    size = PyLong_AsLong(shape);
    dims[0] = size;
  }

  // buffer surdimensionne pour l'alignement
  char* buf = (char*)malloc((8*size+align)*sizeof(char));
  if (!buf) return PyErr_NoMemory();
  
  double* mymem = (double*)buf;
  E_LONG addr = (E_LONG)mymem; // adresse en bytes
  //printf("%ld < %d\n", addr%align, align);
  PyObject* arr;
  if (align == 0 || addr%align == 0)
  {
    arr = PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, mymem);
    Py_INCREF(arr);
    FLAG(arr);
  }
  else
  {
    // Premier numpy (possesseur)
    PyObject* arr1 = PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, mymem);
    Py_INCREF(arr1);
    FLAG(arr1);

    // Le deuxieme non-possesseur aligne
    E_LONG addr2 = addr + addr%align;
    double* data = (double*)addr2;
    arr = PyArray_SimpleNewFromData(nd, dims, NPY_DOUBLE, data);
    Py_INCREF(arr);
    FLAG2(arr);
  }
  //free(buf);
  return arr;
}

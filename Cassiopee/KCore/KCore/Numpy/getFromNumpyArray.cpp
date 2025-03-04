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

#ifdef NPY_1_7_API_VERSION
#define GETDIMS E_Int isFortran = PyArray_IS_F_CONTIGUOUS(a); \
  if (dim == 1) { nfld = 1; size = PyArray_DIMS(a)[0]; } \
  else if (dim == 2) { \
    if (isFortran == 0) { nfld = PyArray_DIMS(a)[0]; size = PyArray_DIMS(a)[1]; } \
    else { nfld = PyArray_DIMS(a)[1]; size = PyArray_DIMS(a)[0]; } } \
  else { size = PyArray_SIZE(a); nfld = 1; }
#else
#define GETDIMS E_Int isFortran = PyArray_CHKFLAGS(a, NPY_F_CONTIGUOUS); \
  if (dim == 1) { nfld = 1; size = PyArray_DIMS(a)[0]; } \
  else if (dim == 2) { \
    if (isFortran == 0) { nfld = PyArray_DIMS(a)[0]; size = PyArray_DIMS(a)[1]; } \
    else { nfld = PyArray_DIMS(a)[1]; size = PyArray_DIMS(a)[0]; } } \
  else { size = PyArray_SIZE(a); nfld = 1; }
#endif

#ifdef NPY_1_7_API_VERSION
#define GETDIMS2 E_Int isFortran = PyArray_IS_F_CONTIGUOUS(a); \
  if (dim == 1) { nfld = 1; size = PyArray_DIMS(a)[0]; } \
  else if (dim == 2) { \
    if (isFortran == 0) { nfld = PyArray_DIMS(a)[0]; size = PyArray_DIMS(a)[1]; } \
    else { nfld = PyArray_DIMS(a)[1]; size = PyArray_DIMS(a)[0]; } \
    if (size == 1) { size = nfld; nfld = 1; } } \
  else { size = PyArray_SIZE(a); nfld = 1; }
#else
#define GETDIMS2 E_Int isFortran = PyArray_CHKFLAGS(a, NPY_F_CONTIGUOUS); \
  if (dim == 1) { nfld = 1; size = PyArray_DIMS(a)[0]; }                          \
  else if (dim == 2) { \
    if (isFortran == 0) { nfld = PyArray_DIMS(a)[0]; size = PyArray_DIMS(a)[1]; } \
    else { nfld = PyArray_DIMS(a)[1]; size = PyArray_DIMS(a)[0]; } \
    if (size == 1) { size = nfld; nfld = 1; } } \
  else { size = PyArray_SIZE(a); nfld = 1; }
#endif

//=============================================================================
// Retourne un FldArrayI a partir d'un numpy d'entiers
// Si shared=True, le FldArrayI pointe sur le tableau python. Il faut ensuite
// utiliser RELEASESHAREDN.
// Sinon, c'est une copie.
// Le FldArrayI map les dimensions si elles sont egales a 2 (npts, nfld).
// sinon, retourne un FldArrayI (npts, 1).
// Retourne 0 (FAIL), 1 (SUCCESS)
//=============================================================================
E_Int K_NUMPY::getFromNumpyArray(PyObject*o , FldArrayI*& f, E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  E_Int dim = PyArray_NDIM(a);
  E_Int size = 0; E_Int nfld = 1;
  GETDIMS;
  if (shared == false) // copy du numpy
    f = new FldArrayI(size, nfld, (E_Int*)PyArray_DATA(a), false);
  else // partage
  {
    Py_INCREF(o);
    f = new FldArrayI(size, nfld, (E_Int*)PyArray_DATA(a), true);
  }
  return 1;
}

//=============================================================================
// Identical to getFromNumpArray but returns a flat array for (1,nb) arrays
// for PointList
//=============================================================================
E_Int K_NUMPY::getFromPointList(PyObject*o , FldArrayI*& f, E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  E_Int dim = PyArray_NDIM(a);
  E_Int size = 0; E_Int nfld = 1;
  GETDIMS2;
  if (shared == false) // copy du numpy
    f = new FldArrayI(size, nfld, (E_Int*)PyArray_DATA(a), false);
  else // partage
  {
    Py_INCREF(o);
    f = new FldArrayI(size, nfld, (E_Int*)PyArray_DATA(a), true);
  }
  return 1;
}

//=============================================================================
// Retourne un FldArray a partir d'un numpy
// IN: o: objet numpy array
// OUT: f: FldArrayF
// IN: shared: 1 (partage avec python), 0 (copie)
// Dans le cas shared, il faut utiliser RELEASESHAREDN.
// Le FldArrayF map les dimensions si elles sont egales a 2 (npts, nfld).
// sinon, retourne un FldArrayF (npts, 1).
// Retourne 0 (FAIL), 1 (SUCCESS)
//=============================================================================
E_Int K_NUMPY::getFromNumpyArray(PyObject* o, FldArrayF*& f, E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  E_Int dim = PyArray_NDIM(a);
  E_Int size = 0; E_Int nfld = 1;
  GETDIMS;
  if (shared == false) // copy du numpy
    f = new FldArrayF(size, nfld, (E_Float*)PyArray_DATA(a), false);
  else // partage
  {
    Py_INCREF(o);
    f = new FldArrayF(size, nfld, (E_Float*)PyArray_DATA(a), true);
  }
  return 1;
}

//=============================================================================
// Retourne un E_Float* a partir d'un numpy
// IN: o: objet numpy array
// OUT: f: ptr sur le tableau (size, nfld)
// OUT: size: taille du tableau
// OUT: nfld: taille du tableau
// IN: shared: 1 (partage avec python), 0 (copie)
// Si shared, il faut utiliser Py_DECREF sur le numpy a la fin.
// Les dimensions (size,nfld) sont corrects si le numpy a 2 dimensions.
// Sinon, retourne (size,1).
// Retourne 0 (FAIL), 1 (SUCCESS)
//=============================================================================
E_Int K_NUMPY::getFromNumpyArray(PyObject* o, E_Float*& f, E_Int& size,
                                 E_Int& nfld, E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  E_Int dim = PyArray_NDIM(a);
  size = 0; nfld = 1;
  GETDIMS;
  if (shared == false) // copie du numpy
  {
    f = new E_Float [size*nfld];
    E_Float* ptr = (E_Float*)PyArray_DATA(a);
#pragma omp parallel for
    for (E_Int i = 0; i < size*nfld; i++) f[i] = ptr[i];
  }
  else // partage
  {
    Py_INCREF(o);
    f = (E_Float*)PyArray_DATA(a);
  }
  return 1;
}

//=============================================================================
// Identique au precedent, mais retourne uniquement la taille a plat
E_Int K_NUMPY::getFromNumpyArray(PyObject* o, E_Float*& f, E_Int& size,
                                 E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  size = PyArray_SIZE(a);
  if (shared == false) // copie du numpy
  {
    f = new E_Float [size];
    E_Float* ptr = (E_Float*)PyArray_DATA(a);
#pragma omp parallel for
    for (E_Int i = 0; i < size; i++) f[i] = ptr[i];
  }
  else // partage
  {
    Py_INCREF(o);
    f = (E_Float*)PyArray_DATA(a);
  }
  return 1;
}

//=============================================================================
// Retourne un E_Int* a partir d'un numpy
// IN: o: objet numpy array
// OUT: f: ptr sur le tableau (size, nfld)
// OUT: size: taille du tableau
// OUT: nfld: taille du tableau
// IN: shared: 1 (partage avec python), 0 (copie)
// Si shared, il faut utiliser Py_DECREF sur le numpy a la fin.
// Les dimensions (size,nfld) sont corrects si le numpy a 2 dimensions.
// Sinon, retourne (size,1).
// Retourne 0 (FAIL), 1 (SUCCESS)
//=============================================================================
E_Int K_NUMPY::getFromNumpyArray(PyObject* o, E_Int*& f, E_Int& size,
                                 E_Int& nfld, E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  E_Int dim = PyArray_NDIM(a);
  size = 0; nfld = 1;
  GETDIMS;
  if (shared == false) // copie du numpy
  {
    f = new E_Int [size*nfld];
    E_Int* ptr = (E_Int*)PyArray_DATA(a);
#pragma omp parallel for
    for (E_Int i = 0; i < size*nfld; i++) f[i] = ptr[i];
  }
  else // partage
  {
    Py_INCREF(o);
    f = (E_Int*)PyArray_DATA(a);
  }
  return 1;
}

//=============================================================================
// Identique au precedent mais retourne la taille a plat
E_Int K_NUMPY::getFromNumpyArray(PyObject* o, E_Int*& f, E_Int& size,
                                 E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  size = PyArray_SIZE(a);
  if (shared == false) // copie du numpy
  {
    f = new E_Int [size];
    E_Int* ptr = (E_Int*)PyArray_DATA(a);
#pragma omp parallel for
    for (E_Int i = 0; i < size; i++) f[i] = ptr[i];
  }
  else // partage
  {
    Py_INCREF(o);
    f = (E_Int*)PyArray_DATA(a);
  }
  return 1;
}

// identical to getFromNumpyArray, but if array is (1,nb) return a 
// flat array (specific for pointList)
E_Int K_NUMPY::getFromPointList(PyObject* o, E_Int*& f, E_Int& size,
                                E_Int& nfld, E_Boolean shared)
{
  IMPORTNUMPY;
  if (PyArray_Check(o) == 0) return 0;
  PyArrayObject* a = (PyArrayObject*)o;
  E_Int dim = PyArray_NDIM(a);
  size = 0; nfld = 1;
  GETDIMS2;
  if (shared == false) // copie du numpy
  {
    f = new E_Int [size*nfld];
    E_Int* ptr = (E_Int*)PyArray_DATA(a);
#pragma omp parallel for
    for (E_Int i = 0; i < size*nfld; i++) f[i] = ptr[i];
  }
  else // partage
  {
    Py_INCREF(o);
    f = (E_Int*)PyArray_DATA(a);
  }
  return 1;
}

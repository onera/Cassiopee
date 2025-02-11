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
# include "Array/Array.h"
using namespace K_FLD;

//=============================================================================
// Extrait les donnees utiles d'un objet python o
// Cet objet python peut etre:
// - une liste python d'entiers
// - un numpy array d'entiers
// Retourne une FldArrayI alloue et rempli a plat (une copie)
// Retourne 0 (FAIL), 1 (SUCCESS)
//=============================================================================
E_Int K_ARRAY::getFromList(PyObject* o, FldArrayI& out)
{
  E_Int val; E_Float valf;
  IMPORTNUMPY;
  if (PyList_Check(o) == true)
  {
    E_Int n = PyList_Size(o);
    if (n == 0) return 0; // nothing in list
    out.malloc(n);
    PyObject* first = PyList_GetItem(o, 0);

    if (PyLong_Check(first) == true || PyInt_Check(first) == true)
    {
      for (E_Int i = 0; i < n; i++)
      {
        val = PyLong_AsLong(PyList_GetItem(o, i));
        out[i] = val;
      }
      return 1;
    }
    
    val = PyArray_PyIntAsInt(first);
    if ((val != -1) || (not PyErr_Occurred()))
    {
      for (E_Int i = 0; i < n; i++)
      {
        val = PyArray_PyIntAsInt(PyList_GetItem(o, i));
        out[i] = val;
      }
      return 1;
    }

    if (PyFloat_Check(first) == true)
    {
      for (E_Int i = 0; i < n; i++)
      {
        valf = PyFloat_AsDouble(PyList_GetItem(o, i));
        out[i] = (E_Int)valf;
      } 
      return 1;
    }
    return 0;
  }
  else if (PyArray_Check(o) == true)
  {
    //PyArrayObject* ac = (PyArrayObject*)
    //  PyArray_ContiguousFromObject(o, NPY_INT, 1, 10000000);
    PyArrayObject* ac = (PyArrayObject*)o; Py_INCREF(ac);
    if (ac == NULL) return 0;
    E_Int nd = PyArray_NDIM(ac);
    if (nd < 1) return 0;
    E_Int n = 1;
    for (E_Int i = 0; i < nd; i++) n *= PyArray_DIMS(ac)[i];
    out.malloc(n);
    E_Int* ptr = (E_Int*)PyArray_DATA(ac);
    for (E_Int i = 0; i < n; i++) out[i] = ptr[i];
    Py_DECREF(ac);
    return 1;
  }
  else return 0;
}

//=============================================================================
// Extrait les donnees utiles d'un objet python o
// Cet objet python peut etre :
// - une liste python de double
// - un numpy array de double
// Retourne une FldArrayF alloue et rempli a plat
// Retourne 0 (FAIL), 1 (SUCCESS)
//=============================================================================
E_Int K_ARRAY::getFromList(PyObject* o, FldArrayF& out)
{
  E_Float val; E_Int vali;
  IMPORTNUMPY;
  if (PyList_Check(o) == true)
  {
    E_Int n = PyList_Size(o);
    out.malloc(n);
    if (PyFloat_Check(PyList_GetItem(o, 0)) == true)
    {
      for (E_Int i = 0; i < n; i++)
      {
        val = PyFloat_AsDouble(PyList_GetItem(o, i));
        out[i] = val;
      }
    }
    else // suppose int
    {
      for (E_Int i = 0; i < n; i++)
      {
        vali = PyLong_AsLong(PyList_GetItem(o, i));
        out[i] = (E_Float)vali;
      } 
    }
    return 1;
  }
  else if (PyArray_Check(o) == true)
  {
    //PyArrayObject* ac = (PyArrayObject*)
    //  PyArray_ContiguousFromObject(o, NPY_DOUBLE, 1, 10000000);
    PyArrayObject* ac = (PyArrayObject*)o; Py_INCREF(ac);
    if (ac == NULL) return 0;
    E_Int nd = PyArray_NDIM(ac);
    if (nd < 1) return 0;
    E_Int n = 1;
    for (E_Int i = 0; i < nd; i++) n *= PyArray_DIMS(ac)[i];
    out.malloc(n);
    E_Float* ptr = (E_Float*)PyArray_DATA(ac);
    for (E_Int i = 0; i < n; i++) out[i] = ptr[i];
    Py_DECREF(ac);
    return 1;
  }
  else return 0;
}

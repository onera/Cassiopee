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
#include <string.h>
#include "Array/Array.h"
#include "String/kstring.h"

using namespace K_FLD;

E_Int __check_array(PyObject* o, PyArrayObject*& a, char*& varString)
{
  PyObject* tpl;
  if (PyList_Check(o) == false)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  E_Int size = PyList_Size(o);
  if (size != 4 && size != 5)
  {
    PyErr_Warn(PyExc_Warning, 
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -2;
  }
 
  // -- varString --
  PyObject* l = PyList_GetItem(o,0);
  if (PyString_Check(l))
  {
    // pointeur sur la chaine python
    varString = PyString_AsString(l);  
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(l))
  {
    varString = (char*)PyUnicode_AsUTF8(l); 
  }
#endif
  else
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. First element must be a string.");
    return -3;
  }

  // -- field --
  tpl = PyList_GetItem(o, 1);
  if (PyArray_Check(tpl) == false)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: second arg in array must be a numpy array.");
    return -4;
  }
  //a = (PyArrayObject*)PyArray_ContiguousFromObject(tpl, NPY_DOUBLE,
  //                                                 1, 10000000);
  a = (PyArrayObject*)tpl; Py_INCREF(a);

  if (a == NULL)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: second arg must be a numpy array.");
    return -4;
  }
  
  if (PyArray_NDIM(a) != 2)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: field must have two dimensions.");
    Py_DECREF(a);
    return -4;
  }
  return 0;
}

E_Int __get_connectivity(PyObject*o, PyArrayObject*& a, 
                         E_Int& ni, E_Int& nj, E_Int& nk,
                         DynArray<E_Int>& c,
                         char*& eltType)
{
  PyObject* tpl;
  PyArrayObject* ac;

  E_Int size = PyList_Size(o);
  if (size == 4) // unstruct array
  {
    // -- element type --
    PyObject* l = PyList_GetItem(o,3);
    if (PyString_Check(l))
    {
      eltType = PyString_AsString(l);  
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l))
    {
      eltType = (char*)PyUnicode_AsUTF8(l); 
    }
#endif
    else
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: an unstruct array must be of list of type ['vars', a, c, 'ELTTYPE']. Last element must be a string.");
      Py_DECREF(a);
      return -7;
    }

    if (K_STRING::cmp(eltType, "NODE") != 0 &&
        K_STRING::cmp(eltType, "BAR") != 0 &&
        K_STRING::cmp(eltType, "TRI") !=0 &&
        K_STRING::cmp(eltType, "QUAD") != 0 &&
        K_STRING::cmp(eltType, "TETRA") !=0 &&
        K_STRING::cmp(eltType, "PYRA") != 0 &&
        K_STRING::cmp(eltType, "PENTA") != 0 &&
        K_STRING::cmp(eltType, "HEXA") != 0 &&
        K_STRING::cmp(eltType, "NGON") != 0 &&
        K_STRING::cmp(eltType, "NODE*") != 0 &&
        K_STRING::cmp(eltType, "BAR*") != 0 &&
        K_STRING::cmp(eltType, "TRI*") !=0 &&
        K_STRING::cmp(eltType, "QUAD*") != 0 &&
        K_STRING::cmp(eltType, "TETRA*") !=0 &&
        K_STRING::cmp(eltType, "PYRA*") != 0 &&
        K_STRING::cmp(eltType, "PENTA*") != 0 &&
        K_STRING::cmp(eltType, "HEXA*") != 0 &&
        K_STRING::cmp(eltType, "NGON*") != 0)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: element type unknown: %s. Must be in NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON or NODE*, BAR*, TRI*, QUAD*, TETRA*, PYRA*, PENTA*, HEXA*, NGON*.");
      Py_DECREF(a);
      return -7;
    }
    
    // temporary hack to preserve exisitng code : need a rethink of indexation start (0 and/or 1)
    E_Int shift=1;
    if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0 )
      shift=0;
    
    // -- connectivity --
    tpl = PyList_GetItem(o, 2);
    if (PyArray_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning, 
                 "getFromArray: third arg in array must be a numpy array.");
      Py_DECREF(a);
      return -6;
    }
    //ac = (PyArrayObject*)PyArray_ContiguousFromObject(tpl, E_NPY_INT,
    //                                                  1, 10000000);
    ac = (PyArrayObject*)tpl; Py_INCREF(ac);

    if (ac == NULL)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: third arg must be a numpy array.");
      Py_DECREF(a);
      return -6;
    }
  
    if (PyArray_NDIM(ac) != 2)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: connectivity must have two dimensions.");
      Py_DECREF(a); Py_DECREF(ac);
      return -4;
    }
    E_Int s = PyArray_DIMS(ac)[1];
    E_Int nfld = PyArray_DIMS(ac)[0];
    c.resize(nfld,s);
    E_Int* d = (E_Int*)PyArray_DATA(ac);
    DynArray<E_Int>::iterator it = c.begin();
    for (E_Int i = 0; i < s; i++)
      for (E_Int n = 0; n < nfld; n++)
      {
        *it = d[i + s*n]-shift; it++;
      }
    
    Py_DECREF(a); Py_DECREF(ac);
    return 2;
  }
  else // struct array
  {
    tpl = PyList_GetItem(o,2);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: third arg must be an integer.");
      Py_DECREF(a);
      return -5;
    }
    ni = PyLong_AsLong(tpl);
    tpl = PyList_GetItem(o,3);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: fourth arg must be an integer.");
      Py_DECREF(a);
      return -5;
    }
    nj = PyLong_AsLong(tpl);
    tpl = PyList_GetItem(o,4);
    if (PyLong_Check(tpl) == false && PyInt_Check(tpl) == false)
    {
      PyErr_Warn(PyExc_Warning,
                 "getFromArray: fifth arg must be an integer.");
      Py_DECREF(a);
      return -5;
    }
    nk = PyLong_AsLong(tpl);
    Py_DECREF(a);
    return 1;
  }
}
 
//=============================================================================
// Extrait les donnees utiles d'un objet python struct array
// defini par: [ 'vars, a, ni, nj, nk ]
// ou d'un objet python unstruct array
// defini par: [ 'vars', a, c, "ELTTYPE"].
// ou ELTTYPE vaut: NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON
// avec ou sans star.
// return 1: valid struct array
//            f: field en stockage Dyn (champ de v1,v2,v3,...)
//            ni, nj, nk: number of points
//            varString
// return 2: valid unstruct array
//            f: field en stockage Dyn (champ de v1,v2,v3,...)
//            c: connectivity en stockage Dyn (champ c1,c2,c3 avec indices
//               commencants a zero).
//            eltType: type of elements
//            varString
// return -1: given object is not a list.
// return -2: not a valid number of elts in list.
// return -3: first element is not a var string.
// return -4: a is not a valid numpy array.
// return -5: array is structured but ni, nj, nk unvalid.
// return -6: array is unstructured but connectivity is unvalid.
// return -7: array is unstructured but elt type is unknown.
// C'est la responsabilite de l'appelant de liberer la memoire de f et 
// eventuellement de c
//=============================================================================
E_Int K_ARRAY::getFromArray(PyObject* o,
                            char*& varString,
                            DynArray<E_Float>*& f,
                            E_Int& ni, E_Int& nj, E_Int& nk,
                            DynArray<E_Int>*& c,
                            char*& eltType)
{
  IMPORTNUMPY;
  
  PyArrayObject* a;
  E_Int res = __check_array(o, a, varString);
  if (res != 0) return res;
  
  E_Int s = PyArray_DIMS(a)[1];
  E_Int nfld = PyArray_DIMS(a)[0];
  E_Int nvar = getNumberOfVariables(varString);

  if (nfld != nvar)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: number of variables different in varString and field.");
    Py_DECREF(a);
    return -4;
  }
  
  f = new DynArray<E_Float>(nfld, s);
  E_Float* d = (E_Float*)PyArray_DATA(a);
  DynArray<E_Float>::iterator it = f->begin();
  for (E_Int i = 0; i < s; i++)
    for (E_Int n = 0; n < nfld; n++)
    {
      *it = d[i + s*n]; it++;
    }

  c = new DynArray<E_Int>(1, 1);//will be resized inside
  return __get_connectivity(o, a, ni, nj, nk, *c, eltType);
}

// same as above ignoring fields other than coordinates
E_Int K_ARRAY::getFromArray(PyObject* o,
                                char*& varString,
                                DynArray<E_Float>& f,
                                E_Int& ni, E_Int& nj, E_Int& nk,
                                DynArray<E_Int>& c,
                                char*& eltType)
{
  IMPORTNUMPY;
  
  PyArrayObject* a;
  E_Int res = __check_array(o, a, varString);
  if (res != 0) return res;
  
  E_Int s = PyArray_DIMS(a)[1];
  E_Int nfld = PyArray_DIMS(a)[0];
  E_Int nvar = getNumberOfVariables(varString);

  if (nfld != nvar)
  {
    PyErr_Warn(PyExc_Warning,
               "getFromArray: number of variables different in varString and field.");
    Py_DECREF(a);
    return -4;
  }
  
  f.resize(nfld, s);
  
  E_Float* d = (E_Float*)PyArray_DATA(a);
  DynArray<E_Float>::iterator it = f.begin();
  for (E_Int i = 0; i < s; i++)
    for (E_Int n = 0; n < nfld; n++)
    {
      *it = d[i + s*n]; it++;
    }

  return __get_connectivity(o, a, ni, nj, nk, c, eltType);
}

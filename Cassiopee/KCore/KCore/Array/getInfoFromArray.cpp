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
#include "Array/Array.h"
#include "String/kstring.h"

//=============================================================================
/* 
   Retourne toutes les infos contenues dans un array,
   a l'exception des tableaux eux-memes.
   Retourne le type (1: structure, 2: non structure).
   Si type=1, ni, nj, nk sont renseignes.
   Si type=2, nvertex, nelt, eltType, sizeConect sont renseignes.
   Cette routine ne fait pas de verification.
*/
//=============================================================================
/*
E_Int K_ARRAY::getInfoFromArray(PyObject* o, char*& varString,
                                E_Int& ni, E_Int& nj, E_Int& nk,
                                E_Int& nvertex, E_Int& nelt, 
                                E_Int& sizeConnect, char*& eltType)
{
  if (PyList_Check(o) == 0)
  {
    PyErr_Warn(PyExc_Warning,
               "getInfoFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  E_Int s = PyList_Size(o);
  if (s == 5) // structure 
  {
    varString = PyString_AsString(PyList_GetItem(o,0));
    ni = PyLong_AsLong(PyList_GetItem(o, 2));
    nj = PyLong_AsLong(PyList_GetItem(o, 3));
    nk = PyLong_AsLong(PyList_GetItem(o, 4));
    return 1;
  }
  else if (s == 4) // non structure
  {
    varString = PyString_AsString(PyList_GetItem(o,0));
    IMPORTNUMPY;
    nvertex = PyArray_DIM((PyArrayObject*)PyList_GetItem(o,1), 1);
    
    eltType = PyString_AsString(PyList_GetItem(o,3));
    PyObject* c = PyList_GetItem(o,2);

    //PyArrayObject* ac = 
    //  (PyArrayObject*)PyArray_ContiguousFromObject(c, E_NPY_INT,
    //						   1, 10000000);
    PyArrayObject* ac = (PyArrayObject*)c; Py_INCREF(ac);
    if (K_STRING::cmp(eltType, "NGON") == 0)
    {
      E_Int* ptr = (E_Int*)PyArray_DATA(ac);
      nelt = ptr[ptr[1]+2];
      sizeConnect = PyArray_DIM(ac, 1);
      Py_DECREF(ac);
    }
    else
    {
      nelt = PyArray_DIM(ac, 1);
      sizeConnect = nelt*PyArray_DIM(ac, 0);
      Py_DECREF(ac);
    }
    return 2;
  }
  else
  {
    PyErr_Warn(PyExc_Warning, 
               "getInfoFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -1;
  }
}
*/
//=============================================================================
// marche sur array et array2 (la seule a conserver)
E_Int K_ARRAY::getInfoFromArray(PyObject* o, char*& varString,
                                E_Int& ni, E_Int& nj, E_Int& nk,
                                E_Int& nvertex, E_Int& nelt, 
                                E_Int& sizeConnect, char*& eltType)
{
  if (PyList_Check(o) == 0)
  {
    PyErr_Warn(PyExc_Warning,
               "getInfoFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  E_Int s = PyList_Size(o);
  if (s == 5) // structure 
  {
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
    else
    {
      PyErr_Warn(PyExc_Warning,
                 "getInfoFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
      return -1;
    }
#endif
    ni = PyLong_AsLong(PyList_GetItem(o, 2));
    nj = PyLong_AsLong(PyList_GetItem(o, 3));
    nk = PyLong_AsLong(PyList_GetItem(o, 4));
    return 1;
  }
  else if (s == 4) // non structure
  {
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
    //varString = PyString_AsString(PyList_GetItem(o,0));
    IMPORTNUMPY;
    PyObject* f = PyList_GetItem(o,1);
    if (PyList_Check(f) == true) // array2
      nvertex = PyArray_DIM((PyArrayObject*)PyList_GetItem(f, 0), 0);
    else
      nvertex = PyArray_DIM((PyArrayObject*)f, 1);
    
    PyObject* e = PyList_GetItem(o,3);
    if (PyString_Check(e))
    {
      // pointeur sur la chaine python
      eltType = PyString_AsString(e);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l))
    {
      eltType = (char*)PyUnicode_AsUTF8(e); 
    }
#endif
    else
    {
      PyErr_Warn(PyExc_Warning,
                 "getInfoFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
      return -1;
    }
    PyObject* c = PyList_GetItem(o,2);

    //PyArrayObject* ac = 
    //  (PyArrayObject*)PyArray_ContiguousFromObject(c, E_NPY_INT,
    //               1, 10000000);  
    if (K_STRING::cmp(eltType, 4, "NGON") == 0)
    {
      if (PyList_Check(c) == true) // array2
      {
        PyArrayObject* acl = (PyArrayObject*)PyList_GetItem(c, 2);
        nelt = PyArray_DIM(acl, 0);
        sizeConnect = 0;
      }
      else // array
      {
        PyArrayObject* ac = (PyArrayObject*)c; Py_INCREF(ac);
        E_Int* ptr = (E_Int*)PyArray_DATA(ac);
        nelt = ptr[ptr[1]+2];
        sizeConnect = PyArray_DIM(ac, 1);
        Py_DECREF(ac);
      }  
    }
    else // elements basiques
    {
      if (PyList_Check(c) == true) // array2
      {
        if (PyList_Size(c) == 0) // NODE
        {
          nelt = 0;
          sizeConnect = 0;
        }
        else
        {
          PyArrayObject* acl = (PyArrayObject*)PyList_GetItem(c, 0);
          nelt = PyArray_DIM(acl, 1);
          sizeConnect = nelt*PyArray_DIM(acl, 0);
        }
      }
      else // array
      {
        PyArrayObject* ac = (PyArrayObject*)c; Py_INCREF(ac);
        nelt = PyArray_DIM(ac, 1);
        sizeConnect = nelt*PyArray_DIM(ac, 0);
        Py_DECREF(ac);
      }
    }
    return 2;
  }
  else
  {
    PyErr_Warn(PyExc_Warning, 
               "getInfoFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -1;
  }
}

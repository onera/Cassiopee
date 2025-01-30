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

//=============================================================================
/* 
   Retourne le nombre de pts dans un array et le type (1: structure,
   2: non structure).
   Routine rapide.
*/
//=============================================================================
E_Int K_ARRAY::getSizeFromArray(PyObject* o, E_Int& size, E_Int& type)
{
  type = 0; size = 0;
  if (PyList_Check(o) == 0)
  {
    PyErr_Warn(PyExc_Warning,
               "getSizeFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  E_Int s = PyList_Size(o);
  if (s == 5) // structure 
  {
    type = 1;
    E_Int ni = PyLong_AsLong(PyList_GetItem(o,2));
    E_Int nj = PyLong_AsLong(PyList_GetItem(o,3));
    E_Int nk = PyLong_AsLong(PyList_GetItem(o,4));
    size = ni * nj * nk;
    return 1;
  }
  else if (s == 4) // non structure
  {
    IMPORTNUMPY;
    type = 2;
    size = PyArray_DIM((PyArrayObject*)PyList_GetItem(o,1), 1);
    return 1;
  }
  else
  {
    PyErr_Warn(PyExc_Warning, 
               "getSizeFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -2;
  }
}

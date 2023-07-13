/*    
    Copyright 2013-2023 Onera.

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

# include "xcore.h"

// Chunk2part for NGON2

//============================================================================
/* IN: Chunk of NGON2 */
//============================================================================
PyObject* K_XCORE::chunk2part(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  if (!PyArg_ParseTuple(args, "O", &arrays)) return NULL;
  
  E_Int nzones = PyList_Size(arrays);

  PyObject* o; PyObject* l;
  E_Float* ptrf; E_Int size; E_Int nfld;
  E_Int* ptri;
  
  for (E_Int i = 0; i < nzones; i++)
  {
    l = PyList_GetItem(arrays, i);

    // 1 must be coordinateX chunk
    o = PyList_GetItem(l, 0);
    K_NUMPY::getFromNumpyArray(o, ptrf, size, nfld, true);
    // 2 must be coordinateY chunk
    o = PyList_GetItem(l, 1);
    K_NUMPY::getFromNumpyArray(o, ptrf, size, nfld, true);
    // 3 must be coordinateZ chunk
    o = PyList_GetItem(l, 2);
    K_NUMPY::getFromNumpyArray(o, ptrf, size, nfld, true);
    // 4 must be ngon chunk
    o = PyList_GetItem(l, 3);
    K_NUMPY::getFromNumpyArray(o, ptri, size, nfld, true);
    // 5 must be ngon so chunk
    o = PyList_GetItem(l, 4);
    K_NUMPY::getFromNumpyArray(o, ptri, size, nfld, true);
    // 6 must be nface chunk
    o = PyList_GetItem(l, 5);
    K_NUMPY::getFromNumpyArray(o, ptri, size, nfld, true);
    // 6 must be nface so chunk
    o = PyList_GetItem(l, 6);
    K_NUMPY::getFromNumpyArray(o, ptri, size, nfld, true);
    
    // PE a venir...

  }

  // ..
  
  // export with buildNumpyArray
  
  // Release numpys
  for (E_Int i = 0; i < nzones; i++) 
  {
    l = PyList_GetItem(arrays, i);
    Py_DECREF(PyList_GetItem(l, 0));
    Py_DECREF(PyList_GetItem(l, 1));
    Py_DECREF(PyList_GetItem(l, 2));
    Py_DECREF(PyList_GetItem(l, 3));
    Py_DECREF(PyList_GetItem(l, 4));
    Py_DECREF(PyList_GetItem(l, 5));
    Py_DECREF(PyList_GetItem(l, 6));
    
  }

  Py_INCREF(Py_None);
  return Py_None;
      
}

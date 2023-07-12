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
  PyObject* coords; // list of coord numpys (chunks)
  PyObject* ngons; // list of ngons numpys (chunks)
  PyObject* nfaces; // list of nfaces numpys (chunks)
  PyObject* pes; // list of pe numpys (chunks)
  if (!PyArg_ParseTuple(args, "OOOO", &coords, &ngons, &nfaces, &pes)) return NULL;
  
  E_Int ncoords = PyList_Size(coords);
  PyObject* o;
  E_Float* ptrf; E_Int size; E_Int nfld;
  E_Int* ptri;
  for (E_Int i = 0; i < ncoords; i++)
  {
    o = PyList_GetItem(coords, i);
    K_NUMPY::getFromNumpyArray(o, ptrf, size, nfld, true);
  }

  E_Int nngons = PyList_Size(ngons);
  for (E_Int i = 0; i < nngons; i++)
  {
    o = PyList_GetItem(ngons, i);
    K_NUMPY::getFromNumpyArray(o, ptri, size, nfld, true);
  }
  

  // Release numpys
  for (E_Int i = 0; i < ncoords; i++) Py_DECREF(PyList_GetItem(coords, i));
  for (E_Int i = 0; i < nngons; i++) Py_DECREF(PyList_GetItem(ngons, i));
  

  Py_INCREF(Py_None);
  return Py_None;
      
}

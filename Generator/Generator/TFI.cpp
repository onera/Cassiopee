/*    
    Copyright 2013-2019 Onera.

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

// TFI generator

# include "generator.h"

using namespace std;
using namespace K_FLD;

// ============================================================================
/* TFI */
// ============================================================================
PyObject* K_GENERATOR::TFIMesh(PyObject* self, PyObject* args)
{
  PyObject* arrays; 
  if (!PyArg_ParseTuple(args, "O", &arrays)) return NULL;
  
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString( PyExc_TypeError, 
                     "TFI: argument must be a list of arrays.");
    return NULL;
  }
  
  E_Int size = PyList_Size(arrays);

  switch ( size ) 
  {
    case 6:
      return TFI3D(arrays);

    case 5:
      return TFIPENTA(arrays);

    case 4:
      for (int i = 0; i < size; i++)
      {
        PyObject* tpl = PyList_GetItem(arrays, i);
        if (PyList_Size(tpl) != 5)
        {
          return TFITETRA(arrays);          
        }
      }
      return TFI2D(arrays);

    case 3:
      return TFITRI(arrays);
          
    default:
      PyErr_SetString(PyExc_TypeError,
                      "TFI: wrong arguments.");
      return NULL;
  }
}

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
#define K_ARRAY_UNIQUE_SYMBOL
#include "rigidMotion.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef PyrigidMotion [] =
{
  {"move", K_RIGIDMOTION::move, METH_VARARGS},
  {"moveN", K_RIGIDMOTION::moveN, METH_VARARGS},
  {NULL, NULL}
};

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
  PyMODINIT_FUNC initrigidMotion();
  PyMODINIT_FUNC initrigidMotion()
  {
    Py_InitModule("rigidMotion", PyrigidMotion);
    import_array();
  }
}

/*    
    Copyright 2013-2018 Onera.

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
#include "initiator.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyinitiator [] =
{
  {"initLamb", K_INITIATOR::initLamb, METH_VARARGS},
  {"initVisbal", K_INITIATOR::initVisbal, METH_VARARGS},
  {"initYee", K_INITIATOR::initYee, METH_VARARGS},
  {"initScully", K_INITIATOR::initScully, METH_VARARGS},
  {"overlayField", K_INITIATOR::overlayField, METH_VARARGS},
  {NULL, NULL}
};

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
  PyMODINIT_FUNC initinitiator();
  PyMODINIT_FUNC initinitiator()
  {
    Py_InitModule("initiator", Pyinitiator);
    import_array();
  }
}

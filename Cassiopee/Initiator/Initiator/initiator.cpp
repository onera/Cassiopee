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
#define K_ARRAY_UNIQUE_SYMBOL
#include "initiator.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyinitiator [] =
{
  {"initLamb", K_INITIATOR::initLamb, METH_VARARGS},
  {"initVisbal", K_INITIATOR::initVisbal, METH_VARARGS},
  {"initWissocq", K_INITIATOR::initWissocq, METH_VARARGS},
  {"initYee", K_INITIATOR::initYee, METH_VARARGS},
  {"initScully", K_INITIATOR::initScully, METH_VARARGS},
  {"overlayField", K_INITIATOR::overlayField, METH_VARARGS},
  {"applyGaussianAL", K_INITIATOR::applyGaussianAL, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "initiator",
        NULL,
        -1,
        Pyinitiator
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_initiator();
  PyMODINIT_FUNC PyInit_initiator()
#else
  PyMODINIT_FUNC initinitiator();
  PyMODINIT_FUNC initinitiator()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("initiator", Pyinitiator);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

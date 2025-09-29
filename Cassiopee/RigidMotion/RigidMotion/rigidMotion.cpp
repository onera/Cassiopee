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
#include "rigidMotion.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef PyrigidMotion [] =
{
  {"move", K_RIGIDMOTION::move, METH_VARARGS},
  {"moveN", K_RIGIDMOTION::moveN, METH_VARARGS},
  {"evalGridMotionN", K_RIGIDMOTION::evalGridMotionN, METH_VARARGS},
  {"evalSpeed3", K_RIGIDMOTION::evalSpeed3, METH_VARARGS},
  {"_computeRotorMotionZ", K_RIGIDMOTION::_computeRotorMotionZ, METH_VARARGS},
  {"_computeRotorMotionInfo", K_RIGIDMOTION::_computeRotorMotionInfo, METH_VARARGS},
  {"copyCoords", K_RIGIDMOTION::copyCoords, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "rigidMotion",
        NULL,
        -1,
        PyrigidMotion
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_rigidMotion();
  PyMODINIT_FUNC PyInit_rigidMotion()
#else
  PyMODINIT_FUNC initrigidMotion();
  PyMODINIT_FUNC initrigidMotion()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("rigidMotion", PyrigidMotion);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

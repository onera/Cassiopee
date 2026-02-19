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
#define K_ARRAY_UNIQUE_SYMBOL
#include "modeler.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef PyModeler[] =
{
  {"exportCAD", K_MODELER::exportCAD, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "modeler",
        NULL,
        -1,
        PyModeler
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_modeler();
  PyMODINIT_FUNC PyInit_modeler()
#else
  PyMODINIT_FUNC initmodeler();
  PyMODINIT_FUNC initmodeler()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("modeler", PyTransform);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

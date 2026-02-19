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
#include "dist2walls.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pydist2walls [] =
{
  {"distance2Walls", K_DIST2WALLS::distance2Walls, METH_VARARGS},
  {"distance2WallsSigned", K_DIST2WALLS::distance2WallsSigned, METH_VARARGS},
  {"distance2WallsOrtho", K_DIST2WALLS::distance2WallsOrtho, METH_VARARGS},
  {"distance2WallsOrthoSigned", K_DIST2WALLS::distance2WallsOrthoSigned, METH_VARARGS},
  {"eikonal", K_DIST2WALLS::eikonal, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "dist2walls",
        NULL,
        -1,
        Pydist2walls
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_dist2walls();
  PyMODINIT_FUNC PyInit_dist2walls()
#else
  PyMODINIT_FUNC initdist2walls();
  PyMODINIT_FUNC initdist2walls()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("dist2walls", Pydist2walls);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

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
#include "geom.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pygeom [] =
{
  {"naca", K_GEOM::nacaMesh, METH_VARARGS},
  {"line", K_GEOM::lineMesh, METH_VARARGS},
  {"circle", K_GEOM::circleMesh, METH_VARARGS},
  {"sphere", K_GEOM::sphereMesh, METH_VARARGS},
  {"cone", K_GEOM::coneMesh, METH_VARARGS},
  {"torus", K_GEOM::torus, METH_VARARGS},
  {"triangle", K_GEOM::triangleMesh, METH_VARARGS},
  {"quadrangle", K_GEOM::quadrangleMesh, METH_VARARGS},
  {"bezier", K_GEOM::bezier, METH_VARARGS},
  {"lineGenerate", K_GEOM::lineGenerateMesh, METH_VARARGS},
  {"lineGenerate2", K_GEOM::lineGenerate2, METH_VARARGS},
  {"addSeparationLine", K_GEOM::addSeparationLineMesh, METH_VARARGS},
  {"axisym", K_GEOM::axisym, METH_VARARGS},
  {"volumeFromCrossSections", K_GEOM::volumeFromCrossSections, METH_VARARGS},
  {"getLength", K_GEOM::getLength, METH_VARARGS},
  {"getDistantIndex", K_GEOM::getDistantIndex, METH_VARARGS},
  {"getNearestPointIndex", K_GEOM::getNearestPointIndex, METH_VARARGS},
  {"getCurvatureAngle", K_GEOM::getCurvatureAngle, METH_VARARGS},
  {"getCurvatureRadius", K_GEOM::getCurvatureRadius, METH_VARARGS},
  {"getCurvatureHeight", K_GEOM::getCurvatureHeight, METH_VARARGS},
  {"getCurvilinearAbscissa", K_GEOM::getCurvilinearAbscissa, METH_VARARGS},
  {"polyline", K_GEOM::polyline, METH_VARARGS},
  {"spline", K_GEOM::spline, METH_VARARGS},
  {"nurbs", K_GEOM::nurbs, METH_VARARGS},
  {"getSharpestAngle", K_GEOM::getSharpestAngleForVertices, METH_VARARGS},
  {"getUV", K_GEOM::getUV, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "geom",
        NULL,
        -1,
        Pygeom
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_geom();
  PyMODINIT_FUNC PyInit_geom()
#else
  PyMODINIT_FUNC initgeom();
  PyMODINIT_FUNC initgeom()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("geom", Pygeom);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

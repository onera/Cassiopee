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
#include "cplot.h"
#include "Data.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef PycplotOSMesa [] =
{
  {"render", K_CPLOT::render, METH_VARARGS},
  {"delete", K_CPLOT::deletez, METH_VARARGS},
  {"add", K_CPLOT::add, METH_VARARGS},
  {"replace", K_CPLOT::replace, METH_VARARGS},
  {"pressKey", K_CPLOT::pressKey, METH_VARARGS},
  {"displayNew", K_CPLOT::displayNew, METH_VARARGS},
  {"displayAgain", K_CPLOT::displayAgain, METH_VARARGS},
  {"setFileName", K_CPLOT::setFileName, METH_VARARGS},
  {"getState", K_CPLOT::getState, METH_VARARGS},
  {"getSelectedZone", K_CPLOT::getSelectedZone, METH_VARARGS},
  {"getSelectedZones", K_CPLOT::getSelectedZones, METH_VARARGS},
  {"getSelectedStatus", K_CPLOT::getSelectedStatus, METH_VARARGS},
  {"getActiveZones", K_CPLOT::getActiveZones, METH_VARARGS},
  {"getActiveStatus", K_CPLOT::getActiveStatus, METH_VARARGS},
  {"getActivePoint", K_CPLOT::getActivePoint, METH_VARARGS},
  {"getActivePointIndex", K_CPLOT::getActivePointIndex, METH_VARARGS},
  {"getActivePointF", K_CPLOT::getActivePointF, METH_VARARGS},
  {"getMouseState", K_CPLOT::getMouseState, METH_VARARGS},
  {"getKeyboard", K_CPLOT::getKeyboard, METH_VARARGS},
  {"resetKeyboard", K_CPLOT::resetKeyboard, METH_VARARGS},
  {"setState", K_CPLOT::setState, METH_VARARGS},
  {"setMode", K_CPLOT::setMode, METH_VARARGS},
  {"setActivePoint", K_CPLOT::setActivePoint, METH_VARARGS},
  {"setShaderPath", K_CPLOT::setShaderPath, METH_VARARGS},
  {"setWindowTitle", K_CPLOT::setWindowTitle, METH_VARARGS},
  {"changeVariable", K_CPLOT::changeVariable, METH_VARARGS},
  {"changeStyle", K_CPLOT::changeStyle, METH_VARARGS},
  {"changeInfoDisplay", K_CPLOT::changeInfoDisplay, METH_VARARGS},
  {"changeBlanking", K_CPLOT::changeBlanking, METH_VARARGS},
  {"setDim", K_CPLOT::setDim, METH_VARARGS},
  {"setSelectedZones", K_CPLOT::setSelectedZones, METH_VARARGS},
  {"unselectAllZones", K_CPLOT::unselectAllZones, METH_VARARGS},
  {"setActiveZones", K_CPLOT::setActiveZones, METH_VARARGS},
  {"setZoneNames", K_CPLOT::setZoneNames, METH_VARARGS},
  {"lookFor", K_CPLOT::lookFor, METH_VARARGS},
  {"fitView", K_CPLOT::fitView, METH_VARARGS},
  {"isDisplayRunning", K_CPLOT::isDisplayRunning, METH_VARARGS},
  {"finalizeExport", K_CPLOT::finalizeExport, METH_VARARGS},
  {"hide", K_CPLOT::hide, METH_VARARGS},
  {"show", K_CPLOT::show, METH_VARARGS},
  {"display1D", K_CPLOT::display1D, METH_VARARGS},
  {"configure", K_CPLOT::configure, METH_VARARGS},
  {"panorama", K_CPLOT::panorama, METH_VARARGS},
  {"panoramaODS", K_CPLOT::panoramaODS, METH_VARARGS},
  {"blur", K_CPLOT::blur, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledefOSMesa = {
        PyModuleDef_HEAD_INIT,
        "cplotOSMesa",
        NULL,
        -1,
        PycplotOSMesa
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_cplotOSMesa();
  PyMODINIT_FUNC PyInit_cplotOSMesa()
#else
  PyMODINIT_FUNC initcplotOSMesa();
  PyMODINIT_FUNC initcplotOSMesa()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* m = PyModule_Create(&moduledefOSMesa);
#else
    PyObject* m = Py_InitModule("cplotOSMesa", PycplotOSMesa);
#endif
    PyModule_AddIntConstant(m, "useDirect", long(Data::Direct));
    PyModule_AddIntConstant(m, "useDL",     long(Data::DL));
    PyModule_AddIntConstant(m, "useVBO",    long(Data::VBO));
#if PY_MAJOR_VERSION >= 3
    return m;
#endif
  }
}

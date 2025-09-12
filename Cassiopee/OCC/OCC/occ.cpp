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
#include "occ.h"
#include "Standard_Version.hxx"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyocc [] =
{
  {"convertCAD2Arrays0", K_OCC::convertCAD2Arrays0, METH_VARARGS},
  {"convertCAD2Arrays1", K_OCC::convertCAD2Arrays1, METH_VARARGS},
  {"convertCAD2Arrays2", K_OCC::convertCAD2Arrays2, METH_VARARGS},

  {"readCAD", K_OCC::readCAD, METH_VARARGS},
  {"writeCAD", K_OCC::writeCAD, METH_VARARGS},
  {"createEmptyCAD", K_OCC::createEmptyCAD, METH_VARARGS},
  {"mergeCAD", K_OCC::mergeCAD, METH_VARARGS},
  {"freeHook", K_OCC::freeHook, METH_VARARGS},
  
  {"printOCAF", K_OCC::printOCAF, METH_VARARGS},
  {"printShapeOCAF", K_OCC::printShapeOCAF, METH_VARARGS},
  {"getFaceNameInOCAF", K_OCC::getFaceNameInOCAF, METH_VARARGS},
  {"getFaceNameInOCAF2", K_OCC::getFaceNameInOCAF2, METH_VARARGS},
  {"getEdgeNameInOCAF2", K_OCC::getEdgeNameInOCAF2, METH_VARARGS},
  
  {"bottle", K_OCC::bottle, METH_VARARGS},
  {"addSphere", K_OCC::addSphere, METH_VARARGS},
  {"addCylinder", K_OCC::addCylinder, METH_VARARGS},
  {"addBox", K_OCC::addBox, METH_VARARGS},
  {"addSquare", K_OCC::addSquare, METH_VARARGS},
  {"addLine", K_OCC::addLine, METH_VARARGS},
  {"addCircle", K_OCC::addCircle, METH_VARARGS},
  {"addSpline", K_OCC::addSpline, METH_VARARGS},
  {"addArc", K_OCC::addArc, METH_VARARGS},

  {"getNbFaces", K_OCC::getNbFaces, METH_VARARGS},
  {"getNbEdges", K_OCC::getNbEdges, METH_VARARGS},
  {"getFileAndFormat", K_OCC::getFileAndFormat, METH_VARARGS},
  
  {"meshGlobalEdges1", K_OCC::meshGlobalEdges1, METH_VARARGS},
  {"meshGlobalEdges2", K_OCC::meshGlobalEdges2, METH_VARARGS},
  {"meshGlobalEdges3", K_OCC::meshGlobalEdges3, METH_VARARGS},
  {"meshGlobalEdges4", K_OCC::meshGlobalEdges4, METH_VARARGS},
  {"meshEdgesByFace", K_OCC::meshEdgesByFace, METH_VARARGS},
  {"meshEdgesByFace2", K_OCC::meshEdgesByFace2, METH_VARARGS},
  {"meshEdgesByFace3", K_OCC::meshEdgesByFace3, METH_VARARGS},
  {"getEdgeNoByFace", K_OCC::getEdgeNoByFace, METH_VARARGS},
  {"identifyLoopsInEdges", K_OCC::identifyLoopsInEdges, METH_VARARGS},
  {"evalEdge", K_OCC::evalEdge, METH_VARARGS},
  {"evalFace", K_OCC::evalFace, METH_VARARGS},
  {"projectOnFaces", K_OCC::projectOnFaces, METH_VARARGS},
  {"projectOnEdges", K_OCC::projectOnEdges, METH_VARARGS},
  {"linkNodes2CAD", K_OCC::linkNodes2CAD, METH_VARARGS},
  {"updateFcadidFromNcadid", K_OCC::updateFcadidFromNcadid, METH_VARARGS},
  {"updateNcadidFromFcadid", K_OCC::updateNcadidFromFcadid, METH_VARARGS},
  {"getNodalParameters", K_OCC::getNodalParameters, METH_VARARGS},
  {"trimesh", K_OCC::trimesh, METH_VARARGS},

  {"meshOneEdge", K_OCC::meshOneEdge, METH_VARARGS},
  {"meshEdgesOfFace", K_OCC::meshEdgesOfFace, METH_VARARGS},

  {"analyseEdges", K_OCC::analyseEdges, METH_VARARGS},
  {"getFaceArea", K_OCC::getFaceArea, METH_VARARGS},
  {"getFaceOrientation", K_OCC::getFaceOrientation, METH_VARARGS},
  {"areEdgeIdentical", K_OCC::areEdgeIdentical, METH_VARARGS},

  {"splitFaces", K_OCC::splitFaces, METH_VARARGS},
  {"fixShape", K_OCC::fixShape, METH_VARARGS},
  {"trimFaces", K_OCC::trimFaces, METH_VARARGS},
  {"sewing", K_OCC::sewing, METH_VARARGS},
  {"removeFaces", K_OCC::removeFaces, METH_VARARGS},
  {"fillHole", K_OCC::fillHole, METH_VARARGS},
  {"addFillet", K_OCC::addFillet, METH_VARARGS},
  {"mergeFaces", K_OCC::mergeFaces, METH_VARARGS},
  {"loft", K_OCC::loft, METH_VARARGS},
  {"revolve", K_OCC::revolve, METH_VARARGS},

  {"translate", K_OCC::translate, METH_VARARGS},
  {"scale", K_OCC::scale, METH_VARARGS},
  {"rotate", K_OCC::rotate, METH_VARARGS},

  {"getOppData", K_OCC::getOppData, METH_VARARGS},
  {"identifyTags", K_OCC::identifyTags, METH_VARARGS},

  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "occ",
        NULL,
        -1,
        Pyocc
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_occ();
  PyMODINIT_FUNC PyInit_occ()
#else
  PyMODINIT_FUNC initocc();
  PyMODINIT_FUNC initocc()
#endif
  {
    import_array();
    printf("INFO: OCC version " OCC_VERSION_STRING "\n");

#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("occ", Pyocc);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

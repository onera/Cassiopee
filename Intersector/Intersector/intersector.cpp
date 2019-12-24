/*    
    Copyright 2013-2020 Onera.

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
#include "intersector.h"

int __activation__;

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyintersector [] =
{
  
  {"conformUnstr", K_INTERSECTOR::conformUnstr, METH_VARARGS},
  
  {"booleanIntersection", K_INTERSECTOR::booleanIntersection, METH_VARARGS},
  {"booleanUnion", K_INTERSECTOR::booleanUnion, METH_VARARGS},
  {"booleanUnionMZ", K_INTERSECTOR::booleanUnionMZ, METH_VARARGS},
  {"booleanMinus", K_INTERSECTOR::booleanMinus, METH_VARARGS},
  {"booleanIntersectionBorder", K_INTERSECTOR::booleanIntersectionBorder, METH_VARARGS},
  {"booleanModifiedSolid", K_INTERSECTOR::booleanModifiedSolid, METH_VARARGS},
  {"DiffSurf", K_INTERSECTOR::DiffSurf, METH_VARARGS},
  {"unify", K_INTERSECTOR::unify, METH_VARARGS},
  {"P1ConservativeChimeraCoeffs", K_INTERSECTOR::P1ConservativeChimeraCoeffs, METH_VARARGS},
  {"selfX", K_INTERSECTOR::selfX, METH_VARARGS},
  {"triangulateExteriorFaces", K_INTERSECTOR::triangulateExteriorFaces, METH_VARARGS},
  {"triangulateSpecifiedFaces", K_INTERSECTOR::triangulateSpecifiedFaces, METH_VARARGS},
  {"triangulateNFaces", K_INTERSECTOR::triangulateNFaces, METH_VARARGS},
  {"convexifyFaces", K_INTERSECTOR::convexifyFaces, METH_VARARGS},
  {"prepareCellsSplit", K_INTERSECTOR::prepareCellsSplit, METH_VARARGS},
  {"simplifyCells", K_INTERSECTOR::simplifyCells, METH_VARARGS},
  {"simplifySurf", K_INTERSECTOR::simplifySurf, METH_VARARGS},
  {"splitNonStarCells", K_INTERSECTOR::splitNonStarCells, METH_VARARGS},
  {"collapseUncomputableFaces", K_INTERSECTOR::collapseUncomputableFaces, METH_VARARGS},
  {"removeNonManifoldExternalCells", K_INTERSECTOR::removeNonManifoldExternalCells, METH_VARARGS},
  {"agglomerateSmallCells", K_INTERSECTOR::agglomerateSmallCells, METH_VARARGS},
  {"agglomerateNonStarCells", K_INTERSECTOR::agglomerateNonStarCells, METH_VARARGS},
  //{"agglomerateUncomputableCells", K_INTERSECTOR::agglomerateUncomputableCells, METH_VARARGS},
  {"agglomerateCellsWithSpecifiedFaces", K_INTERSECTOR::agglomerateCellsWithSpecifiedFaces, METH_VARARGS},
  {"adaptCells", K_INTERSECTOR::adaptCells, METH_VARARGS},
  {"adaptCellsNodal", K_INTERSECTOR::adaptCellsNodal, METH_VARARGS},
  {"adaptBox", K_INTERSECTOR::adaptBox, METH_VARARGS},
  {"createHMesh", K_INTERSECTOR::createHMesh, METH_VARARGS},
  {"deleteHMesh", K_INTERSECTOR::deleteHMesh, METH_VARARGS},
  {"conformizeHMesh", K_INTERSECTOR::conformizeHMesh, METH_VARARGS},
  {"closeOctalCells", K_INTERSECTOR::closeOctalCells, METH_VARARGS},
  {"extractUncomputables", K_INTERSECTOR::extractUncomputables, METH_VARARGS},
  {"extractPathologicalCells", K_INTERSECTOR::extractPathologicalCells, METH_VARARGS},
  {"extractOuterLayers", K_INTERSECTOR::extractOuterLayers, METH_VARARGS},
  {"extractNthCell", K_INTERSECTOR::extractNthCell, METH_VARARGS},
  {"extractNthFace", K_INTERSECTOR::extractNthFace, METH_VARARGS},
  {"removeNthCell", K_INTERSECTOR::removeNthCell, METH_VARARGS},

  {"getOverlappingFaces", K_INTERSECTOR::getOverlappingFaces, METH_VARARGS},
  {"getAnisoInnerFaces", K_INTERSECTOR::getAnisoInnerFaces, METH_VARARGS},

  {"statsUncomputableFaces", K_INTERSECTOR::statsUncomputableFaces, METH_VARARGS},
  {"statsSize", K_INTERSECTOR::statsSize, METH_VARARGS},
  
  {"computeAspectRatio", K_INTERSECTOR::computeAspectRatio, METH_VARARGS},
  {"centroids", K_INTERSECTOR::centroids, METH_VARARGS},
  
  {"diffMesh", K_INTERSECTOR::diffMesh, METH_VARARGS},

  { "checkCellsClosure", K_INTERSECTOR::checkCellsClosure, METH_VARARGS },
  { "checkForDegenCells", K_INTERSECTOR::checkForDegenCells, METH_VARARGS },
  { "edgeLengthExtrema", K_INTERSECTOR::edgeLengthExtrema, METH_VARARGS },
  { "removeBaffles", K_INTERSECTOR::removeBaffles, METH_VARARGS },
  { "convert2Polyhedron", K_INTERSECTOR::convert2Polyhedron, METH_VARARGS },
  { "oneZonePerCell", K_INTERSECTOR::oneZonePerCell, METH_VARARGS },
  
  { "extrudeBC", K_INTERSECTOR::extrudeBC, METH_VARARGS },
  { "extrudeSurf", K_INTERSECTOR::extrudeSurf, METH_VARARGS },
  { "extrudeRevolSurf", K_INTERSECTOR::extrudeRevolSurf, METH_VARARGS },

  { "reorientExternalFaces", K_INTERSECTOR::reorientExternalFaces, METH_VARARGS },
  { "reorientSpecifiedFaces", K_INTERSECTOR::reorientSpecifiedFaces, METH_VARARGS },
  { "reorientSurf", K_INTERSECTOR::reorientSurf, METH_VARARGS },

  { "convertNGON2DToNGON3D", K_INTERSECTOR::convertNGON2DToNGON3D, METH_VARARGS },
  { "convertBasic2NGONFaces", K_INTERSECTOR::convertBasic2NGONFaces, METH_VARARGS },
  { "oneph", K_INTERSECTOR::oneph, METH_VARARGS },

  /////////// syncronizing the tree ///////////
  { "updatePointLists", K_INTERSECTOR::updatePointLists, METH_VARARGS },
  /////////////////////////////////////////////
  { "merge", K_INTERSECTOR::merge, METH_VARARGS },
  { "concatenate", K_INTERSECTOR::concatenate, METH_VARARGS },

  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
struct module_state {
    PyObject *error;
};
static int myextension_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}
static int myextension_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "intersector",
        NULL,
        sizeof(struct module_state),
        Pyintersector,
        NULL,
        myextension_traverse,
        myextension_clear,
        NULL
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_intersector();
  PyMODINIT_FUNC PyInit_intersector()
#else
  PyMODINIT_FUNC initintersector();
  PyMODINIT_FUNC initintersector()
#endif
  {
    __activation__ = K_KCORE::activation("0");
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("intersector", Pyintersector);
#endif
    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

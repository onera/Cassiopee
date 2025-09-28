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
#include "transform.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef PyTransform[] =
{
  {"_cart2CylA", K_TRANSFORM::_cart2CylA, METH_VARARGS},
  {"_cart2CylZ", K_TRANSFORM::_cart2CylZ, METH_VARARGS},
  {"_cyl2CartA", K_TRANSFORM::_cyl2CartA, METH_VARARGS},
  {"_cyl2CartZ", K_TRANSFORM::_cyl2CartZ, METH_VARARGS},
  {"translate", K_TRANSFORM::translate, METH_VARARGS},
  //{"rotateA1", K_TRANSFORM::rotateA1, METH_VARARGS},
  //{"rotateA2", K_TRANSFORM::rotateA2, METH_VARARGS},
  //{"rotateA3", K_TRANSFORM::rotateA3, METH_VARARGS},
  {"_rotateA1", K_TRANSFORM::_rotateA1, METH_VARARGS},
  {"_rotateA2", K_TRANSFORM::_rotateA2, METH_VARARGS},
  {"_rotateA3", K_TRANSFORM::_rotateA3, METH_VARARGS},
  {"homothety", K_TRANSFORM::homothety, METH_VARARGS},
  {"contract", K_TRANSFORM::contract, METH_VARARGS},
  {"symetrize", K_TRANSFORM::symetrize, METH_VARARGS},
  {"perturbate", K_TRANSFORM::perturbate, METH_VARARGS},
  {"smooth", K_TRANSFORM::smooth, METH_VARARGS},
  {"_smoothField", K_TRANSFORM::_smoothField, METH_VARARGS},
  {"deform", K_TRANSFORM::deform, METH_VARARGS},
  {"deform2", K_TRANSFORM::deform2, METH_VARARGS},
  {"deformPoint", K_TRANSFORM::deformPoint, METH_VARARGS},
  {"deformMeshStruct", K_TRANSFORM::deformMeshStruct, METH_VARARGS},
  {"computeDeformationVector", K_TRANSFORM::computeDeformationVector, METH_VARARGS},
  {"_freeForm", K_TRANSFORM::_freeForm, METH_VARARGS},
  {"projectAllDirs", K_TRANSFORM::projectAllDirs, METH_VARARGS},
  {"projectDir", K_TRANSFORM::projectDir, METH_VARARGS},
  {"projectOrtho", K_TRANSFORM::projectOrtho, METH_VARARGS},
  {"projectOrthoSmooth", K_TRANSFORM::projectOrthoSmooth, METH_VARARGS},
  {"projectRay", K_TRANSFORM::projectRay, METH_VARARGS},
  {"projectSmoothDir", K_TRANSFORM::projectSmoothDir, METH_VARARGS},
  {"_alignVectorFieldWithRadialCylindricProjection", K_TRANSFORM::_alignVectorFieldWithRadialCylindricProjection, METH_VARARGS},
  {"join", K_TRANSFORM::join, METH_VARARGS},
  {"joinBoth", K_TRANSFORM::joinBoth, METH_VARARGS},
  {"joinAll", K_TRANSFORM::joinAll, METH_VARARGS},
  {"joinAllBoth", K_TRANSFORM::joinAllBoth, METH_VARARGS},
  {"patch", K_TRANSFORM::patch, METH_VARARGS},
  {"patch2", K_TRANSFORM::patch2, METH_VARARGS},
  {"subzoneStruct", K_TRANSFORM::subzoneStruct, METH_VARARGS},
  {"subzoneStructInt", K_TRANSFORM::subzoneStructInt, METH_VARARGS},
  {"subzoneStructIntBoth", K_TRANSFORM::subzoneStructIntBoth, METH_VARARGS},
  {"subzoneUnstruct", K_TRANSFORM::subzoneUnstruct, METH_VARARGS},
  {"subzoneUnstructBoth",K_TRANSFORM::subzoneUnstructBoth, METH_VARARGS},
  {"subzoneElements",K_TRANSFORM::subzoneElements, METH_VARARGS},
  {"subzoneElementsBoth",K_TRANSFORM::subzoneElementsBoth, METH_VARARGS},
  {"subzoneFaces",K_TRANSFORM::subzoneFaces, METH_VARARGS},
  {"subzoneFacesBoth",K_TRANSFORM::subzoneFacesBoth, METH_VARARGS},
  {"oneovern", K_TRANSFORM::oneovern, METH_VARARGS},
  {"reorder", K_TRANSFORM::reorder, METH_VARARGS},
  {"reorderAll", K_TRANSFORM::reorderAll, METH_VARARGS},
  {"reorderAllUnstr", K_TRANSFORM::reorderAllUnstr, METH_VARARGS},
  {"addkplane", K_TRANSFORM::addkplane, METH_VARARGS},
  {"addkplaneCenters", K_TRANSFORM::addkplaneCenters, METH_VARARGS},
  {"splitCurvatureAngle", K_TRANSFORM::splitCurvatureAngle, METH_VARARGS},
  {"splitCurvatureRadius", K_TRANSFORM::splitCurvatureRadius, METH_VARARGS},
  {"splitConnexity", K_TRANSFORM::splitConnexity, METH_VARARGS},
  {"splitSharpEdges", K_TRANSFORM::splitSharpEdges, METH_VARARGS},
  {"splitSharpEdgesList", K_TRANSFORM::splitSharpEdgesList, METH_VARARGS},
  {"splitBAR", K_TRANSFORM::splitBAR, METH_VARARGS},
  {"splitTBranches",K_TRANSFORM::splitTBranches, METH_VARARGS},
  {"splitTRI", K_TRANSFORM::splitTRI, METH_VARARGS},
  {"splitManifold", K_TRANSFORM::splitManifold, METH_VARARGS},
  {"collapse", K_TRANSFORM::collapse, METH_VARARGS},
  {"mergeCart", K_TRANSFORM::mergeCartGrids, METH_VARARGS},
  {"merge", K_TRANSFORM::mergeStructGrids, METH_VARARGS},
  {"breakElements", K_TRANSFORM::breakElements, METH_VARARGS},
  {"splitNGon", K_TRANSFORM::splitNGon, METH_VARARGS},
  {"splitNGon2", K_TRANSFORM::splitNGon2, METH_VARARGS},
  {"splitElement", K_TRANSFORM::splitElement, METH_VARARGS},
  {"dualNGon", K_TRANSFORM::dualNGon, METH_VARARGS},
  {"flipEdges", K_TRANSFORM::flipEdges, METH_VARARGS},
  {"contractEdges", K_TRANSFORM::contractEdges, METH_VARARGS},
  {"checkTriMesh", K_TRANSFORM::checkTriMesh, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "transform",
        NULL,
        -1,
        PyTransform
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_transform();
  PyMODINIT_FUNC PyInit_transform()
#else
  PyMODINIT_FUNC inittransform();
  PyMODINIT_FUNC inittransform()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("transform", PyTransform);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

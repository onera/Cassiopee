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
#include "generator.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pygenerator [] =
{
  {"cart", K_GENERATOR::cartStruct, METH_VARARGS},
  {"cartr1", K_GENERATOR::cartr1, METH_VARARGS},
  {"cartr2", K_GENERATOR::cartr2, METH_VARARGS},
  {"cartHexa", K_GENERATOR::cartHexa, METH_VARARGS},
  {"cartTetra", K_GENERATOR::cartTetra, METH_VARARGS},
  {"cartPenta", K_GENERATOR::cartPenta, METH_VARARGS},
  {"cartPyra", K_GENERATOR::cartPyra, METH_VARARGS},
  {"cartNGon", K_GENERATOR::cartNGon, METH_VARARGS},
  {"cylinder", K_GENERATOR::cylinderMesh, METH_VARARGS},
  {"cylinder2", K_GENERATOR::cylinderMesh2, METH_VARARGS},
  {"cylinder3", K_GENERATOR::cylinderMesh3, METH_VARARGS},
  {"delaunay", K_GENERATOR::delaunay, METH_VARARGS},
  {"grow", K_GENERATOR::growMesh, METH_VARARGS},
  {"stack", K_GENERATOR::stackMesh, METH_VARARGS},
  {"TFI", K_GENERATOR::TFIMesh, METH_VARARGS},
  {"map", K_GENERATOR::mapMesh, METH_VARARGS},
  {"TTM", K_GENERATOR::TTMMesh, METH_VARARGS},
  {"getVolumeMap", K_GENERATOR::getVolumeMapOfMesh, METH_VARARGS},
  {"getCellCenters", K_GENERATOR::getCellCenters, METH_VARARGS},
  {"getFaceCentersAndAreas", K_GENERATOR::getFaceCentersAndAreas, METH_VARARGS},
  {"getOrthogonalityMap", K_GENERATOR::getOrthogonalityMap, METH_VARARGS},
  {"getRegularityMap", K_GENERATOR::getRegularityMap, METH_VARARGS},
  {"getAngleRegularityMap", K_GENERATOR::getAngleRegularityMap, METH_VARARGS},
  {"getNormalMap", K_GENERATOR::getNormalMapOfMesh, METH_VARARGS},
  {"getCircumCircleMap", K_GENERATOR::getCircumCircleMap, METH_VARARGS},
  {"getInCircleMap", K_GENERATOR::getInCircleMap, METH_VARARGS},
  {"closeMesh", K_GENERATOR::closeMesh, METH_VARARGS},
  {"closeMeshLegacy", K_GENERATOR::closeMeshLegacy, METH_VARARGS},
  {"closeBorders", K_GENERATOR::closeBorders, METH_VARARGS},
  {"pointedHat", K_GENERATOR::pointedHat, METH_VARARGS},
  {"stitchedHat", K_GENERATOR::stitchedHat, METH_VARARGS},
  {"hyper2D", K_GENERATOR::hyper2DMesh, METH_VARARGS},
  {"hyper2D2", K_GENERATOR::hyper2D2Mesh, METH_VARARGS},
  {"hyper2D3", K_GENERATOR::hyper2D3Mesh, METH_VARARGS},
  {"hyper2D4", K_GENERATOR::hyper2D4Mesh, METH_VARARGS},
  {"enforceLine", K_GENERATOR::enforceLineMesh, METH_VARARGS},
  {"enforce", K_GENERATOR::enforceMesh, METH_VARARGS},
  {"enforceY", K_GENERATOR::enforceYMesh, METH_VARARGS},
  {"enforceMoinsY", K_GENERATOR::enforceMoinsYMesh, METH_VARARGS},
  {"enforcePlusY", K_GENERATOR::enforcePlusYMesh, METH_VARARGS},
  {"enforceX", K_GENERATOR::enforceXMesh, METH_VARARGS},
  {"enforceMoinsX", K_GENERATOR::enforceMoinsXMesh, METH_VARARGS},
  {"enforcePlusX", K_GENERATOR::enforcePlusXMesh, METH_VARARGS},
  {"enforcePoint", K_GENERATOR::enforcePoint, METH_VARARGS},
  {"enforceCurvature", K_GENERATOR::enforceCurvature, METH_VARARGS},
  {"addPointInDistribution", K_GENERATOR::addPointInDistribution, METH_VARARGS},
  {"check", K_GENERATOR::checkMesh, METH_VARARGS},
  {"checkPointInCEBB", K_GENERATOR::checkPointInCEBBOfMesh, METH_VARARGS},
  {"bboxOfCells", K_GENERATOR::getBBOfCells, METH_VARARGS},
  {"obbox", K_GENERATOR::obbox, METH_VARARGS},
  {"barycenter", K_GENERATOR::barycenter, METH_VARARGS},
  {"CEBBIntersection", K_GENERATOR::getCEBBIntersectionOfArrays, METH_VARARGS},
  {"bboxIntersection", K_GENERATOR::bboxIntersection, METH_VARARGS},
  {"_bboxIntersectionZ", K_GENERATOR::_bboxIntersectionZ, METH_VARARGS},  
  {"obboxIntersection", K_GENERATOR::obboxIntersection, METH_VARARGS},
  {"_obboxIntersectionZ", K_GENERATOR::_obboxIntersectionZ, METH_VARARGS},
  {"crossIntersection", K_GENERATOR::crossIntersection, METH_VARARGS},
  {"_crossIntersectionZ", K_GENERATOR::_crossIntersectionZ, METH_VARARGS},
  {"getCellPlanarity", K_GENERATOR::computeCellPlanarity, METH_VARARGS},
  {"checkDelaunay", K_GENERATOR::checkDelaunay, METH_VARARGS},
  {"selectInsideElts", K_GENERATOR::selectInsideElts, METH_VARARGS},
  {"densify", K_GENERATOR::densifyMesh, METH_VARARGS},
  {"T3mesher2D", K_GENERATOR::T3mesher2D, METH_VARARGS},
  {"fittingPlaster", K_GENERATOR::fittingPlaster, METH_VARARGS},
  {"gapfixer", K_GENERATOR::gapfixer, METH_VARARGS},
  {"gapsmanager", K_GENERATOR::gapsmanager, METH_VARARGS},
  {"front2Hexa", K_GENERATOR::front2Hexa, METH_VARARGS},
  {"front2Struct", K_GENERATOR::front2Struct, METH_VARARGS},
  {"snapFront", K_GENERATOR::snapFront, METH_VARARGS},
  {"snapSharpEdges", K_GENERATOR::snapSharpEdges, METH_VARARGS},
  {"fillWithStruct", K_GENERATOR::fillWithStruct, METH_VARARGS},
  {"octree", K_GENERATOR::octree, METH_VARARGS},
  {"balanceOctree", K_GENERATOR::balanceOctree, METH_VARARGS},
  {"octree3", K_GENERATOR::octree3, METH_VARARGS},
  {"octree2AMR", K_GENERATOR::octree2AMR, METH_VARARGS},
  {"octree2Struct", K_GENERATOR::octree2Struct, METH_VARARGS},
  {"extendCartGrids",K_GENERATOR::extendCartGrids, METH_VARARGS},
  {"adaptOctree", K_GENERATOR::adaptOctree, METH_VARARGS},
  {"adaptOctree3", K_GENERATOR::adaptOctree3, METH_VARARGS}, 
  {"conformOctree3", K_GENERATOR::conformOctree3, METH_VARARGS},
  {"modifyIndicToExpandLayer", K_GENERATOR::modifyIndicToExpandLayer, METH_VARARGS},
  {"straightenVector", K_GENERATOR::straightenVector, METH_VARARGS},
  {"computeEta", K_GENERATOR::computeEta, METH_VARARGS},
  {"getLocalStepFactor", K_GENERATOR::getLocalStepFactor, METH_VARARGS},
  {"getLocalStepFactor2", K_GENERATOR::getLocalStepFactor2, METH_VARARGS},
  {"getEdgeRatio", K_GENERATOR::getEdgeRatio, METH_VARARGS},
  {"getMaxLength", K_GENERATOR::getMaxLength, METH_VARARGS},
  {"getTriQualityMap", K_GENERATOR::getTriQualityMap, METH_VARARGS},
  {"netgen1", K_GENERATOR::netgen1, METH_VARARGS},
  {"netgen2", K_GENERATOR::netgen2, METH_VARARGS},
  {"tetgen", K_GENERATOR::tetgen, METH_VARARGS},
  {"mmgs", K_GENERATOR::mmgs, METH_VARARGS},
  {"quad2Pyra", K_GENERATOR::quad2Pyra, METH_VARARGS},
  {"blankSelf", K_GENERATOR::blankSelf, METH_VARARGS},
  {"blankFirst", K_GENERATOR::blankFirst, METH_VARARGS},
  {"blankExt", K_GENERATOR::blankExt, METH_VARARGS},
  {"blankPrev", K_GENERATOR::blankPrev, METH_VARARGS},
  {"extrapWithCellN", K_GENERATOR::extrapWithCellN, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "generator",
        NULL,
        -1,
        Pygenerator
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_generator();
  PyMODINIT_FUNC PyInit_generator()
#else
  PyMODINIT_FUNC initgenerator();
  PyMODINIT_FUNC initgenerator()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("generator", Pygenerator);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

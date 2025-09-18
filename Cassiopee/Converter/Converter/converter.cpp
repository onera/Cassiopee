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
#include "converter.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyconverter [] =
{
  {"copy", K_CONVERTER::copy, METH_VARARGS},
  {"setPartialFields", K_CONVERTER::setPartialFields, METH_VARARGS},
  {"_setPartialFields", K_CONVERTER::_setPartialFields, METH_VARARGS},
  {"setPartialFieldsPT", K_CONVERTER::setPartialFieldsPT, METH_VARARGS},
  {"updatePartialFields", K_CONVERTER::updatePartialFields, METH_VARARGS},
  {"updatePartialFieldsPT", K_CONVERTER::updatePartialFieldsPT, METH_VARARGS},
  {"_updatePartialFields", K_CONVERTER::_updatePartialFields, METH_VARARGS},
  {"filterPartialFields", K_CONVERTER::filterPartialFields, METH_VARARGS},
  {"_setPartialFieldsAverage", K_CONVERTER::_setPartialFieldsAverage, METH_VARARGS},
  {"extractVars", K_CONVERTER::extractVars, METH_VARARGS},
  {"addVars", K_CONVERTER::addVars, METH_VARARGS},
  {"addVar", K_CONVERTER::addVar, METH_VARARGS},
  {"randomizeVar", K_CONVERTER::randomizeVar, METH_VARARGS},
  {"convertFile2Arrays", K_CONVERTER::convertFile2Arrays, METH_VARARGS},
  {"convertArrays2File", K_CONVERTER::convertArrays2File, METH_VARARGS},
  {"initVars", K_CONVERTER::initVars, METH_VARARGS},
  {"diffArrays", K_CONVERTER::diffArrays, METH_VARARGS},
  {"getArgMin", K_CONVERTER::getArgMin, METH_VARARGS},
  {"getArgMax", K_CONVERTER::getArgMax, METH_VARARGS},
  {"getMeanValue", K_CONVERTER::getMeanValue, METH_VARARGS},
  {"getMeanRangeValue", K_CONVERTER::getMeanRangeValue, METH_VARARGS},
  {"normL0", K_CONVERTER::normL0, METH_VARARGS},
  {"normL2", K_CONVERTER::normL2, METH_VARARGS},
  {"normalize", K_CONVERTER::normalize, METH_VARARGS},
  {"magnitude", K_CONVERTER::magnitude, METH_VARARGS},
  {"isFinite", K_CONVERTER::isFinite, METH_VARARGS},
  {"setNANValuesAt", K_CONVERTER::setNANValuesAt, METH_VARARGS},
  {"convertBAR2Struct", K_CONVERTER::convertBAR2Struct, METH_VARARGS},
  {"convertStruct2Tetra", K_CONVERTER::convertStruct2Tetra, METH_VARARGS},
  {"convertStruct2TetraBary", K_CONVERTER::convertStruct2TetraBary, METH_VARARGS},
  {"convertStruct2TetraBaryBoth", K_CONVERTER::convertStruct2TetraBaryBoth, METH_VARARGS},
  {"convertStruct2Hexa", K_CONVERTER::convertStruct2Hexa, METH_VARARGS},
  {"convertStruct2NGon", K_CONVERTER::convertStruct2NGon, METH_VARARGS},
  {"convertHexa2Struct", K_CONVERTER::convertHexa2Struct, METH_VARARGS},
  {"convertUnstruct2NGon", K_CONVERTER::convertUnstruct2NGon, METH_VARARGS},
  {"convertUnstruct2Hexa", K_CONVERTER::convertUnstruct2Hexa, METH_VARARGS},
  {"convertHexa2Tetra", K_CONVERTER::convertHexa2Tetra, METH_VARARGS},
  {"convertPenta2Tetra", K_CONVERTER::convertPenta2Tetra, METH_VARARGS},
  {"convertPyra2Tetra", K_CONVERTER::convertPyra2Tetra, METH_VARARGS},
  {"convertArray2TetraBary", K_CONVERTER::convertArray2TetraBary, METH_VARARGS},
  {"convertArray2TetraBaryBoth", K_CONVERTER::convertArray2TetraBaryBoth, METH_VARARGS},
  {"convertNGon2TetraBary", K_CONVERTER::convertNGon2TetraBary, METH_VARARGS},
  {"convertNGon2TetraBaryBoth", K_CONVERTER::convertNGon2TetraBaryBoth, METH_VARARGS},
  {"convertHO2LO", K_CONVERTER::convertHO2LO, METH_VARARGS},
  {"convertLO2HO", K_CONVERTER::convertLO2HO, METH_VARARGS},
  {"convertTri2Quad", K_CONVERTER::convertTri2Quad, METH_VARARGS},
  {"convertQuad2Tri", K_CONVERTER::convertQuad2Tri, METH_VARARGS},
  {"convertMix2BE", K_CONVERTER::convertMix2BE, METH_VARARGS},
  {"convertArray2Node", K_CONVERTER::convertArray2Node, METH_VARARGS},
  {"convertStrand2Penta", K_CONVERTER::convertStrand2Penta, METH_VARARGS},
  {"convertPenta2Strand", K_CONVERTER::convertPenta2Strand, METH_VARARGS},
  {"node2Center", K_CONVERTER::node2Center, METH_VARARGS},
  {"node2Center_OLD", K_CONVERTER::node2Center_OLD, METH_VARARGS},
  {"center2Node", K_CONVERTER::center2Node, METH_VARARGS},
  {"center2Node_OLD", K_CONVERTER::center2Node_OLD, METH_VARARGS},
  {"node2ExtCenter", K_CONVERTER::node2ExtCenter, METH_VARARGS},
  {"extCenter2Node", K_CONVERTER::extCenter2Node, METH_VARARGS},
  {"center2ExtCenter", K_CONVERTER::center2ExtCenter, METH_VARARGS},
  {"convertFile2PyTree", K_CONVERTER::convertFile2PyTree, METH_VARARGS},
  {"convertFile2PartialPyTree", K_CONVERTER::convertFile2PartialPyTree, METH_VARARGS},
  {"convertPyTree2File", K_CONVERTER::convertPyTree2File, METH_VARARGS},
  {"convertFile2PyTreeFromPath", K_CONVERTER::convertFile2PyTreeFromPath, METH_VARARGS},
  {"convertPyTree2FilePartial", K_CONVERTER::convertPyTree2FilePartial, METH_VARARGS},
  {"readPyTreeFromPaths", K_CONVERTER::readPyTreeFromPaths, METH_VARARGS},
  {"writePyTreePaths", K_CONVERTER::writePyTreePaths, METH_VARARGS},
  {"deletePyTreePaths", K_CONVERTER::deletePyTreePaths, METH_VARARGS},

  {"convertFile2PyTreeTau", K_CONVERTER::convertFile2PyTreeTau, METH_VARARGS},
  {"convertPyTree2FileTau", K_CONVERTER::convertPyTree2FileTau, METH_VARARGS},
  {"convertFile2PyTreeFsdm", K_CONVERTER::convertFile2PyTreeFsdm, METH_VARARGS},
  {"convertPyTree2FileFsdm", K_CONVERTER::convertPyTree2FileFsdm, METH_VARARGS},
  {"convertPyTree2FFD", K_CONVERTER::convertPyTree2FFD, METH_VARARGS},
  
  {"addGhostCellsNGonNodes", K_CONVERTER::addGhostCellsNGonNodes, METH_VARARGS},
  {"addGhostCellsNGonCenters", K_CONVERTER::addGhostCellsNGonCenters, METH_VARARGS},
  {"addGhostCellsNGonBoth", K_CONVERTER::addGhostCellsNGonBoth, METH_VARARGS},
  {"rmGhostCellsNGonNodes", K_CONVERTER::rmGhostCellsNGonNodes, METH_VARARGS},
  {"rmGhostCellsNGonCenters", K_CONVERTER::rmGhostCellsNGonCenters, METH_VARARGS},
  {"rmGhostCellsNGonBoth", K_CONVERTER::rmGhostCellsNGonBoth, METH_VARARGS},
  {"cpyGhost2Real", K_CONVERTER::cpyGhost2Real, METH_VARARGS},
  {"cpyReal2Ghost", K_CONVERTER::cpyReal2Ghost, METH_VARARGS},
  {"_setBCDataInGhostCellsStruct", K_CONVERTER::setBCDataInGhostCellsStruct, METH_VARARGS},
  {"extrapInterior2BCFaceStruct", K_CONVERTER::extrapInterior2BCFaceStruct, METH_VARARGS},
  {"nullifyVectorAtBCFaceStruct", K_CONVERTER::nullifyVectorAtBCFaceStruct, METH_VARARGS},
  {"cpyConnectA2ConnectP", K_CONVERTER::cpyConnectA2ConnectP, METH_VARARGS},
  {"cpyConnectP2ConnectA", K_CONVERTER::cpyConnectP2ConnectA, METH_VARARGS},
  {"cpyConnectP2ConnectA2", K_CONVERTER::cpyConnectP2ConnectA2, METH_VARARGS},
  {"cpyValueByField", K_CONVERTER::cpyValueByField, METH_VARARGS},
  {"detectEmptyBC", K_CONVERTER::detectEmptyBC, METH_VARARGS},
  {"tagDefinedBC", K_CONVERTER::tagDefinedBC, METH_VARARGS},
  {"fillJoin", K_CONVERTER::fillJoin, METH_VARARGS},
  {"fillJoinNMNodes", K_CONVERTER::fillJoinNMNodes, METH_VARARGS},
  {"fillJoinNMCenters", K_CONVERTER::fillJoinNMCenters, METH_VARARGS},
  {"fillCornerGhostCells", K_CONVERTER::fillCornerGhostCells, METH_VARARGS},
  {"fillCornerGhostCells2", K_CONVERTER::fillCornerGhostCells2, METH_VARARGS},
  {"getJoinBorderIndices", K_CONVERTER::getJoinBorderIndices, METH_VARARGS},
  {"getJoinDonorIndices", K_CONVERTER::getJoinDonorIndices, METH_VARARGS},
  {"conformizeNGon", K_CONVERTER::conformizeNGon, METH_VARARGS},
  {"registerFaces", K_CONVERTER::registerFaces, METH_VARARGS},
  {"registerCells", K_CONVERTER::registerCells, METH_VARARGS},
  {"registerNodes", K_CONVERTER::registerNodes, METH_VARARGS},
  {"registerElements", K_CONVERTER::registerElements, METH_VARARGS},
  {"registerAllFaces", K_CONVERTER::registerAllFaces, METH_VARARGS},
  {"registerAllNodes", K_CONVERTER::registerAllNodes, METH_VARARGS},
  {"registerAllElements", K_CONVERTER::registerAllElements, METH_VARARGS},
  {"freeHook", K_CONVERTER::freeHook, METH_VARARGS},
  {"identifyElements", K_CONVERTER::identifyElements, METH_VARARGS},
  {"identifyFaces", K_CONVERTER::identifyFaces, METH_VARARGS},
  {"identifyNodes", K_CONVERTER::identifyNodes, METH_VARARGS},
  {"identifySolutions", K_CONVERTER::identifySolutions, METH_VARARGS},
  {"nearestElements", K_CONVERTER::nearestElements, METH_VARARGS},
  {"nearestFaces", K_CONVERTER::nearestFaces, METH_VARARGS},
  {"nearestNodes", K_CONVERTER::nearestNodes, METH_VARARGS},
  {"createGlobalIndex", K_CONVERTER::createGlobalIndex, METH_VARARGS},
  {"recoverGlobalIndex", K_CONVERTER::recoverGlobalIndex, METH_VARARGS},
  {"adaptPE2NFace", K_CONVERTER::adaptPE2NFace, METH_VARARGS},
  {"adaptNFace2PE", K_CONVERTER::adaptNFace2PE, METH_VARARGS},
  {"adaptNGon2Index", K_CONVERTER::adaptNGon2Index, METH_VARARGS},
  {"adaptNFace2Index", K_CONVERTER::adaptNFace2Index, METH_VARARGS},
  {"adaptBCFace2BCC", K_CONVERTER::adaptBCFace2BCC, METH_VARARGS},
  {"adaptBCFacePL2VertexPL", K_CONVERTER::adaptBCFacePL2VertexPL, METH_VARARGS},
  {"adaptBCVertexPL2FacePL", K_CONVERTER::adaptBCVertexPL2FacePL, METH_VARARGS},
  {"adaptNGon42NGon3", K_CONVERTER::adaptNGon42NGon3, METH_VARARGS},
  {"adaptNGon32NGon4", K_CONVERTER::adaptNGon32NGon4, METH_VARARGS},
  {"adaptShiftedPE2PE", K_CONVERTER::adaptShiftedPE2PE, METH_VARARGS},
  {"signNGonFaces", K_CONVERTER::signNGonFaces, METH_VARARGS},
  {"unsignNGonFaces", K_CONVERTER::unsignNGonFaces, METH_VARARGS},
  {"sliceNGonFaces", K_CONVERTER::sliceNGonFaces, METH_VARARGS},
  {"makeParentElements", K_CONVERTER::makeParentElements, METH_VARARGS},
  {"adaptSurfaceNGon", K_CONVERTER::adaptSurfaceNGon, METH_VARARGS},
  {"adapt2FastP", K_CONVERTER::adapt2FastP, METH_VARARGS},
  {"createElsaHybrid", K_CONVERTER::createElsaHybrid, METH_VARARGS},
  {"diffIndex", K_CONVERTER::diffIndex, METH_VARARGS},
  {"pointList2Ranges", K_CONVERTER::pointList2Ranges, METH_VARARGS},
  {"pointList2SPL", K_CONVERTER::pointList2SPL, METH_VARARGS},
  {"range2PointList", K_CONVERTER::range2PointList, METH_VARARGS},
  {"PR2VL", K_CONVERTER::PR2VL, METH_VARARGS},
  {"addGhostCellsNG", K_CONVERTER::addGhostCellsNG, METH_VARARGS},
  {"extractBCMatchStruct", K_CONVERTER::extractBCMatchStruct, METH_VARARGS},
  {"extractBCMatchNG", K_CONVERTER::extractBCMatchNG, METH_VARARGS},
  {"buildBCMatchFieldStruct", K_CONVERTER::buildBCMatchFieldStruct, METH_VARARGS},
  {"buildBCMatchFieldNG", K_CONVERTER::buildBCMatchFieldNG, METH_VARARGS},
  {"extractBCFields", K_CONVERTER::extractBCFields, METH_VARARGS},
  {"extractFields", K_CONVERTER::extractFields, METH_VARARGS},
  {"iSend", K_CONVERTER::iSend, METH_VARARGS},
  {"recv", K_CONVERTER::recv, METH_VARARGS},
  {"waitAll", K_CONVERTER::waitAll, METH_VARARGS},
  {"createBBTree", K_CONVERTER::createBBTree, METH_VARARGS},
  {"intersect", K_CONVERTER::intersect, METH_VARARGS},
  {"intersect2", K_CONVERTER::intersect2, METH_VARARGS},
  {"deleteBBTree", K_CONVERTER::deleteBBTree, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "converter",
        NULL,
        -1,
        Pyconverter
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_converter();
  PyMODINIT_FUNC PyInit_converter()
#else
  PyMODINIT_FUNC initconverter();
  PyMODINIT_FUNC initconverter()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("converter", Pyconverter);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}


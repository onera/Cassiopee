/*    
    Copyright 2013-2019 Onera.

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
#include "connector.h"

int __activation__;

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyconnector [] =
{
  {"getIBMPtsBasic", K_CONNECTOR::getIBMPtsBasic, METH_VARARGS},
  {"getIBMPtsWithFront", K_CONNECTOR::getIBMPtsWithFront, METH_VARARGS},
  {"getIBMPtsWithoutFront", K_CONNECTOR::getIBMPtsWithoutFront, METH_VARARGS},
  {"optimizeOverlap", K_CONNECTOR::optimizeOverlap, METH_VARARGS},
  {"maximizeBlankedCells", K_CONNECTOR::maximizeBlankedCells, METH_VARARGS},
  {"blankCells", K_CONNECTOR::blankCells, METH_VARARGS},
  {"_blankCells", K_CONNECTOR::_blankCells, METH_VARARGS},
  {"blankCellsTetra", K_CONNECTOR::blankCellsTetra, METH_VARARGS},
  {"createTetraMask", K_CONNECTOR::createTetraMask, METH_VARARGS},
  {"deleteTetraMask", K_CONNECTOR::deleteTetraMask, METH_VARARGS},
  {"createTriMask", K_CONNECTOR::createTriMask, METH_VARARGS},
  {"deleteTriMask", K_CONNECTOR::deleteTriMask, METH_VARARGS},
  {"maskXRay", K_CONNECTOR::maskXRay, METH_VARARGS},
  {"getIntersectingDomainsAABB", K_CONNECTOR::getIntersectingDomainsAABB, METH_VARARGS},
  {"setDoublyDefinedBC", K_CONNECTOR::setDoublyDefinedBC, METH_VARARGS},
  {"getOversetHolesInterpCellCenters", K_CONNECTOR::getOversetHolesInterpCellCenters, METH_VARARGS},
  {"getOversetHolesInterpNodes", K_CONNECTOR::getOversetHolesInterpNodes, METH_VARARGS},
  {"_getOversetHolesInterpCellCenters", K_CONNECTOR::_getOversetHolesInterpCellCenters, METH_VARARGS},
  {"_getOversetHolesInterpNodes", K_CONNECTOR::_getOversetHolesInterpNodes, METH_VARARGS},
  {"getEXPoints", K_CONNECTOR::getEXPoints, METH_VARARGS},
  {"getInterpolatedPoints", K_CONNECTOR::getInterpolatedPoints, METH_VARARGS},
  {"getInterpolatedPointsZ", K_CONNECTOR::getInterpolatedPointsZ, METH_VARARGS},
  {"setInterpolations", K_CONNECTOR::setInterpolations, METH_VARARGS},
  {"setInterpData", K_CONNECTOR::setInterpData, METH_VARARGS},
  {"setInterpDataDW", K_CONNECTOR::setInterpDataDW, METH_VARARGS},
  {"setInterpDataForGC", K_CONNECTOR::setInterpDataForGC, METH_VARARGS},
  {"setInterpDataLS", K_CONNECTOR::setInterpDataLS, METH_VARARGS},
  {"setInterpDataCons", K_CONNECTOR::setInterpDataCons, METH_VARARGS},
  {"setInterpTransfers", K_CONNECTOR::setInterpTransfers, METH_VARARGS},
  {"initNuma", K_CONNECTOR::initNuma, METH_VARARGS},
  {"_setInterpTransfers", K_CONNECTOR::_setInterpTransfers, METH_VARARGS},
  {"__setInterpTransfers", K_CONNECTOR::__setInterpTransfers, METH_VARARGS},
  {"___setInterpTransfers", K_CONNECTOR::___setInterpTransfers, METH_VARARGS},
  {"setInterpTransfersD", K_CONNECTOR::setInterpTransfersD, METH_VARARGS},
  {"_setInterpTransfersD", K_CONNECTOR::_setInterpTransfersD, METH_VARARGS},
  {"__setInterpTransfersD", K_CONNECTOR::__setInterpTransfersD, METH_VARARGS},
  {"writeCoefs", K_CONNECTOR::writeCoefs, METH_VARARGS},
  {"chimeraTransfer", K_CONNECTOR::chimeraTransfer, METH_VARARGS},
  {"transferFields", K_CONNECTOR::transferFields, METH_VARARGS},
  {"changeWall", K_CONNECTOR::changeWall, METH_VARARGS},
  {"changeWallEX", K_CONNECTOR::changeWallEX, METH_VARARGS},
  {"blankIntersectingCells", K_CONNECTOR::blankIntersectingCells, METH_VARARGS},
  {"cellN2OversetHolesStruct", K_CONNECTOR::cellN2OversetHolesStruct, METH_VARARGS},
  {"cellN2OversetHolesUnStruct", K_CONNECTOR::cellN2OversetHolesUnStruct, METH_VARARGS},
  {"identifyMatching", K_CONNECTOR::identifyMatching, METH_VARARGS},
  {"identifyMatchingP", K_CONNECTOR::identifyMatchingP, METH_VARARGS},
  {"identifyMatchingNM", K_CONNECTOR::identifyMatchingNM, METH_VARARGS},
  {"identifyDegenerated", K_CONNECTOR::identifyDegenerated, METH_VARARGS},
  {"gatherMatching", K_CONNECTOR::gatherMatching, METH_VARARGS},
  {"gatherMatchingNM", K_CONNECTOR::gatherMatchingNM, METH_VARARGS},
  {"gatherMatchingNGon", K_CONNECTOR::gatherMatchingNGon, METH_VARARGS},
  {"gatherDegenerated", K_CONNECTOR::gatherDegenerated, METH_VARARGS},
  {"setIBCTransfers", K_CONNECTOR::setIBCTransfers, METH_VARARGS},
  {"setIBCTransfersD", K_CONNECTOR::setIBCTransfersD, METH_VARARGS},
  {"_setIBCTransfers", K_CONNECTOR::_setIBCTransfers, METH_VARARGS},
  {"_setIBCTransfersD", K_CONNECTOR::_setIBCTransfersD, METH_VARARGS},
  {"modifyBorders", K_CONNECTOR::modifyBorders, METH_VARARGS},
  {"applyBCOverlapsNG", K_CONNECTOR::applyBCOverlapsNG, METH_VARARGS},
  {"applyBCOverlapStruct", K_CONNECTOR::applyBCOverlapStruct, METH_VARARGS},
  {"getExtrapAbsCoefs", K_CONNECTOR::getExtrapAbsCoefs, METH_VARARGS},
  {"_getEmptyBCInfoNGON", K_CONNECTOR::_getEmptyBCInfoNGON, METH_VARARGS},
  {"_updateNatureForIBM",K_CONNECTOR::_updateNatureForIBM, METH_VARARGS},//on a zone, in place
  {"indiceToCoord2",K_CONNECTOR::indiceToCoord2, METH_VARARGS},//on a zone, in place
  {"correctCoeffList",K_CONNECTOR::correctCoeffList, METH_VARARGS},//on a zone, in place
  {"_blankClosestTargetCells",K_CONNECTOR::_blankClosestTargetCells, METH_VARARGS},
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
        "connector",
        NULL,
        sizeof(struct module_state),
        Pyconnector,
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
  PyMODINIT_FUNC PyInit_connector();
  PyMODINIT_FUNC PyInit_connector()
#else
  PyMODINIT_FUNC initconnector();
  PyMODINIT_FUNC initconnector()
#endif
  {
    __activation__ = K_KCORE::activation("0");
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("connector", Pyconnector);
#endif
    import_array();
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

/*
    Copyright 2013-2024 Onera.

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
#include "xcore.h"
#include "SplitElement/splitter.h"
#include "test/xmpi_t1.hpp"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pyxcore [] =
{
    {"test_all", xcore::test_all, METH_VARARGS}, // all xmpi tests
    {"splitElements", splitElements, METH_VARARGS},
    {"chunk2partNGon", K_XCORE::chunk2partNGon, METH_VARARGS},
    {"chunk2partElt", K_XCORE::chunk2partElt, METH_VARARGS},
    {"exchangeFields", K_XCORE::exchangeFields, METH_VARARGS},

    {"AdaptMesh_Init", K_XCORE::AdaptMesh_Init, METH_VARARGS},
    {"AdaptMesh_AssignRefData", K_XCORE::AdaptMesh_AssignRefData, METH_VARARGS},
    {"AdaptMesh_LoadBalance", K_XCORE::AdaptMesh_LoadBalance, METH_VARARGS},
    {"AdaptMesh_Adapt", K_XCORE::AdaptMesh_Adapt, METH_VARARGS},
    {"AdaptMesh_ExtractMesh", K_XCORE::AdaptMesh_ExtractMesh, METH_VARARGS},
    {"AdaptMesh_Exit", K_XCORE::AdaptMesh_Exit, METH_VARARGS},
    {"AdaptMesh_ExtractOwners", K_XCORE::AdaptMesh_ExtractOwners, METH_VARARGS},
    {"AdaptMesh_ExtractNeighbours", K_XCORE::AdaptMesh_ExtractNeighbours, METH_VARARGS},
    {"AdaptMesh_ExtractCellLevels", K_XCORE::AdaptMesh_ExtractCellLevels, METH_VARARGS},
    {"AdaptMesh_ExtractCellRanges", K_XCORE::AdaptMesh_ExtractCellRanges, METH_VARARGS},
    {"AdaptMesh_ExtractHaloCellLevels", K_XCORE::AdaptMesh_ExtractHaloCellLevels, METH_VARARGS},
    {"AdaptMesh_TagFaces", K_XCORE::AdaptMesh_TagFaces, METH_VARARGS},
    {"AdaptMesh_TriangulateFaces", K_XCORE::AdaptMesh_TriangulateFaces, METH_VARARGS},
    {"AdaptMesh_GeneratePrisms", K_XCORE::AdaptMesh_GeneratePrisms, METH_VARARGS},
    {"AdaptMesh_AdaptGeom", K_XCORE::AdaptMesh_AdaptGeom, METH_VARARGS},
    {"AdaptMesh_ExtractTaggedFaces", K_XCORE::AdaptMesh_ExtractTaggedFaces, METH_VARARGS},

    {"removeIntersectingKPlanes", K_XCORE::removeIntersectingKPlanes, METH_VARARGS},

    {"IntersectMesh_Init", K_XCORE::IntersectMesh_Init, METH_VARARGS},
    {"IntersectMesh_TriangulateFaceSet", K_XCORE::IntersectMesh_TriangulateFaceSet, METH_VARARGS},
    {"IntersectMesh_ExtractMesh", K_XCORE::IntersectMesh_ExtractMesh, METH_VARARGS},
    {"IntersectMesh_Exit", K_XCORE::IntersectMesh_Exit, METH_VARARGS},
    {"IntersectMesh_ExtractFaceSet", K_XCORE::IntersectMesh_ExtractFaceSet, METH_VARARGS},
    
    {"icapsule_init", K_XCORE::icapsule_init, METH_VARARGS},
    {"icapsule_adapt", K_XCORE::icapsule_adapt, METH_VARARGS},
    {"icapsule_intersect", K_XCORE::icapsule_intersect, METH_VARARGS},
    
    {"icapsule_extract_master", K_XCORE::icapsule_extract_master, METH_VARARGS},
    {"icapsule_extract_slave", K_XCORE::icapsule_extract_slave, METH_VARARGS},
    {"icapsule_extract_slaves", K_XCORE::icapsule_extract_slaves, METH_VARARGS},

    {"write_im", K_XCORE::write_im, METH_VARARGS},
    {"write_bim", K_XCORE::write_bim, METH_VARARGS},
    {"write_bim_s", K_XCORE::write_bim_s, METH_VARARGS},

    {"triangulate_skin", K_XCORE::triangulate_skin, METH_VARARGS},

    {"extractCell", K_XCORE::extractCell, METH_VARARGS},
    
    {"extractFacesFromPointTag", K_XCORE::extractFacesFromPointTag, METH_VARARGS},

    {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "xcore",
        NULL,
        -1,
        Pyxcore
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_xcore();
  PyMODINIT_FUNC PyInit_xcore()
#else
  PyMODINIT_FUNC initxcore();
  PyMODINIT_FUNC initxcore()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("xcore", Pyxcore);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

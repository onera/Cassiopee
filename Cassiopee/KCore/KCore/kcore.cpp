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
#include "kcore.h"

// ============================================================================
/* Dictionnary of all functions of the python module */
// ============================================================================
static PyMethodDef Pykcore [] =
{
  {"isNamePresent", K_KCORE::isNamePresent, METH_VARARGS},
  {"isCoordinateXPresent", K_KCORE::isCoordinateXPresent, METH_VARARGS},
  {"isCoordinateYPresent", K_KCORE::isCoordinateYPresent, METH_VARARGS},
  {"isCoordinateZPresent", K_KCORE::isCoordinateZPresent, METH_VARARGS},
  {"indiceStruct2Unstr", K_KCORE::indiceStruct2Unstr, METH_VARARGS},
  {"indiceStruct2Unstr2", K_KCORE::indiceStruct2Unstr2, METH_VARARGS},
  {"indiceFace2Connect", K_KCORE::indiceFace2Connect, METH_VARARGS},
  {"setOmpMaxThreads", K_KCORE::setOmpMaxThreads, METH_VARARGS},
  {"getOmpMaxThreads", K_KCORE::getOmpMaxThreads, METH_VARARGS},
  {"empty", K_KCORE::empty, METH_VARARGS},
  {"tester", K_KCORE::tester, METH_VARARGS},
  {"testerAcc", K_KCORE::testerAcc, METH_VARARGS},
  {"copyto", K_KCORE::copyto, METH_VARARGS},
  {"copyfrom", K_KCORE::copyfrom, METH_VARARGS},
  {NULL, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "kcore",
        NULL,
        -1,
        Pykcore
};
#endif

// ============================================================================
/* Init of module */
// ============================================================================
extern "C"
{
#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_kcore();
  PyMODINIT_FUNC PyInit_kcore()
#else
  PyMODINIT_FUNC initkcore();
  PyMODINIT_FUNC initkcore()
#endif
  {
    import_array();
#if PY_MAJOR_VERSION >= 3
    PyObject* module = PyModule_Create(&moduledef);
#else
    Py_InitModule("kcore", Pykcore);
#endif
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
  }
}

//=============================================================================
extern "C"
{
  void k6conv2center1_(const E_Int& ni, const E_Int& nj, const E_Int& nk, 
                       const E_Int& nfld, E_Float* fieldnode, 
                       E_Float* fieldcenter);
  void k6compunstrmetric_(E_Int& npts, E_Int& nbcell, 
                          E_Int& nedges, E_Int& nnodes, 
                          E_Int* cn, E_Float* x, E_Float* y, E_Float* z, 
                          E_Float* xint, E_Float* yint, E_Float* zint,
                          E_Float* surfx, E_Float* surfy, E_Float* surfz,
                          E_Float* snorm, E_Float* vol);
  void k6unstructsurf_(E_Int& npts, E_Int& nelts, E_Int& nedges, 
                       E_Int& nnodes, E_Int* cn, 
                       E_Float* xt, E_Float* yt, E_Float* zt,
                       E_Float* snx, E_Float* sny, E_Float* snz,
                       E_Float* surface);
  void k6compstructcellcenter_(E_Int& im, E_Int& jm, E_Int& km, 
                               E_Int& nbNode, E_Int& nbcell, 
                               const E_Float* xt, const E_Float* yt, 
                               const E_Float* zt, E_Float* bary);
  void k6comptetracellcenter_(const E_Int& npts, const E_Int& nelts, 
                              const E_Int* cn, const E_Float* xt,
                              const E_Float* yt, const E_Float* zt,
                              E_Float* bary);
}
//=============================================================================
/* Fonctions fortran declarees dans KCore mais non appelees dans KCore   
   Used to force some functions to belong to kcore library  */
//=============================================================================
void K_KCORE::testFooKCore()
{
  E_Int i=0; E_Float f=0.;
  
  k6conv2center1_(i, i, i, i, NULL, NULL);
  k6unstructsurf_(i, i, i, i, NULL, NULL, NULL, NULL,
                  NULL, NULL, NULL, NULL); 
  k6compunstrmetric_(i, i, i, i, 
                     NULL, NULL, NULL, NULL, 
                     NULL, NULL, NULL,
                     NULL, NULL, NULL,
                     NULL, NULL);
  k6compstructcellcenter_(i, i, i, i, i, 
                          NULL, NULL, 
                          NULL, NULL);
  k6comptetracellcenter_(i, i, 
                         NULL, NULL,
                         NULL, NULL,
                         NULL);
}


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
  
  void k6compminlengthoftetracell_(E_Int& npts, E_Int& indA, E_Int& indB,
                                   E_Int& indC, E_Int& indD, 
                                   E_Float* xt, E_Float* yt, E_Float* zt, 
                                   E_Float& minl);

  void k6compminlengthofcell_(E_Int& ni, E_Int& nj, E_Int& nk, E_Int& indA, 
                              E_Float* xt, E_Float* yt, E_Float* zt, 
                              E_Float& minl);

  void k6compintsurfofcell_(E_Int& ind, E_Int& ni, E_Int& nj, E_Int& nk, 
                            E_Float* xt, E_Float* yt, E_Float* zt, 
                            E_Float* surf);
  void k6compintsurf_(E_Int& ni, E_Int& nj, E_Int& nk, E_Int& npts, 
                      E_Int& nintt, E_Int& ninti, E_Int& nintij, 
                      E_Float* xt, E_Float* yt, E_Float* zt, 
                      E_Float* sx, E_Float* sy, E_Float* sz, E_Float* surfno);
  void k6structsurft_(E_Int& ni, E_Int& nj, E_Int& nk,
                      E_Int& ncells, E_Float* xt, E_Float* yt, E_Float* zt, 
                      E_Float* length);

  void k6structsurf1dt_(E_Int& ni, E_Int& nj, E_Int& nk,
                        E_Float* xt, E_Float* yt, E_Float* zt, 
                        E_Float* length);

  void k6compsurfofstructcell_(
    E_Int& ni, E_Int& nj, E_Int& nk,
    E_Int& indcell, E_Float* xt, E_Float* yt, E_Float* zt, 
    E_Float& surface);

  void k6compstructmetric_(
    E_Int& ni, E_Int& nj, E_Int& nk, E_Int& nbcell, 
    E_Int& nbInt, E_Int& nbInti, E_Int& nbIntj, 
    E_Int& nbIntk, E_Float* x, E_Float* y, E_Float* z, 
    E_Float* vol, E_Float* surfx, E_Float* surfy, E_Float* surfz,
    E_Float* snorm, E_Float* cix, E_Float* ciy, E_Float* ciz);

  void k6compunstrmetric_(E_Int& npts, E_Int& nbcell, 
                          E_Int& nedges, E_Int& nnodes, 
                          E_Int* cn, E_Float* x, E_Float* y, E_Float* z, 
                          E_Float* xint, E_Float* yint, E_Float* zint,
                          E_Float* surfx, E_Float* surfy, E_Float* surfz,
                          E_Float* snorm, E_Float* vol);

  void k6normstructsurft_(
    const E_Int& ni, const E_Int& nj, const E_Int& npts, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* nxt, E_Float* nyt, E_Float* nzt);

  void k6unstructsurf_(E_Int& npts, E_Int& nelts, E_Int& nedges, 
                       E_Int& nnodes, E_Int* cn, 
                       E_Float* xt, E_Float* yt, E_Float* zt,
                       E_Float* snx, E_Float* sny, E_Float* snz,
                       E_Float* surface);
  
  void k6normunstructsurf_(E_Int& nt, E_Int& nv, E_Int* cn, 
                           E_Float* coordx, E_Float* coordy, E_Float* coordz,
                           E_Float* nsurf);

  void k6compstructcellcenter_(E_Int& im, E_Int& jm, E_Int& km, 
                               E_Int& nbNode, E_Int& nbcell, 
                               const E_Float* xt, const E_Float* yt, 
                               const E_Float* zt, E_Float* bary);

  void k6comptetracellcenter_(const E_Int& npts, const E_Int& nelts, 
                              const E_Int* cn, const E_Float* xt,
                              const E_Float* yt, const E_Float* zt,
                              E_Float* bary);

  void k6compmindist_(const E_Int& ni1, const E_Int& nj1, 
                      const E_Float* x1, const E_Float* y1, 
                      const E_Float* z1, 
                      const E_Int& ni2, const E_Int& nj2, 
                      const E_Float* x2, const E_Float* y2, 
                      const E_Float* z2,
                      E_Int& ind1, E_Int& ind2, E_Float& dmin);
  
  void k6rectifynormals_(const E_Int& ni1, const E_Int& nj1, 
                         const E_Int& ind1,
                         const E_Float* x1, const E_Float* y1,
                         const E_Float* z1,
                         const E_Int& ni2, const E_Int& nj2,
                         const E_Int& ind2,
                         const E_Float* x2, const E_Float* y2, 
                         const E_Float* z2, const E_Float& distmin,
                         E_Int& notvalid, E_Int& isopp);

}
//=============================================================================
/* Fonctions fortran declarees dans KCore mais non appelees dans KCore   
   Used to force some functions to belong to kcore library  */
//=============================================================================
void K_KCORE::testFooKCore()
{
  E_Int i=0; E_Float f=0.;
  
  k6conv2center1_(i, i, i, i, NULL, NULL);
  k6compminlengthofcell_(i, i, i, i, 
                         NULL, NULL, NULL, f);
  k6compminlengthoftetracell_(i, i, i, i, i, 
                              NULL, NULL, NULL, f);
  k6compintsurfofcell_(i, i, i, i, 
                       NULL, NULL, NULL, NULL);
  k6compintsurf_(i, i, i, i, i, i, i, NULL, NULL, NULL,
                 NULL, NULL, NULL, NULL);
  k6structsurft_(i, i, i, i, NULL, NULL, NULL, NULL);
  k6structsurf1dt_(i, i, i, NULL, NULL, NULL, NULL);
  k6compsurfofstructcell_(i, i, i, i, NULL, NULL, NULL, f);

  k6unstructsurf_(i, i, i, i, NULL, NULL, NULL, NULL,
                  NULL, NULL, NULL, NULL); 

  k6normstructsurft_(i, i, i,  NULL, NULL, NULL,
                     NULL, NULL, NULL);
  k6normunstructsurf_(i, i, NULL, NULL, NULL, NULL, NULL);
  k6compstructmetric_(i, i, i, i, i, i, i, 
                      i, NULL, NULL, NULL, 
                      NULL, NULL, NULL, NULL,
                      NULL, NULL, NULL, NULL);
  k6compunstrmetric_(i, i, i, i, 
                     NULL, NULL, NULL, NULL, 
                     NULL, NULL, NULL,
                     NULL, NULL, NULL,
                     NULL, NULL);
  k6compmindist_(i, i, NULL, NULL, NULL, i, i, 
                 NULL, NULL, NULL, i, i, f);
  k6rectifynormals_(i, i, i, NULL, NULL,
                    NULL, i, i, i,
                    NULL, NULL, NULL, f,
                    i, i);
  k6compstructcellcenter_(i, i, i, i, i, 
                          NULL, NULL, 
                          NULL, NULL);

  k6comptetracellcenter_(i, i, 
                         NULL, NULL,
                         NULL, NULL,
                         NULL);
}


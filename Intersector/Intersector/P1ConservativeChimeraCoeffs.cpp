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

# include "intersector.h"
#include "Nuga/include/conservative_chimera.h"
#include <memory>
//# include "Numpy/Numpy.h"
//#include <sstream>
//#include <fstream>
//#include <iostream>
//#include "chrono.h"

using namespace std;
using namespace K_FLD;


//============================================================================
/* Blank cells defined in arrays by a Tetra mesh mask (as an input hook) */
//============================================================================
PyObject* K_INTERSECTOR::P1ConservativeChimeraCoeffs(PyObject* self, PyObject* args)
{
  typedef K_FLD::FldArrayF crd_t;
  typedef K_FLD::FldArrayI cnt_t;
  //typedef NUGA::NGON_BooleanOperator<crd_t, cnt_t> boolean_t;
  //typedef NUGA::P1_Conservative_Chimera<crd_t, cnt_t> chimera_t;
   
  PyObject *meshD, *meshR, *cellNR;
  
  char *varString1, *varString2, *varString3, *eltType1, *eltType2, *eltType3;
  E_Int ni, nj, nk;
  crd_t *fldR(0), *fldD(0), *fCelln(0);
  cnt_t *cnR(0), *cnD(0), *cCelln(0);
  
  if (!PyArg_ParseTuple(args, "OOO", &meshR, &cellNR, &meshD)) return NULL;
    
  /////////////////////////////////////////////////////////////////////////
  // Extraction des donnees
  
  E_Int res = K_ARRAY::getFromArray(meshD, varString1, fldD, ni, nj, nk, cnD, eltType1);

  if (res != 2 || strcmp(eltType1, "NGON") != 0 )
  {
    PyErr_SetString(PyExc_ValueError,
		   "P1ConservativeChimeraCoeffs: the donnor zone must be NGON.");
      return NULL;
  }
  
  std::unique_ptr<crd_t> afD(fldD); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<cnt_t> acD(cnD); // to avoid to call explicit delete at several places in the code.
  
  //std::cout << "bgm : " << crd->cols() << "/" << cnt->cols() << std::endl;
  //std::cout << "res : " << res << std::endl;
  
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                      "P1ConservativeChimeraCoeffs: donnor zone must contain coordinates.");
    return NULL;
  }

  res = K_ARRAY::getFromArray(meshR, varString2, fldR, ni, nj, nk, cnR, eltType2);

  std::unique_ptr<crd_t> afR(fldR); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<cnt_t> acR(cnR); // to avoid to call explicit delete at several places in the code.

  if (res != 2 || strcmp(eltType2, "NGON") != 0 )
  {
  	std::cout << "res : " << res << std::endl;
  	std::cout << "eltType2 : " << eltType2 << std::endl;

    PyErr_SetString(PyExc_ValueError,
		   "P1ConservativeChimeraCoeffs: the receiver zone must be NGON.");
      return NULL;
  }

  res = K_ARRAY::getFromArray(cellNR, varString3, fCelln, ni, nj, nk, cCelln, eltType3);
  
  std::unique_ptr<crd_t> afC(fCelln); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<cnt_t> acC(cCelln); // to avoid to call explicit delete at several places in the code.
   
  if (res == -1 || strcmp(eltType3, "NGON*") != 0 )
  {
    PyErr_SetString(PyExc_TypeError,
                      "P1ConservativeChimeraCoeffs: cellN must be specified at centers.");
    return NULL;
  }

  E_Int posc = K_ARRAY::isCellNatureField2Present(varString3);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                      "P1ConservativeChimeraCoeffs: celln variable not found for one structured array.");
    return NULL;
  }

  // Coordinates positions
  E_Int posDx = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posDy = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posDz = K_ARRAY::isCoordinateZPresent(varString1);
  E_Int posRx = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posRy = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posRz = K_ARRAY::isCoordinateZPresent(varString2);
  

  //  
  size_t sz = fCelln->getSize();
  std::vector<E_Int> cN(sz);
  for (size_t i = 0; i < sz; ++i) cN[i]=(*fCelln)[i];

  Vector_t<E_Int> dindices, xr, roids;
  Vector_t<E_Float> dcoeffs;

  /*std::cout << "nb points donnor : " << fldD->getSize() << std::endl;
  std::cout << "nb points receiver : " << fldR->getSize() << std::endl;
  std::cout << "nb cells donnor : " << (*cnD)[2 + (*cnD)[1]] << std::endl;
  std::cout << "nb cells receiver : " << (*cnR)[2 + (*cnR)[1]] << std::endl;
  std::cout << "pos R : " << posRx << " " << posRy << " " << posRz << std::endl;
  std::cout << "pos D : " << posDx << " " << posDy << " " << posDz << std::endl;
  std::cout << "cellN size : " << cN.size() << std::endl;*/
  
  E_Int err = NUGA::P1_CONSERVATIVE::compute_chimera_coeffs(*fldR, posRx, posRy, posRz, *cnR, 
  	                                                        *fldD, posDx, posDy, posDz, *cnD, 
  	                                                        cN,
  	                                                        dindices, dcoeffs, xr, roids);

  //std::cout << "return val : " << err << std::endl;
  
  if (err)
  {
    PyErr_SetString(PyExc_TypeError, "P1ConservativeChimeraCoeffs: failed.");
    return NULL;
  }

  PyObject *l(PyList_New(0)), *tpl;

  tpl = K_NUMPY::buildNumpyArray(&dindices[0], dindices.size(), 1, 0);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  tpl = K_NUMPY::buildNumpyArray(&dcoeffs[0], dcoeffs.size(), 1, 0);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  tpl = K_NUMPY::buildNumpyArray(&xr[0], xr.size(), 1, 0);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);

  tpl = K_NUMPY::buildNumpyArray(&roids[0], roids.size(), 1, 0);
  PyList_Append(l, tpl);
  Py_DECREF(tpl);
 
  return l;
}

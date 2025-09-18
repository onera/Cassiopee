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

# include "intersector.h"
#include "Nuga/include/conservative_chimera.h"
#include <memory>
#include "Nuga/include/mesh_t.hxx"


using namespace std;
using namespace K_FLD;
using namespace NUGA;

//============================================================================
/* Conservative interpolator*/
//============================================================================
PyObject* K_INTERSECTOR::P1ConservativeInterpolation(PyObject* self, PyObject* args)
{
  PyObject *meshD, *meshR, *fldD;
  
  if (!PYPARSETUPLE_(args, OOO_, &meshR, &meshD, &fldD)) return NULL;
  
  using crd_t = K_FLD::FloatArray;
  using cnt_t = K_FLD::IntArray;
  using zmesh_t = ph_mesh_t; 
  
  crd_t crdR;
  cnt_t cntR;
  char* varStringR, *eltTypeR;
  // Check array # 1
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray(meshR, varStringR, crdR, ni, nj, nk, cntR, eltTypeR);

  if (res != 2 || strcmp(eltTypeR, "NGON") != 0 )
  {
    PyErr_SetString(PyExc_ValueError,
       "P1ConservativeChimeraCoeffs: the donnor zone must be NGON.");
      return NULL;
  }

  //std::cout << "nb of fields : " << fR->getNfld() << std::endl;
  //std::cout << "vars : " << varStringR << std::endl;

  crd_t crdD;
  cnt_t cntD;
  char* varStringD, *eltTypeD;
  res = K_ARRAY::getFromArray(meshD, varStringD, crdD, ni, nj, nk, cntD, eltTypeD);

  if (res != 2 || strcmp(eltTypeR, "NGON") != 0 )
  {
    PyErr_SetString(PyExc_ValueError,
       "P1ConservativeChimeraCoeffs: the donnor zone must be NGON.");
      return NULL;
  }

  K_FLD::FldArrayF* fldsC;
  K_FLD::FldArrayI* cn;
  char* fvarStringsC, *feltType;
  res = K_ARRAY::getFromArray(fldD, fvarStringsC, fldsC, ni, nj, nk, cn, feltType);

  std::unique_ptr<K_FLD::FldArrayF> afldC(fldsC); // to avoid to call explicit delete at several places in the code.
  std::unique_ptr<K_FLD::FldArrayI> acn(cn); // to avoid to call explicit delete at several places in the code.

  E_Int nfields = fldsC->getNfld();
  
  std::vector<field> don_fields(nfields);

  for (E_Int j = 0; j < nfields; ++j)
    don_fields[j].f = fldsC->begin(j+1);
  //fixme : get gradients 

  /*std::cout << "res : " << res << std::endl;
  std::cout << "var : " << fvarStringsC << std::endl;
  std::cout << "field C : " << fldsC.rows() << "/" << fldsC.cols() << std::endl;
  std::cout << "cn : " << cn.rows() << "/" << cn.cols() << std::endl;*/


  zmesh_t mR(crdR, cntR), mD(crdD, cntD);
  E_Float RTOL{1.e-12};
  
  std::vector<std::vector<E_Float>> rec_fields;
  bool do_omp = true;

  NUGA::interpolate(mR, mD, RTOL, don_fields, rec_fields, do_omp);

  //  pushing out the received fields
  K_FLD::FloatArray farr(nfields, rec_fields[0].size());
  for (size_t i=0; i < rec_fields.size(); ++i)
  {
    std::vector<E_Float>& fld = rec_fields[i];
    for (size_t j = 0; j < fld.size(); ++j)farr(i, j) = fld[j];
  }

  PyObject* tpl = K_ARRAY::buildArray(farr, fvarStringsC, cntR, -1, feltType, false);
  //PyList_Append(l, tpl);
  //Py_DECREF(tpl);

  
  return tpl;

}

//============================================================================
/* Computes chimera  coeefs in a conservative manner*/
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
  
  if (!PYPARSETUPLE_(args, OOO_, &meshR, &cellNR, &meshD)) return NULL;
    
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

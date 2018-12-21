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
#include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert unstructured array to a hexaedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertUnstruct2Hexa(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, 
                                    f, nil, njl, nkl, cnl, eltType, true);
  if (res != 2) 
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2Hexa: array must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType, "QUAD") == 0 || strcmp(eltType, "HEXA") == 0 ||
      strcmp(eltType, "QUAD*") == 0 || strcmp(eltType, "HEXA*") == 0 ||
      strcmp(eltType, "BAR") == 0 || strcmp(eltType, "BAR*") == 0)
  {
    RELEASESHAREDU(array, f, cnl);
    return array;
  }
  E_Int nelts = cnl->getSize();
  char etString[64];

  if (strcmp(eltType, "TRI") == 0 || strcmp(eltType, "TRI*") == 0) 
  {
    E_Boolean loc = false; strcpy(etString, "QUAD");
    if (strcmp(eltType, "TRI*") == 0 ) {strcpy(etString, "QUAD*"); loc = true;}
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);

    E_Int eltt = 3; //QUAD
    E_Int nvert = 4;
    PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, f->getSize(), 
                                        nelts, eltt, etString, loc);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fn(f->getSize(), f->getNfld(), fnp, true); fn = *f;
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI ct(nelts, nvert, cnnp, true);
    E_Int* ct1 = ct.begin(1);
    E_Int* ct2 = ct.begin(2);
    E_Int* ct3 = ct.begin(3);
    E_Int* ct4 = ct.begin(4);

#pragma omp parallel for default(shared) if (nelts > 50)
    for (E_Int et = 0; et < nelts; et++)
    {
      ct1[et] = cn1[et];
      ct2[et] = cn2[et];
      ct3[et] = cn3[et];
      ct4[et] = cn3[et];
    }
    RELEASESHAREDU(array, f, cnl);
    return tpl;
  }
  else if (strcmp(eltType, "TETRA") == 0 || strcmp(eltType, "TETRA*") == 0) 
  {
    E_Boolean loc = false; strcpy(etString, "HEXA");
    if (strcmp(eltType, "TETRA*") == 0) {strcpy(etString, "HEXA*"); loc = true;}
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    E_Int* cn4 = cnl->begin(4);

    E_Int eltt = 7; //HEXA
    E_Int nvert = 8;
    PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, f->getSize(), nelts, eltt, etString, loc);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fn(f->getSize(), f->getNfld(), fnp, true); fn = *f;
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI ct(nelts, nvert, cnnp, true);
    E_Int* ct1 = ct.begin(1);
    E_Int* ct2 = ct.begin(2);
    E_Int* ct3 = ct.begin(3);
    E_Int* ct4 = ct.begin(4);
    E_Int* ct5 = ct.begin(5);
    E_Int* ct6 = ct.begin(6);
    E_Int* ct7 = ct.begin(7);
    E_Int* ct8 = ct.begin(8);

#pragma omp parallel for default(shared)
    for (E_Int et = 0; et < nelts; et++)
    {
      ct1[et] = cn1[et];
      ct2[et] = cn2[et];
      ct3[et] = cn3[et];
      ct4[et] = cn3[et];
      ct5[et] = cn4[et];
      ct6[et] = cn4[et];
      ct7[et] = cn4[et];
      ct8[et] = cn4[et];
    }
    RELEASESHAREDU(array, f, cnl);
    return tpl;
  }
  else if (strcmp(eltType, "PENTA") == 0 || strcmp(eltType, "PENTA*") == 0)
  {
    E_Boolean loc = false; strcpy(etString, "HEXA");
    if (strcmp(eltType, "PENTA*") == 0) {strcpy(etString, "HEXA*"); loc = true;}
    E_Int* cn1 = cnl->begin(1);
    E_Int* cn2 = cnl->begin(2);
    E_Int* cn3 = cnl->begin(3);
    E_Int* cn4 = cnl->begin(4);
    E_Int* cn5 = cnl->begin(5);
    E_Int* cn6 = cnl->begin(6);
    E_Int eltt = 7;//HEXA
    E_Int nvert = 8;
    PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, f->getSize(), nelts, eltt, etString, loc);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fn(f->getSize(), f->getNfld(), fnp, true); fn = *f;
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI ct(nelts, nvert, cnnp, true);
    E_Int* ct1 = ct.begin(1);
    E_Int* ct2 = ct.begin(2);
    E_Int* ct3 = ct.begin(3);
    E_Int* ct4 = ct.begin(4);
    E_Int* ct5 = ct.begin(5);
    E_Int* ct6 = ct.begin(6);
    E_Int* ct7 = ct.begin(7);
    E_Int* ct8 = ct.begin(8);

#pragma omp parallel for default(shared)
    for (E_Int et = 0; et < nelts; et++)
    {
      ct1[et] = cn1[et];
      ct2[et] = cn2[et];
      ct3[et] = cn3[et];
      ct4[et] = cn3[et];
      ct5[et] = cn4[et];
      ct6[et] = cn5[et];
      ct7[et] = cn6[et];
      ct8[et] = cn6[et];
    }
    RELEASESHAREDU(array, f, cnl);
    return tpl;
  }
  else 
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2Hexa: invalid element type.");
    return NULL;
  }
}

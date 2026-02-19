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

// fonctions optimisee de modification du cellN

#include "connector.h"

// ============================================================================
// modCellN1
// 0 -> -1; 2 -> 1
//=============================================================================
PyObject* K_CONNECTOR::_modCellN1(PyObject* self, PyObject* args)
{
  PyObject* array; char* cellNName;
  if (!PYPARSETUPLE_(args, O_ S_, &array, &cellNName))
      return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "modCellN1: invalid array.");
    return NULL;
  }

  E_Int posCellN = K_ARRAY::isNamePresent(cellNName, varString);
   
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "modCellN1: can't find cellN in array.");
    return NULL;
  }
  posCellN++;
    
  E_Int npts = f->getSize();
  E_Float* ft = f->begin(posCellN);

 #pragma omp parallel default(shared)
  {
    E_Float val;
#pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      val = ft[ind];
      if (K_FUNC::fEqualZero(val)) ft[ind] = -1.;
      else if (K_FUNC::fEqualZero(val - 2.)) ft[ind] = 1;  
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}

// ============================================================================
// modCellN2
// -1 -> 0
//=============================================================================
PyObject* K_CONNECTOR::_modCellN2(PyObject* self, PyObject* args)
{
  PyObject* array; char* cellNName;
  if (!PYPARSETUPLE_(args, O_ S_,
                    &array, &cellNName))
      return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "modCellN2: invalid array.");
    return NULL;
  }

  E_Int posCellN = K_ARRAY::isNamePresent(cellNName, varString);
   
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "modCellN2: can't find cellN in array.");
    return NULL;
  }
  posCellN++;
    
  E_Int npts = f->getSize();
  E_Float* ft = f->begin(posCellN);

 #pragma omp parallel default(shared)
  {
    E_Float val;
#pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      val = ft[ind];
      if (val == -1.) ft[ind] = 0;
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}

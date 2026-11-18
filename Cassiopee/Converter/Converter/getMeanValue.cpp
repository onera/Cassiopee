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
#include "converter.h"

//=============================================================================
/* Get mean value of an array */
//=============================================================================
PyObject* K_CONVERTER::getMeanValue(PyObject* self, PyObject* args)
{
  PyObject* array; char* varName;
  if (!PYPARSETUPLE_(args, O_ S_, &array, &varName)) return NULL;

  // Check array
  E_Int im, jm, km;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);

  E_Int posvar = -1;
  if (res == 1 || res == 2)
  {
    posvar = K_ARRAY::isNamePresent(varName, varString);
    if (posvar == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getMeanValue: variable %s is not present in array.", 
              varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }
    posvar++;

    E_Int fsize = f->getSize();
    if (fsize == 0)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getArgMean: array for variable %s is empty.", varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }

    E_Float* fp = f->begin(posvar);
    E_Float meanVal = 0.;

#pragma omp parallel for default(shared) reduction(+:meanVal) 
    for (E_Int i = 0; i < fsize; i++) meanVal += fp[i];

    meanVal = meanVal / E_Float(fsize);
    RELEASESHAREDB(res, array, f, cn);
    return Py_BuildValue(R_, meanVal);
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getMeanValue: argument is not an array.");
    return NULL;
  }
}

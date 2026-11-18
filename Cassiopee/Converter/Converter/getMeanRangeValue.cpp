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
/* Get mean value of a range in an array */
//=============================================================================
PyObject* K_CONVERTER::getMeanRangeValue(PyObject* self, PyObject* args)
{
  PyObject* array; char* varName;
  E_Float rmin, rmax;
  if (!PYPARSETUPLE_(args, O_ S_ RR_, &array, &varName, &rmin, &rmax)) return NULL;

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
      sprintf(error, "getMeanRangeValue: variable %s is not present in array.", 
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
      sprintf(error, "getMeanRangeValue: array for variable %s is empty.", varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }

    E_Float* fp = f->begin(posvar);
    E_Float meanVal = 0.;
    
    std::vector<E_Float> vals;
    vals.insert(vals.end(), fp, fp+fsize);
    
    std::sort(vals.begin(), vals.end());
    
    size_t sz1 = size_t(E_Float(fsize)*rmin);
    size_t sz2 = size_t(E_Float(fsize)*rmax);
    
    for (size_t i = sz1; i < sz2; ++i) meanVal += vals[i];

    meanVal = meanVal / E_Float(1+sz2-sz1);

    RELEASESHAREDB(res, array, f, cn);
    return Py_BuildValue(R_, meanVal);
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getMeanRangeValue: argument is not an array.");
    return NULL;
  }
}


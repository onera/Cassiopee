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

// set NAN value to a given value in numpy

# include "converter.h"

using namespace K_FLD;
using namespace std;

#include <cmath>

//==============================================================================
PyObject* K_CONVERTER::setNANValuesAt(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float val;
  if (!PYPARSETUPLE_(args, O_ R_, &array, &val)) return NULL;
  
  FldArrayF* f; 
  E_Int res = K_NUMPY::getFromNumpyArray(array, f);

  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setNANValuesAt: numpy is invalid.");
    return NULL;
  }
  
  E_Int npts = f->getSize();
  E_Float* fp = f->begin();

  #pragma omp parallel
  {
    #pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      if (std::isfinite(fp[i]) == false) fp[i] = val;
    }
  }

  RELEASESHAREDN(array, f);

  Py_INCREF(Py_None);
  return Py_None;
}
    

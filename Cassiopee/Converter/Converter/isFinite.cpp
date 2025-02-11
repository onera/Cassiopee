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

// test if NAN or infinite value in numpy

# include "converter.h"

using namespace K_FLD;
using namespace std;

#include <cmath>

//==============================================================================
PyObject* K_CONVERTER::isFinite(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;
  
  FldArrayF* f; 
  E_Int res = K_NUMPY::getFromNumpyArray(array, f, true);

  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "isFinite: numpy is invalid.");
    return NULL;
  }
  
  E_Int npts = f->getSize();
  E_Float* fp = f->begin();
  E_Int nthreads = __NUMTHREADS__;
  E_Int* ret = new E_Int [nthreads];
  for (E_Int i = 0; i < nthreads; i++) ret[i] = 0;

  #pragma omp parallel
  {
    E_Int ithread = __CURRENT_THREAD__;
    #pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      if (std::isfinite(fp[i]) == false) ret[ithread] = 1;
    }
  }
  E_Int found = 0;
  for (E_Int i = 0; i < nthreads; i++) found = max(found, ret[i]);

  RELEASESHAREDN(array, f);
  delete [] ret;

  return Py_BuildValue("l", found);
}
    

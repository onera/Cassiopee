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

// copy a numpy from device (openmp5)
// numpy must already exists in cpu memory
#include "kcore.h"

PyObject* K_KCORE::copyfrom(PyObject* self, PyObject* args)
{
  PyObject* numpyArray; 
  if (!PYPARSETUPLE_(args, O_ II_, &numpyArray))
  {
    return NULL;
  }

  FldArrayF* f;
  K_NUMPY::getFromNumpyArray(numpyArray, f); 
  //E_Float* ipttarget = f->begin();
  //E_Int sizetot = f->getSize();

#ifdef _OPENACC
//#pragma omp target update from (ipttarget[:sizetot])
#endif

  RELEASESHAREDN(numpyArray, f);

  Py_INCREF(Py_None);
  return Py_None;
}

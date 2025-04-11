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
#include <limits>
#include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert a shifted Parent Element numpy to a Parent Element numpy where the
   first element has an index of 1 */
//=============================================================================
PyObject* K_CONVERTER::adaptShiftedPE2PE(PyObject* self, PyObject* args)
{
  PyObject* arrayFE; 
  E_Int nelts;
  if (!PYPARSETUPLE_(args, O_ I_,
                     &arrayFE, &nelts)) return NULL;

  // Check numpy (FE)
  FldArrayI* cFE;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayFE, cFE, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptShiftedPE2PE: cFE numpy is invalid.");
    return NULL;
  }
  
  E_Int nfaces = cFE->getSize();
  E_Int* cFEg = cFE->begin(1);
  E_Int* cFEd = cFEg + nfaces;
  
  const E_Int numThreads = __NUMTHREADS__;
  // Use padding to avoid False Sharing
  const E_Int PAD = CACHELINE/sizeof(E_Int);  // Elements per cache line
  
  E_Int eIdxMin = std::numeric_limits<E_Int>::max();
  E_Int threadMin[numThreads*PAD];
  for (E_Int i = 0; i < numThreads; i++) threadMin[i*PAD] = eIdxMin;
  
  // Find the smallest element index in the FE connectivity, eIdxMin
  #pragma omp parallel num_threads(numThreads)
  {
    E_Int threadPos = __CURRENT_THREAD__;
    threadPos *= PAD;
    E_Int e1, e2;
    #pragma omp for
    for (E_Int i = 0; i < nfaces; i++)
    {
      e1 = cFEg[i]; e2 = cFEd[i];
      if (e1 < 1) e1 = eIdxMin;
      if (e2 < 1) e2 = eIdxMin;
      threadMin[threadPos] = K_FUNC::E_min(threadMin[threadPos], e1, e2);
    }
  }
  for (E_Int i = 0; i < numThreads; i++)
  {
      eIdxMin = K_FUNC::E_min(threadMin[i*PAD], eIdxMin);
  }
  
  if (eIdxMin > 1 and eIdxMin == nelts+1)
  {
    // Shift element indices such that they start at 1
    PyObject* tpl = K_NUMPY::buildNumpyArray(nfaces, 2, 1, 1);
    E_Int* pE = K_NUMPY::getNumpyPtrI(tpl);
    
    //#pragma omp parallel for
    for (E_Int i = 0; i < 2*nfaces; i++)
    {
      E_Int elt = cFEg[i];
      if (elt > 0) elt += 1-eIdxMin;
      pE[i] = elt;
    }
   
    RELEASESHAREDN(arrayFE, cFE);
    return tpl;
  }
  
  // The smallest element index is 1, do nothing
  RELEASESHAREDN(arrayFE, cFE);
  Py_INCREF(Py_None);
  return Py_None;
}

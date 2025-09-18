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
#include "converter.h"

// ============================================================================
/* Convert a strand array to a penta mesh */
// ============================================================================
PyObject* K_CONVERTER::convertStrand2Penta(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // input strand array is stored as a TRI array
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, f,
                               nil, njl, nkl, cnl, eltType);

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertStrand2Penta: array must be unstructured TRI (strand).");
    return NULL; 
  }

  if (res == 2)
  {
    if (strcmp(eltType, "TRI") != 0)
    {
      RELEASESHAREDU(array, f, cnl);
      PyErr_SetString(PyExc_TypeError, 
                      "convertStrand2Penta: array must be a single TRI connectivity (strand).");
      return NULL; 
    }
  }

  E_Int npts = f->getSize(); // np * nk
  E_Int nfld = f->getNfld(); // nbre de champs

  E_Int api = f->getApi();
  
  // Get TRI connectivity
  FldArrayI& cm = *(cnl->getConnect(0));
  E_Int ne = cm.getSize();
  E_Int nvpe = cm.getNfld(); // must be 3

  if (nvpe != 3)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "convertStrand2Penta: element type on the skin must be TRI.");
    return NULL; 
  }

  // Guess nk - connectivity reference the first vertices
  const E_Int numThreads = __NUMTHREADS__;
  E_Int* threadVmax = new E_Int [numThreads];
  for (E_Int i = 0; i < numThreads; i++){ threadVmax[i]=0; }

  #pragma omp parallel num_threads(numThreads)
  {
    E_Int threadId = __CURRENT_THREAD__;

    #pragma omp for
    for (E_Int i = 0; i < ne; i++)
    {
      threadVmax[threadId] = std::max(threadVmax[threadId], cm(i,1));
      threadVmax[threadId] = std::max(threadVmax[threadId], cm(i,2));
      threadVmax[threadId] = std::max(threadVmax[threadId], cm(i,3));
    }
  }

  E_Int vmax = 0;
  for (E_Int i = 0; i < numThreads; i++) vmax = std::max(vmax, threadVmax[i]);
  delete [] threadVmax;

  E_Int np = vmax;
  E_Int nk = npts / np;
  if (np*nk != npts)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "convertStrand2Penta: array is not a strand grid.");
    return NULL; 
  }
  printf("detected np=" SF_D_ " nk=" SF_D_ "\n", np, nk);

  // maillage penta : 
  // npts = np*nk
  // nelts = ne*(nk-1)
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, np*nk, ne*(nk-1),
                                       "PENTA", false, api);
  FldArrayF* f2; FldArrayI* cnl2;
  K_ARRAY::getFromArray3(tpl, f2, cnl2);
  FldArrayI& cm2 = *(cnl2->getConnect(0));

  #pragma omp parallel
  {
    E_Int ind1, ind2, ind3;
    
    // copy fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* fp2 = f2->begin(n);

      #pragma omp for
      for (E_Int i = 0; i < npts; i++) fp2[i] = fp[i];
    }

    // copy connect
    for (E_Int k = 0; k < nk-1; k++)
    {
      #pragma omp for
      for (E_Int i = 0; i < ne; i++)
      {
        ind1 = i+k*ne; ind2 = k*np; ind3 = (k+1)*np;
        cm2(ind1,1) = cm(i,1)+ind2;
        cm2(ind1,2) = cm(i,2)+ind2;
        cm2(ind1,3) = cm(i,3)+ind2;
        cm2(ind1,4) = cm(i,1)+ind3;
        cm2(ind1,5) = cm(i,2)+ind3;
        cm2(ind1,6) = cm(i,3)+ind3;
      }
    }
  }

  RELEASESHAREDU(array, f, cnl);
  RELEASESHAREDU(tpl, f2, cnl2);
  return tpl;
}

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
#include "kcore.h"
#include "converter.h"
using namespace K_FLD;

//=============================================================================
/* Return a subset (slice) of an unstructured NGon array. The first numpy array
   contains the vertex indices for all input faces and the second numpy array
   contains the face offsets (shape nfaces+1) */
//=============================================================================
PyObject* K_CONVERTER::sliceNGonFaces(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* arrayFPL;
  if (!PYPARSETUPLE_(args, OO_, &array, &arrayFPL)) return NULL;

  // Check numpy
  FldArrayI* facePL;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayFPL, facePL);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "sliceNGonFaces: facePL numpy is invalid.");
    return NULL;
  }
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "sliceNGonFaces: for unstructured NGon arrays only.");
    RELEASESHAREDN(arrayFPL, facePL);
    RELEASESHAREDS(array, f);
    return NULL;
  }
  else if (res == 2)
  { 
    PyObject* tpl = NULL;
    if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0)
    {
      E_Int* ngon = cn->getNGon();
      E_Int* indPG = cn->getIndPG();
      E_Int nfaces = facePL->getSize();
      E_Int* fplp = facePL->begin();
      
      E_Int nthreads = __NUMTHREADS__;
      std::vector<E_Int> tnpts(nthreads+1, 0);
      
      // Compute the total number of (non-unique) vertices for all input faces
      #pragma omp parallel num_threads(nthreads)
      {
        E_Int nv, fidx;
        E_Int threadId = __CURRENT_THREAD__;
        
        // Manual partitioning of the for loop
        E_Int chunk = nfaces/nthreads;
        E_Int start = threadId*chunk;
        E_Int end = (threadId == nthreads - 1) ? nfaces : (threadId + 1)*chunk;
        
        for (E_Int i = start; i < end; i++)
        {
          fidx = fplp[i]-1;
          //E_Int* face = 
          cn->getFace(fidx, nv, ngon, indPG);
          tnpts[threadId+1] += nv;
        }
      }
      
      for (E_Int i = 1; i < nthreads+1; i++) tnpts[i] += tnpts[i-1];
      E_Int npts = tnpts[nthreads];
      
      // Build the output numpy arrays
      PyObject* tplv = K_NUMPY::buildNumpyArray(npts, 1, 1);
      E_Int* faceVertices = K_NUMPY::getNumpyPtrI(tplv);
      PyObject* tplf = K_NUMPY::buildNumpyArray(nfaces+1, 1, 1);
      E_Int* faceOffset = K_NUMPY::getNumpyPtrI(tplf);
      
      // Fill the output numpy arrays
      #pragma omp parallel num_threads(nthreads)
      {
        E_Int nv, fidx;
        E_Int threadId = __CURRENT_THREAD__;
        E_Int tOffset = tnpts[threadId];
        
        // Manual partitioning of the for loop
        E_Int chunk = nfaces/nthreads;
        E_Int start = threadId*chunk;
        E_Int end = (threadId == nthreads - 1) ? nfaces : (threadId + 1)*chunk;
        
        for (E_Int i = start; i < end; i++)
        {
          fidx = fplp[i]-1;
          E_Int* face = cn->getFace(fidx, nv, ngon, indPG);
          faceOffset[i] = tOffset;
          for (E_Int p = 0; p < nv; p++) faceVertices[tOffset + p] = face[p];
          tOffset += nv;
        }
      }
      faceOffset[nfaces] = npts;
      
      RELEASESHAREDN(arrayFPL, facePL);
      RELEASESHAREDU(array, f, cn);
      return Py_BuildValue("(OO)", tplv, tplf);
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "sliceNGonFaces: for unstructured NGon arrays only.");
    }

    RELEASESHAREDN(arrayFPL, facePL);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "sliceNGonFaces: unrecognised type of array.");
    return NULL;
  }
}

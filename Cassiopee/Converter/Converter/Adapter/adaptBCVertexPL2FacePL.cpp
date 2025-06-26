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
/* Adapt a vertex point list numpy to a face point list numpy for unstructured 
   arrays. All face points must be tagged for a face to be selected. */
//=============================================================================
PyObject* K_CONVERTER::adaptBCVertexPL2FacePL(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* arrayVPL;
  if (!PYPARSETUPLE_(args, OO_, &array, &arrayVPL)) return NULL;

  // Check numpy (vertex point list)
  FldArrayI* vertexPL;
  E_Int res = K_NUMPY::getFromPointList(arrayVPL, vertexPL, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptBCVertexPL2FacePL: vertexPL numpy is invalid.");
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
                    "adaptBCVertexPL2FacePL: for unstructured arrays only.");
    RELEASESHAREDN(arrayVPL, vertexPL);
    RELEASESHAREDS(array, f);
    return NULL;
  }
  else if (res == 2)
  { 
    PyObject* tpl = NULL;
    E_Int npts = f->getSize();
    if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0)
    {
      tpl = adaptBCVertexPL2FacePL_NGON(cn, vertexPL, npts);
    }
    else
    {
      //tpl = adaptBCVertexPL2FacePL_ME(cn, vertexPL, npts);
      PyErr_SetString(PyExc_TypeError,
                      "adaptBCVertexPL2FacePL: not implemented yet for ME arrays.");
      return NULL;
    }

    RELEASESHAREDN(arrayVPL, vertexPL);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptBCVertexPL2FacePL: unrecognised type of array.");
    return NULL;
  }
}

PyObject* K_CONVERTER::adaptBCVertexPL2FacePL_NGON(FldArrayI* cn, FldArrayI* vpl, E_Int npts)
{
  E_Int* ngon = cn->getNGon();
  E_Int* indPG = cn->getIndPG();
  E_Int nfaces = cn->getNFaces();
  E_Int nptsVPL = vpl->getSize();
  E_Int* vplp = vpl->begin();
  
  // Tag vertices that are in the vertex point list
  FldArrayI tmap(npts); tmap.setAllValuesAt(-1); E_Int* tmapP = tmap.begin();
  #pragma omp parallel for
  for (E_Int i = 0; i < nptsVPL; i++) tmapP[vplp[i]-1] = 1;
  
  E_Int nthreads = __NUMTHREADS__;
  std::vector<std::vector<E_Int> > faceList(nthreads);

  // Store faces that are fully covered by the selected nodes
  #pragma omp parallel num_threads(nthreads)
  {
    E_Int threadId = __CURRENT_THREAD__;
    // Ensure found face indices are sorted by manually partitioning the for loop
    // (needed for reproducibility)
    E_Int chunk = nfaces/nthreads;
    E_Int start = threadId*chunk;
    E_Int end = (threadId == nthreads - 1) ? nfaces : (threadId + 1)*chunk;
    
    faceList[threadId].reserve(nptsVPL/(4*nthreads));  // initial guess
    
    E_Int nv, indv, compt;
    for (E_Int i = start; i < end; i++)
    {
      E_Int* face = cn->getFace(i, nv, ngon, indPG);
      compt = 0;
      for (E_Int p = 0; p < nv; p++)
      {
        indv = face[p]-1;
        if (tmapP[indv] != -1) compt++;
      }
      if (nv == compt) faceList[threadId].push_back(i+1);
    }
  }

  // Compute total number of BC faces and offsets
  std::vector<E_Int> ncumFacesPL(nthreads+1); ncumFacesPL[0] = 0;
  for (E_Int i = 0; i < nthreads; i++) ncumFacesPL[i+1] = ncumFacesPL[i] + faceList[i].size();

  // Fill numpy array
  PyObject* tpl = K_NUMPY::buildNumpyArray(ncumFacesPL[nthreads], 1, 1);
  E_Int* facePL = K_NUMPY::getNumpyPtrI(tpl);
  #pragma omp parallel num_threads(nthreads)
  {
    E_Int threadId = __CURRENT_THREAD__;
    E_Int offset = ncumFacesPL[threadId];
    E_Int size = ncumFacesPL[threadId+1] - offset;
    for (E_Int i = 0; i < size; i++) facePL[offset+i] = faceList[threadId][i];
  }

  return tpl;
}

PyObject* K_CONVERTER::adaptBCVertexPL2FacePL_ME(FldArrayI* cn, FldArrayI* vpl, E_Int npts)
{
  return NULL;
}

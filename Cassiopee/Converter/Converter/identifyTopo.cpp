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

// Topological face identification

# include "converter.h"
#include <unordered_map>

using namespace K_FLD;

// ============================================================================
// Identify face indices topologically between two arrays.
// Return indices into the set of faces that were searched.
// array1: reference array containing faces to be indexed. Must be NGON.
// array2: target array whose elements (BE/ME or NGON) are to be matched
//         against faces of array1.
// vertexIndir: renumbering of vertices of array2 into the indexing space of
//              array1 (0-based indexing). Used when both arrays share topology
//              but differ in vertex ordering, eg, after G.close.
// faceIndices: subset of faces from array1 to be considered.
// ============================================================================
PyObject* K_CONVERTER::identifyFacesTopo(PyObject* self, PyObject* args)
{
  PyObject *array1, *array2;
  PyObject* vertexIndir = NULL;
  PyObject* faceIndices = NULL;
  if (!PYPARSETUPLE_(args, OOOO_,
                     &array1, &array2, &vertexIndir, &faceIndices)) return NULL;

  // Check array1
  E_Int ni1, nj1, nk1, res1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  res1 = K_ARRAY::getFromArray3(array1, varString1, 
                                f1, ni1, nj1, nk1, cn1, eltType1);

  // Check array2
  E_Int ni2, nj2, nk2, res2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  res2 = K_ARRAY::getFromArray3(array2, varString2, 
                                f2, ni2, nj2, nk2, cn2, eltType2);

  if (res1 != 2 || res2 != 2)
  {
    if (res1 != 2)
      PyErr_SetString(PyExc_TypeError,
                      "identifyFacesTopo: array1 is invalid.");
    if (res2 != 2) 
      PyErr_SetString(PyExc_TypeError,
                      "identifyFacesTopo: array2 is invalid.");
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);
    return NULL;
  }

  // Get vertex indirection (0-based)
  FldArrayI* listVertexIndir;
  if (vertexIndir != NULL)
  {
    E_Int resi = K_NUMPY::getFromNumpyArray(vertexIndir, listVertexIndir);
    if (resi == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "identifyFacesTopo: vertex indirection must be a numpy of "
                      "integers.");
      RELEASESHAREDU(array1, f1, cn1);
      RELEASESHAREDU(array2, f2, cn2);
      return NULL;
    }
  }
  
  // Get list of face indices of array1
  FldArrayI* listFaceIndices;
  if (faceIndices != NULL)
  {
    E_Int resi = K_NUMPY::getFromNumpyArray(faceIndices, listFaceIndices);
    if (resi == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "identifyFacesTopo: subset of face indices must be a "
                      "numpy of integers.");
      if (vertexIndir != NULL) RELEASESHAREDN(vertexIndir, listVertexIndir);
      RELEASESHAREDU(array1, f1, cn1);
      RELEASESHAREDU(array2, f2, cn2);
      return NULL;
    }
  }

  // Hash faces of array1
  E_Int dim1, dim2;
  std::unordered_map<Topology, E_Int, BernsteinHash<Topology> > faceMap;

  if (K_STRING::cmp(eltType1, 4, "NGON") == 0)
  {
    Topology F;
    E_Int nv, nfaces1, fidx;
    dim1 = cn1->getDim();
    E_Int *ngon1 = cn1->getNGon(), *indPG1 = cn1->getIndPG();
    if (faceIndices == NULL)
    {
      // Hash all faces of array1
      nfaces1 = cn1->getNFaces();
      faceMap.reserve(nfaces1);
      for (E_Int i = 0; i < nfaces1; i++)
      {
        E_Int* face = cn1->getFace(i, nv, ngon1, indPG1);
        F.set(face, nv);
        faceMap.insert(std::make_pair(F, i));
      }
    }
    else
    {
      // Hash a subset of all faces of array1
      E_Int* listFaceIndicesp = listFaceIndices->begin();
      nfaces1 = listFaceIndices->getSize();
      faceMap.reserve(nfaces1);
      for (E_Int i = 0; i < nfaces1; i++)
      {
        fidx = listFaceIndicesp[i]-1;
        E_Int* face = cn1->getFace(fidx, nv, ngon1, indPG1);
        F.set(face, nv);
        faceMap.insert(std::make_pair(F, i));  // position inside the subset
      }
    }
  }
  else  // BE/ME
  {
    PyErr_SetString(PyExc_TypeError,
                    "identifyFacesTopo: BE/ME array1 not implemented yet.");
    RELEASESHAREDU(array1, f1, cn1);
    RELEASESHAREDU(array2, f2, cn2);
    if (vertexIndir != NULL) RELEASESHAREDN(vertexIndir, listVertexIndir);
    if (faceIndices != NULL) RELEASESHAREDN(faceIndices, listFaceIndices);
    return NULL;
  }

  PyObject* fmap = NULL;

  if (K_STRING::cmp(eltType2, 4, "NGON") == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "identifyFacesTopo: NGON array2 not implemented yet.");
    RELEASESHAREDU(array1, f1, cn1);
    RELEASESHAREDU(array2, f2, cn2);
    if (vertexIndir != NULL) RELEASESHAREDN(vertexIndir, listVertexIndir);
    if (faceIndices != NULL) RELEASESHAREDN(faceIndices, listFaceIndices);
    return NULL;
  }
  else  // BE/ME
  {
    dim2 = K_CONNECT::getDimME(eltType2);
    if (dim1 != dim2+1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "identifyFacesTopo: dimension of array2 must be one less "
                      "that of array1.");
      RELEASESHAREDU(array1, f1, cn1);
      RELEASESHAREDU(array2, f2, cn2);
      if (vertexIndir != NULL) RELEASESHAREDN(vertexIndir, listVertexIndir);
      if (faceIndices != NULL) RELEASESHAREDN(faceIndices, listFaceIndices);
      return NULL;
    }

    // Compute the cumulative number of elts per conn. of array2
    // Because dim(array2) = dim(array1) - 1, 'elements' of array2 are faces
    E_Int nc2 = cn2->getNConnect();
    std::vector<E_Int> cumnepc(nc2+1); cumnepc[0] = 0;
    for (E_Int ic = 0; ic < nc2; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn2->getConnect(ic));
      E_Int nelts = cm.getSize();
      cumnepc[ic+1] = cumnepc[ic] + nelts;
    }
    E_Int ntotElts = cumnepc[nc2];

    fmap = K_NUMPY::buildNumpyArray(ntotElts, 1, 1);
    E_Int* fmapp = K_NUMPY::getNumpyPtrI(fmap);

    // Loop over 'elements' (faces) of array2 and find indices of matching
    // faces in array1
    #pragma omp parallel
    {
      for (E_Int ic = 0; ic < nc2; ic++)
      {
        E_Int eidx, vidx;
        FldArrayI& cm = *(cn2->getConnect(ic));
        E_Int nelts = cm.getSize();
        E_Int nvpe = cm.getNfld();
        Topology F;
        std::vector<E_Int> eltloc(nvpe);

        if (vertexIndir == NULL)
        {
          #pragma omp for
          for (E_Int i = 0; i < nelts; i++)
          {
            eidx = cumnepc[ic] + i;
            for (E_Int j = 0; j < nvpe; j++) eltloc[j] = cm(i,j+1);
            F.set(eltloc, nvpe);
            auto resf = faceMap.find(F);
            if (resf != faceMap.end()) fmapp[eidx] = resf->second;  // 0-based
            else fmapp[eidx] = -1;
          }
        }
        else
        {
          E_Int* listVertexIndirp = listVertexIndir->begin();
          #pragma omp for
          for (E_Int i = 0; i < nelts; i++)
          {
            eidx = cumnepc[ic] + i;
            for (E_Int j = 0; j < nvpe; j++)
            {
              vidx = cm(i,j+1);
              eltloc[j] = listVertexIndirp[vidx-1] + 1;  // reindex array2 vertices
            }
            F.set(eltloc, nvpe);
            auto resf = faceMap.find(F);
            if (resf != faceMap.end())
              fmapp[eidx] = resf->second;  // index of matching face in array1
            else fmapp[eidx] = -1;  // no topological match found
          }
        }
      }
    }
  }

  RELEASESHAREDU(array1, f1, cn1);
  RELEASESHAREDU(array2, f2, cn2);
  if (vertexIndir != NULL) RELEASESHAREDN(vertexIndir, listVertexIndir);
  if (faceIndices != NULL) RELEASESHAREDN(faceIndices, listFaceIndices);
  return fmap;
}

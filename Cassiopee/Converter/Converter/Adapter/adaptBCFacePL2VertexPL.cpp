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
#include <unordered_set>
#include <algorithm> // for std::copy

#include "kcore.h"
#include "converter.h"
using namespace K_FLD;

//=============================================================================
/* Adapt a face point list numpy to a vertex point list numpy for unstructured 
   arrays */
//=============================================================================
PyObject* K_CONVERTER::adaptBCFacePL2VertexPL(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* arrayFPL;
  if (!PYPARSETUPLE_(args, OO_, &array, &arrayFPL)) return NULL;

  // Check numpy (face point list)
  FldArrayI* facePL;
  E_Int res = K_NUMPY::getFromPointList(arrayFPL, facePL);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptBCFacePL2VertexPL: facePL numpy is invalid.");
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
                    "adaptBCFacePL2VertexPL: for unstructured arrays only.");
    RELEASESHAREDN(arrayFPL, facePL);
    RELEASESHAREDS(array, f);
    return NULL;
  }
  else if (res == 2)
  { 
    PyObject* tpl = NULL;
    if (K_STRING::cmp(eltType, 4, "NGON") == 0)
    {
      tpl = adaptBCFacePL2VertexPL_NGON(cn, facePL);
    }
    else
    {
      E_Int npts = f->getSize();
      tpl = adaptBCFacePL2VertexPL_ME(cn, eltType, npts, facePL);
    }

    RELEASESHAREDN(arrayFPL, facePL);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptBCFacePL2VertexPL: unrecognised type of array.");
    return NULL;
  }
}

PyObject* K_CONVERTER::adaptBCFacePL2VertexPL_NGON(FldArrayI* cn, FldArrayI* fpl)
{
  E_Int* fplp = fpl->begin();
  E_Int size = fpl->getSize();
  E_Int nv, npts2, fidx;
  E_Int* ngon = cn->getNGon();
  E_Int* indPG = cn->getIndPG();
  
  std::unordered_set<E_Int> vertexSet;
  
  // Loop over all faces of the face point list and fill vertex set
  for (E_Int i = 0; i < size; i++)
  {
    fidx = fplp[i]-1;
    E_Int* face = cn->getFace(fidx, nv, ngon, indPG);
    for (E_Int j = 0; j < nv; j++) vertexSet.insert(face[j]);
  }
  npts2 = vertexSet.size();
  
  PyObject* tpl = K_NUMPY::buildNumpyArray(npts2, 1, 1);
  E_Int* vertexPL = K_NUMPY::getNumpyPtrI(tpl);
  std::copy(vertexSet.begin(), vertexSet.end(), vertexPL);
  return tpl;
}

PyObject* K_CONVERTER::adaptBCFacePL2VertexPL_ME(
  FldArrayI* cn, const char* eltType, const E_Int npts, FldArrayI* fpl
)
{
  E_Int* fplp = fpl->begin();
  E_Int size = fpl->getSize();
  E_Int nc = cn->getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Compute number of faces per element, nfpe
  std::vector<E_Int> nfpe;
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, true);
  if (ierr != 0) return NULL;
  
  // Compute the number of faces per connectivity
  std::vector<E_Int> nfpc(nc);
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
    E_Int nelts = cm.getSize();
    nfpc[ic] = nelts*nfpe[ic];
  }

  // Compute starting offsets in the face list for each connectivity
  std::vector<E_Int> offsets;
  K_CONNECT::computeStartOffsets(fplp, size, nfpc, offsets);
  nfpc.insert(nfpc.begin(), 0);
  // std::cout << "npts = " << npts << std::endl;
  // std::cout << "size = " << size << std::endl;
  // for (size_t i = 0; i < nfpc.size(); i++)
  //   std::cout << "nfpc = " << nfpc[i] << std::endl;
  // for (size_t i = 0; i < offsets.size(); i++)
  //   std::cout << "offsets = " << offsets[i] << std::endl;

  // Tag vertex indices that are visited when looping over the input facets
  std::vector<E_Int> vindir(npts, 0);

  #pragma omp parallel
  {
    for (E_Int ic = 0; ic < nc; ic++)
    {
      E_Int nvpf, fidx, lcfidx, eidx, f, ind;
      E_Int off = offsets[ic];
      E_Int nfaces = offsets[ic+1] - off;
      K_FLD::FldArrayI& cm = *(cn->getConnect(ic));

      std::vector<std::vector<E_Int> > facets;
      K_CONNECT::getEVFacets(facets, eltTypes[ic], false);
    
      // Loop over each input facet belonging to this connectivity and tag
      // vertex indices of this facet
      #pragma omp for
      for (E_Int i = 0; i < nfaces; i++)
      {
        fidx = fplp[i+off] - 1;      // global face index
        lcfidx = fidx - nfpc[ic];    // face index local to this connectivity
        eidx = lcfidx/nfpe[ic];      // element index local to this connectivity
        f = lcfidx - eidx*nfpe[ic];  // face index local to the element eidx
        nvpf = facets[f].size();     // number of vertices of this facet
        // std::cout << "i = " << i << std::endl;
        // std::cout << "fidx = " << fidx << std::endl;
        // std::cout << "lcfidx = " << lcfidx << std::endl;
        // std::cout << "eidx = " << eidx << std::endl;
        // std::cout << "f = " << f << std::endl;
        // std::cout << "nvpf = " << nvpf << std::endl;
        for (E_Int j = 0; j < nvpf; j++)
        {
          ind = cm(eidx, facets[f][j]) - 1;
          vindir[ind] = 1;
        }
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  // Transform the vertex mask of zeros and ones into a vertex map
  // from old to new connectivities, and get the number of unique visited
  // vertices, npts2
  E_Int npts2 = K_CONNECT::prefixSum(vindir);
  // std::cout << "npts2 = " << npts2 << std::endl;
  
  PyObject* tpl = K_NUMPY::buildNumpyArray(npts2, 1, 1);
  E_Int* vertexPL = K_NUMPY::getNumpyPtrI(tpl);

  #pragma omp parallel
  {
    E_Int indv;
    #pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      indv = vindir[i];
      if (indv > 0) vertexPL[indv-1] = i+1;
    }
  }
  return tpl;
}

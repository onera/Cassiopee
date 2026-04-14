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
#include "Connect/connect.h"
#include <unordered_map>

//=============================================================================
// Change a Elts-Vertex connectivity to a Face-Element connectivity.
// IN: eltType: type d'element (TRI, QUAD,...).
// IN: cEV: Elts-Vertex connectivity. For each elt, give vertices index.
// IN: keepDuplicates: whether to keep duplicates of internal faces.
// OUT: cFE: Face-Element connectivity. For a face, give the neighbour elements.
//      If a face has 1 neighbour element only, the second index value is 0.
// cFE does not need to be preallocated to the total number of faces or number
// of unique faces.
// Return 0 (success), 1 (fail)
// ! Marche pour les maillages non structures BE et ME
// ! Marche uniquement pour les maillage conformes
//=============================================================================
E_Int K_CONNECT::connectEV2FE(
  const char* eltType, FldArrayI& cEV, FldArrayI& cFE, E_Bool keepDuplicates
)
{
  const E_Int UNSET = 0;
  
  E_Int nc = cEV.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Compute number of faces per element, nfpe
  std::vector<E_Int> nfpe;
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, true);
  if (ierr != 0) return 1;

  // Compute cumulative number of elements and faces per connectivity
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;
  std::vector<E_Int> cumnfpc(nc+1); cumnfpc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cEV.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int nfaces = nelts*nfpe[ic];
    cumnepc[ic+1] = cumnepc[ic] + nelts;
    cumnfpc[ic+1] = cumnfpc[ic] + nfaces;
  }
  E_Int ntotFaces = cumnfpc[nc];

  // Hash all faces
  E_Int nvpf, fidx, eidx;
  E_Int face[5];
  TopologyOpt F;
  struct FaceAttrs
  {
    E_Int fidx;  // global face index
    E_Int left;  // left element
    E_Int right;  // right element (-1 if not set yet)
  };
  std::unordered_map<TopologyOpt, FaceAttrs, BernsteinHash<TopologyOpt> > faceMap;
  faceMap.reserve(ntotFaces);

  if (keepDuplicates)
  {
    // Pre-allocate the Face-Element connectivity
    cFE.malloc(ntotFaces, 2);
    E_Int* facesp1 = cFE.begin(1);
    E_Int* facesp2 = cFE.begin(2);

    // Loop over all connectivities
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cEV.getConnect(ic));
      E_Int nelts = cm.getSize();
      std::vector<std::vector<E_Int> > facets;
      K_CONNECT::getEVFacets(facets, eltTypes[ic], false, true);

      for (E_Int i = 0; i < nelts; i++)
      {
        eidx = cumnepc[ic] + i;  // global element index
        // Loop over each facet of this element
        for (E_Int f = 0; f < nfpe[ic]; f++)
        {
          fidx = cumnfpc[ic] + i*nfpe[ic] + f;  // global face index
          nvpf = facets[f].size();  // number of vertices per facet
          // Fill face and insert in map
          for (E_Int j = 0; j < nvpf; j++) face[j] = cm(i, facets[f][j]);
          F.set(face, nvpf);
          auto resf = faceMap.emplace(F, FaceAttrs{fidx, eidx, UNSET});
          if (resf.second)  // first time this face is visited
          {
            facesp1[fidx] = eidx;
            facesp2[fidx] = UNSET;
          }
          else  // this face has already been visited
          {
            const FaceAttrs& attrs = resf.first->second;
            facesp2[attrs.fidx] = eidx;  // amend right element of the orignal face
            // register the duplicated face with swapped left and right elements
            facesp1[fidx] = eidx;
            facesp2[fidx] = attrs.left;
          }
        }
      }
    }
  }
  else
  {
    // Loop over all connectivities
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cEV.getConnect(ic));
      E_Int nelts = cm.getSize();
      std::vector<std::vector<E_Int> > facets;
      K_CONNECT::getEVFacets(facets, eltTypes[ic], false, true);

      for (E_Int i = 0; i < nelts; i++)
      {
        eidx = cumnepc[ic] + i;  // global element index
        // Loop over each facet of this element
        for (E_Int f = 0; f < nfpe[ic]; f++)
        {
          fidx = cumnfpc[ic] + i*nfpe[ic] + f;  // global face index
          nvpf = facets[f].size();  // number of vertices per facet
          // Fill face and insert in map
          for (E_Int j = 0; j < nvpf; j++) face[j] = cm(i, facets[f][j]);
          F.set(face, nvpf);
          auto resf = faceMap.emplace(F, FaceAttrs{fidx, eidx, UNSET});
          if (!resf.second)  // this face has already been visited
          {
            FaceAttrs& attrs = resf.first->second;
            attrs.right = eidx;
          }
        }
      }
    }
    E_Int nuniqueFaces = faceMap.size();

    // Fill the Face-Element connectivity
    cFE.malloc(nuniqueFaces, 2);
    E_Int* facesp1 = cFE.begin(1);
    E_Int* facesp2 = cFE.begin(2);

    for (const auto& facet : faceMap)
    {
      const FaceAttrs& attrs = facet.second;
      fidx = attrs.fidx;
      facesp1[fidx] = attrs.left;
      facesp2[fidx] = attrs.right;
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  return 0;
}

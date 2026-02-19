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

#include "kcore.h"
#include "Connect/connect.h"
#include "parallel.h"
#include <algorithm>
#include <string.h>
#include "String/kstring.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
// Change a Elts-Vertex connectivity to a Element-Element neighbours 
// connectivity.
// IN: cEV: Elts-Vertex connectivity. For each elt, give vertices index.
// IN: nv: nombre de noeuds dans le maillage
// IN: eltType: type d'element (TRI, QUAD,...)
// OUT: cEEN: Elts-Element neighbours connectivity. For each element, 
// give the element index of its neighbours.
// cEEN doit deja etre alloue au nombre d'elements.
// Retourne 0 (echec), 1 (succes)
// ! Marche pour les maillages non structures BE et ME
// ! Marche uniquement pour les maillage conformes
// ! Warning: par definition, un voisin a une facette commune avec un elt
//   Algo: for each mesh element, fill a list of neighbour element candidates
//         gotten from the cVE connectivity (incl. diagonal ones and self).
//         Sort candidates and count number of occurences of each. If greater
//         than 3 in 3D, 2 in 2D, 1 in 1D, then self and candidate share a face.
//=============================================================================
E_Int K_CONNECT::connectEV2EENbrs(
  const char* eltType, E_Int nv, 
  FldArrayI& cEV,
  vector<vector<E_Int> >& cEEN
)
{
  E_Int nc = cEV.getNConnect();

  // Get dimensionality
  E_Int dim = K_CONNECT::getDimME(eltType);
  
  // Compute cumulative number of elements per connectivity (offsets)
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cEV.getConnect(ic));
    E_Int nelts = cm.getSize();
    cumnepc[ic+1] = cumnepc[ic] + nelts;
  }

  // Get vertex -> element connectivity
  vector<vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(cEV, cVE);

  #pragma omp parallel default(shared)
  {
    E_Int nelts, nvpe, eidx, n, ind, prev, noccs, nneis;
    std::vector<E_Int> candidates; candidates.reserve(64);
 
    // Loop over all connectivities to fill cEEN
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cEV.getConnect(ic));
      nelts = cm.getSize();
      nvpe = cm.getNfld();

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        // Fill vector of neighbour element candidates
        // (incl. diagonal ones and self)
        candidates.clear();
        eidx = cumnepc[ic] + i;
        vector<E_Int>& cEENi = cEEN[eidx]; cEENi.reserve(8);
        
        for (E_Int j = 1; j <= nvpe; j++)
        {
          n = cm(i, j);
          const auto& cVEn = cVE[n-1];
          candidates.insert(candidates.end(), cVEn.begin(), cVEn.end());
        }

        // Remove eidx from candidates and sort
        candidates.erase(
          std::remove(candidates.begin(), candidates.end(), eidx),
          candidates.end()
        );
        std::sort(candidates.begin(), candidates.end());

        // Count occurrences and append if >= dim
        prev = -1; noccs = 0;
        nneis = candidates.size();
        for (E_Int e = 0; e < nneis; e++)
        {
          ind = candidates[e];
          if (ind == prev) noccs++;
          else
          {
            if (noccs >= dim) cEENi.push_back(prev);
            prev = ind;
            noccs = 1;
          }
        }

        // Check last candidate
        if (noccs >= dim) cEENi.push_back(prev);
      }
    }
  }

  return 1;  // success
}

//=============================================================================
// Identique mais retourne aussi le no local de la face commune
// Algo: pour chaque facette d'un element, on cherche a la matcher a un voisin.
// ! Marche pour les maillages non structures BE uniquement
//=============================================================================
E_Int K_CONNECT::connectEV2EENbrs(
  const char* eltType, E_Int nv, 
  FldArrayI& cEV,
  vector<vector<E_Int> >& cEEN,
  vector<vector<E_Int> >& commonFace
)
{
  // Number of face per element for each connectivity
  vector<E_Int> nfpe;
  E_Int ierr = getNFPE(nfpe, eltType, true);
  if (ierr != 0) return ierr;

  E_Int nc = cEV.getNConnect();
  vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Compute cumulative number of elements per connectivity (offsets)
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cEV.getConnect(ic));
    E_Int nelts = cm.getSize();
    cumnepc[ic+1] = cumnepc[ic] + nelts;
  }

  // Get vertex -> element connectivity
  vector<vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(cEV, cVE);

  // Boucle sur les connectivites pour remplir cEEN
  #pragma omp parallel default(shared) reduction(+:ierr)
  {
    E_Int nmatch, ind0, ind1, ind2, eidx, n, nidx, nneis, nvpf, nelts, nvpe;
    vector<vector<E_Int> > facets;
 
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cEV.getConnect(ic));
      nelts = cm.getSize();
      nvpe = cm.getNfld();
      ierr += getEVFacets(facets, eltTypes[ic], true);

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        eidx = cumnepc[ic] + i;
        vector<E_Int>& cEEN1 = cEEN[eidx];
        cEEN1.reserve(nfpe[ic]);
        vector<E_Int>& commonFace1 = commonFace[eidx];
        commonFace1.reserve(nfpe[ic]);

        // Loop over each facet of element eidx
        for (E_Int f = 0; f < nfpe[ic]; f++)
        {
          // Number of vertices per face
          nvpf = facets[f].size();
          // First vertex of that facet
          ind0 = cm(i, facets[f][0]) - 1;

          // Get all element indices sharing that vertex
          const vector<E_Int>& cVE1 = cVE[ind0];
          nneis = cVE1.size();

          // Loop over all element sharing that vertex to determine the one 
          // with which this facet is shared
          for (E_Int v = 0; v < nneis; v++)
          {
            nidx = cVE1[v];
            // Skip elements belonging to another connectivity (shortcoming)
            if (nidx < cumnepc[ic] || nidx >= cumnepc[ic] + nelts) continue;
            if (nidx == eidx) continue;
            n = nidx - cumnepc[ic];  // neighbour element index local to this conn.
            nmatch = 0;
            for (E_Int k = 0; k < nvpf; k++)
            {
              ind1 = cm(i, facets[f][k]);
              for (E_Int j = 1; j <= nvpe; j++)
              {
                ind2 = cm(n, j);
                if (ind1 == ind2) { nmatch++; break; }
              }
            }
            if (nmatch == nvpf)
            {
              cEEN1.push_back(nidx); commonFace1.push_back(f);
              break;
            }
          }
        }

        std::sort(cEEN1.begin(), cEEN1.end());
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  if (ierr == 0) ierr = 1; // duct tape because success should be 0, not 1
  else ierr = 0;
  return ierr;
}

//=============================================================================
// Same algo as in connectEV2EENbrs but only returns the number of neighbours
//=============================================================================
E_Int K_CONNECT::connectEV2NNbrs(
  const char* eltType, E_Int nv, 
  FldArrayI& cEV,
  vector<E_Int>& cENN
)
{
  E_Int nc = cEV.getNConnect();

  // Get dimensionality
  E_Int dim = K_CONNECT::getDimME(eltType);

  // Compute cumulative number of elements per connectivity (offsets)
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cEV.getConnect(ic));
    E_Int nelts = cm.getSize();
    cumnepc[ic+1] = cumnepc[ic] + nelts;
  }

  // Get vertex -> element connectivity
  vector<vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(cEV, cVE);

  #pragma omp parallel default(shared)
  {
    E_Int nelts, nvpe, eidx, n, ind, prev, noccs, nneis, count;
    std::vector<E_Int> candidates; candidates.reserve(64);
 
    // Loop over all connectivities to fill cENN
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cEV.getConnect(ic));
      nelts = cm.getSize();
      nvpe = cm.getNfld();

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        // Fill vector of neighbour element candidates
        // (incl. diagonal ones and self)
        candidates.clear();
        eidx = cumnepc[ic] + i;
        
        for (E_Int j = 1; j <= nvpe; j++)
        {
          n = cm(i, j);
          const auto& cVEn = cVE[n-1];
          candidates.insert(candidates.end(), cVEn.begin(), cVEn.end());
        }

        // Remove eidx from candidates
        candidates.erase(
          std::remove(candidates.begin(), candidates.end(), eidx),
          candidates.end()
        );
        std::sort(candidates.begin(), candidates.end());

        // Count occurrences
        count = 0; prev = -1; noccs = 0;
        nneis = candidates.size();
        for (E_Int e = 0; e < nneis; e++)
        {
          ind = candidates[e];
          if (ind == prev) noccs++;
          else
          {
            if (noccs >= dim) count++;
            prev = ind;
            noccs = 1;
          }
        }

        // Catch last candidate
        if (noccs >= dim) count++;
        cENN[eidx] = count;
      }
    }
  }

  return 1;  // success
}

//=============================================================================
/* Recherche de l'element voisin (PENTA ou HEXA) contenant aussi la facette 
   QUAD ABCD dans cEEN */
//=============================================================================
E_Int K_CONNECT::getNbrForQuadFace(
  E_Int indA, E_Int indB, E_Int indC, E_Int indD, 
  FldArrayI& cn, vector<E_Int>& cEEN)
{
  E_Int nvoisins = cEEN.size();
  E_Int ind;
  E_Int nfld = cn.getNfld();
  for (E_Int noet = 0; noet < nvoisins; noet++)
  {
    E_Int et = cEEN[noet];
    E_Int c = 0;
    for (E_Int i = 1; i <= nfld; i++)
    {
      ind = cn(et,i)-1;
      if (ind == indA || ind == indB || ind == indC || ind == indD) c++; 
    }
    if (c == 4) return et;
  }
  return -1;
}

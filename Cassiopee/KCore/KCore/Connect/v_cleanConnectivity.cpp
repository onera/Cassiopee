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
#include <iostream>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include "Array/Array.h"
#include <string.h>
#include "kcore.h"
#include "String/kstring.h"
#include "Connect/connect.h"

#include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/merge.h"

using namespace K_FLD;

const E_Int COLLAPSED = std::numeric_limits<E_Int>::min();

// Nettoyage de connectivites non-structurees
PyObject* K_CONNECT::V_cleanConnectivity(
  const char* varString, FldArrayF& f,
  FldArrayI& cn, const char* eltType,
  E_Float tol, E_Bool rmOverlappingPts,
  E_Bool rmOrphanPts, E_Bool rmDuplicatedFaces,
  E_Bool rmDuplicatedElts, E_Bool rmDegeneratedFaces,
  E_Bool rmDegeneratedElts, E_Bool exportIndirPts
)
{
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "cleanConnectivity: coords must be present in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  PyObject* o = NULL;
  if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0)
  {
    o = V_cleanConnectivityNGon(posx, posy, posz, varString, f, cn,
                                tol, rmOverlappingPts, rmOrphanPts,
                                rmDuplicatedFaces, rmDuplicatedElts,
                                rmDegeneratedFaces, rmDegeneratedElts,
                                exportIndirPts);
  }
  else
  {
    o = V_cleanConnectivityME(posx, posy, posz, varString, f, cn,
                              eltType, tol, rmOverlappingPts, rmOrphanPts,
                              rmDuplicatedElts, rmDegeneratedElts,
                              exportIndirPts);
  }
  return o;
}

// Nettoyage de la connectivite NGON
PyObject* K_CONNECT::V_cleanConnectivityNGon(
  E_Int posx, E_Int posy, E_Int posz, const char* varString,
  FldArrayF& f, FldArrayI& cn,
  E_Float tol, E_Bool rmOverlappingPts, E_Bool rmOrphanPts,
  E_Bool rmDuplicatedFaces, E_Bool rmDuplicatedElts,
  E_Bool rmDegeneratedFaces, E_Bool rmDegeneratedElts,
  E_Bool exportIndirPts
)
{
  E_Bool rmDirtyFaces = (rmDuplicatedFaces || rmDegeneratedFaces);
  E_Bool rmDirtyElts = (rmDuplicatedElts || rmDegeneratedElts);
  
  PyObject* tpl = NULL;
  E_Int *ngon = cn.getNGon(), *indPG = cn.getIndPG();
  E_Int *nface = cn.getNFace(), *indPH = cn.getIndPH();
  E_Int sizeFN = cn.getSizeNGon(), sizeEF = cn.getSizeNFace();
  E_Int nfaces = cn.getNFaces(), nelts = cn.getNElts();
  E_Int nfld = f.getNfld(), npts = f.getSize(), api = f.getApi();
  E_Int ngonType = cn.getNGonType();
  E_Bool hasCnOffsets = false;
  if (ngonType == 2 || ngonType == 3) hasCnOffsets = true;
  E_Int shift = 1; if (ngonType == 3) shift = 0;

  // Get dimensionality
  E_Int dim = 3;
  E_Int nvf0; cn.getFace(0, nvf0, ngon, indPG);
  if (nvf0 == 1) dim = 1;
  else if (nvf0 == 2) dim = 2;
  //std::cout<<"dim: " << dim << std::endl;
  
  // --- 1. Points ---
  // 1a. Identify orphan points, ie, initialise indirection table for use in 1b
  E_Int nuniquePts = npts;
  std::vector<E_Int> indir;
  if (dim == 0) rmOrphanPts = false;
  if (rmOrphanPts)
  {
    E_Int nv, vidx;
    indir.resize(npts, -1);
    for (E_Int i = 0; i < nfaces; i++)
    {
      E_Int* face = cn.getFace(i, nv, ngon, indPG);
      for (E_Int p = 0; p < nv; p++)
      {
        vidx = face[p]-1;
        indir[vidx] = vidx;
      }
    }

    rmOrphanPts = false;
    for (E_Int i = 0; i < npts; i++)
    {
      if (indir[i] == -1)
      {
        rmOrphanPts = true;
        break;
      }
    }
  }

  // 1b. Identify overlapping points geometrically
  E_Bool rmDirtyPts = (rmOverlappingPts || rmOrphanPts);
  if (rmOverlappingPts)
  {
    nuniquePts = K_CONNECT::V_identifyDirtyPoints(posx, posy, posz, f, tol,
                                                  indir, rmOverlappingPts);
    if (nuniquePts < 0) return NULL;
    else if (nuniquePts == npts) rmDirtyPts = false;
  }
  //std::cout<<"npts: " << npts << std::endl;
  //std::cout<<"nuniquePts: " << nuniquePts << std::endl;

  // An update is necessary before topological operations in 2. & 3.
  if (rmDirtyPts)
  {
    E_Int j, itrl, nv, vidx;
    E_Int ind = 0;
    // 1.c Reindex vertices in FN (no change in size)
    for (E_Int i = 0; i < nfaces; i++)
    {
      cn.getFace(i, nv, ngon, indPG);
      for (E_Int p = shift; p < nv+shift; p++)
      {
        vidx = ngon[ind+p]-1;
        if (indir[vidx] != vidx) ngon[ind+p] = indir[vidx]+1;
      }
      ind += nv+shift;
    }

    // 1.d Reindex and compress fields
    for (E_Int fld = 1; fld <= nfld; fld++)
    {
      j = 0, itrl = -1;
      E_Float* fp = f.begin(fld);
      for (E_Int i = 0; i < npts; i++)
      {
        // Write if point is neither a duplicate nor an orphan
        if (indir[i] < 0) continue;
        if (itrl == -1 || indir[i] - indir[itrl] == 1)
        {
          fp[j] = fp[i];
          j += 1; itrl = i;
        }
      }
    }
  }
  
  if (!exportIndirPts) indir.clear();

  // --- 2. Identify dirty elements topologically ---
  E_Int nuniqueElts = nelts;
  std::vector<E_Int> indirPH;
  if (rmDirtyElts)
  {
    nuniqueElts = K_CONNECT::V_identifyDirtyElementsNGon(dim, cn, nface, indPH,
                                                         indirPH,
                                                         rmDegeneratedElts);
    if (nuniqueElts != nelts) rmDuplicatedFaces = true;
    else rmDuplicatedElts = false;
  }
  //std::cout<<"nelts: " << nelts << std::endl;
  //std::cout<<"nuniqueElts: " << nuniqueElts << std::endl;

  // --- 3. Identify dirty faces topologically ---
  E_Int nuniqueFaces = nfaces;
  std::vector<E_Int> indirPG;
  if (rmDirtyFaces)
  {
    nuniqueFaces = K_CONNECT::V_identifyDirtyFacesNGon(dim, cn, ngon, indPG,
                                                       indirPG,
                                                       rmDegeneratedFaces);
    if (nuniqueFaces == nfaces) rmDuplicatedFaces = false;
  }
  //std::cout<<"nfaces: " << nfaces << std::endl;
  //std::cout<<"nuniqueFaces: " << nuniqueFaces << std::endl;

  // --- 4. Reindex & Compress connectivities ---
  E_Int j, k; // write pointers
  E_Int ind; // read pointer, same for i
  E_Int itrl, nv, nvins, vidx, nf, nfins, fidx, indirPHi;
  E_Int sizeFN2 = sizeFN, sizeEF2 = sizeEF;

  // 4.a Compress FN connectivity
  if (rmDirtyFaces)
  {
    j = 0; itrl = -1; k = 0; ind = 0;
    for (E_Int i = 0; i < nfaces; i++)
    {
      cn.getFace(i, nv, ngon, indPG);
      
      // Detect duplicated or collapsed faces
      if ((itrl != -1 && abs(indirPG[i]) - abs(indirPG[itrl]) != 1) || indirPG[i] == COLLAPSED)
      {
        // Shift read pointer and skip
        ind += nv+shift;
        continue;
      }

      // Compress
      if (indirPG[i] >= 0) // non-degenerated face
      {
        nvins = nv;
        for (E_Int p = 0; p < nv+shift; p++) ngon[k+p] = ngon[ind+p];
      }
      else // degenerated face
      {
        nvins = 0; // actual number of vertices inserted
        std::unordered_set<E_Int> vertexSet;
        for (E_Int p = shift; p < nv+shift; p++)
        {
          vidx = ngon[ind+p];
          // If the index is not in the set, add it and update FN
          if (vertexSet.insert(vidx).second)
          {
            ngon[k+shift+nvins] = vidx; nvins++;
          }
        }
        if (shift == 1) ngon[k] = nvins;
      }
      if (hasCnOffsets) // array2 or array3
      {
        indPG[j] = k; j += 1;
      }
      
      // Shift read and write pointers
      ind += nv+shift; k += nvins+shift; itrl = i;
    }
    sizeFN2 = k;
  }
  
  // 4.b EF connectivity
  if (rmDirtyFaces || rmDirtyElts)
  {
    j = 0; itrl = -1; k = 0; ind = 0;
    for (E_Int i = 0; i < nelts; i++)
    {
      cn.getElt(i, nf, nface, indPH);
      if (rmDirtyElts) indirPHi = indirPH[i];
      else indirPHi = i;

      // Detect duplicated or collapsed elements
      if ((indirPHi == COLLAPSED) ||
          (rmDirtyElts && itrl != -1 && indirPHi - indirPH[itrl] != 1))
      {
        // Shift read pointer and skip
        ind += nf+shift;
        continue;
      }
      
      // Reindex and /or compress
      if (indirPHi >= 0 && !rmDegeneratedFaces) // non-degenerated face
      {
        nfins = nf;
        nface[k] = nface[ind];
        for (E_Int f = shift; f < nf+shift; f++)
        {
          fidx = nface[ind+f];
          if (rmDuplicatedFaces) fidx = indirPG[fidx-1];
          nface[k+f] = fidx;
        }
      }
      else
      {
        nfins = 0;
        for (E_Int f = shift; f < nf+shift; f++)
        {
          fidx = nface[ind+f];
          if (rmDuplicatedFaces) fidx = indirPG[fidx-1];
          if (fidx == COLLAPSED) continue;
          nface[k+shift+nfins] = abs(fidx); nfins++;
        }
        if (shift == 1) nface[k] = nfins;
      }
      
      if (rmDirtyElts && hasCnOffsets) // array2 or array3
      {
        indPH[j] = k; j += 1;
      }
      
      // Shift read and write pointers
      ind += nf+shift; k += nfins+shift; itrl = i;
    }
    sizeEF2 = k;

    //std::cout<<"sizeFN2: " << sizeFN2 << std::endl;
    //std::cout<<"sizeEF2: " << sizeEF2 << std::endl;
  }

  // --- 5. Create resized connectivity ---
  if (rmOverlappingPts || rmDuplicatedFaces || rmDuplicatedElts)
  {  
    E_Bool center = false;
    tpl = K_ARRAY::buildArray3(nfld, varString, nuniquePts, nuniqueElts,
                               nuniqueFaces, "NGON", sizeFN2, sizeEF2,
                               ngonType, center, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    E_Int *ngon2 = cn2->getNGon(), *nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (hasCnOffsets)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }

    #pragma omp parallel
    {
      // Copy compressed fields to f2
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < nuniquePts; i++) f2p[i] = fp[i];
      }

      // Copy compressed connectivity to cn2
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeFN2; i++) ngon2[i] = ngon[i];
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeEF2; i++) nface2[i] = nface[i];

      if (hasCnOffsets) // array2 or array3
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < nuniqueFaces; i++) indPG2[i] = indPG[i];
        #pragma omp for
        for (E_Int i = 0; i < nuniqueElts; i++) indPH2[i] = indPH[i];
      }
    }

    RELEASESHAREDU(tpl, f2, cn2);
  }
  
  if (exportIndirPts)
  {
    PyObject* vmap = K_NUMPY::buildNumpyArray(npts, 1, 1);
    E_Int* vmapp = K_NUMPY::getNumpyPtrI(vmap);
    #pragma omp parallel for
    for (E_Int i = 0; i < npts; i++) vmapp[i] = indir[i];
    return Py_BuildValue("(OO)", tpl, vmap);
  }
  return tpl;
}

// Nettoyage de la connectivite et du tableau des sommets.
E_Int K_CONNECT::V_identifyDirtyPoints(
  E_Int posx, E_Int posy, E_Int posz, 
  FldArrayF& f, E_Float tol, std::vector<E_Int>& indir,
  E_Bool rmOverlappingPts
)
{
  size_t npts = f.getSize();
  E_Bool rmOrphanPts = (indir.size() == npts);
  if (!(rmOverlappingPts || rmOrphanPts)) return npts; // no orphans nor duplicates
  
  // Local indirection table
  std::vector<E_Int> locIndir(npts);
  
  if (rmOverlappingPts)
  {
    if (posx == -1 || posy == -1 || posz == -1)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "cleanConnectivity: field must have coordinates.");
      return -1;
    }

    // Build the KdTree
    ArrayAccessor<FldArrayF> coordAcc(f, posx, posy, posz);
    #pragma omp parallel for
    for (size_t i = 0; i < npts; ++i) locIndir[i] = i;
    K_SEARCH::KdTree<FldArrayF> moving_tree(coordAcc, locIndir, tol, true);
    
    std::vector<triplet_t> palma;
    std::vector<std::vector<triplet_t> > palma_thrd(__NUMTHREADS__);
    std::vector<std::vector<E_Int> > onodes_thrd(__NUMTHREADS__);
    std::vector<std::vector<E_Float> > dist2_thrd(__NUMTHREADS__);

    // Loop on target nodes and get all moving nodes in sphere of radius tol.
    #pragma omp parallel
    {
      size_t osz;
      E_Int id, fj;
      E_Float d2;

      id = __CURRENT_THREAD__;

      #pragma omp for schedule(static)
      for (size_t fi = 0; fi < npts; fi++)
      {
        onodes_thrd[id].clear(); dist2_thrd[id].clear();

        moving_tree.getInSphere(fi, tol, onodes_thrd[id], dist2_thrd[id]);
        osz = onodes_thrd[id].size();

        for (size_t j = 0; j < osz; ++j)
        {
          d2 = dist2_thrd[id][j];
          fj = onodes_thrd[id][j];
          palma_thrd[id].push_back(triplet_t(d2, fi, fj));
        }
      }
    }

    size_t sz = 0;
    for (E_Int t = 0; t < __NUMTHREADS__; ++t) sz += palma_thrd[t].size();

    palma.reserve(sz);

    for (E_Int t = 0; t < __NUMTHREADS__; ++t)
    {
      if (!palma_thrd[t].empty())
        palma.insert(palma.end(), ALL(palma_thrd[t]));
    }

    // Move nodes
    rmOverlappingPts = !palma.empty(); // update switch
    if (rmOverlappingPts) // duplicates found
    {
      // Sort tuples by increasing distance
      std::stable_sort(palma.begin(), palma.end(), lower_than);

      E_Int fi, fj;
      size_t psz = palma.size();
      for (size_t i = 0; i < psz; ++i)
      {
        fi = palma[i].Ni;
        fj = palma[i].Nj;

        // Merging rule: if target and moving have not been moved already
        if ((locIndir[fi] == fi) && (locIndir[fj] == fj)) locIndir[fj] = fi;
      }
      
      // Remove circular references by instructing the pointers to point to
      // the roots
      if (rmOrphanPts)
      {
        E_Int prev, fi;
        // seen array used to detect cycles
        E_Int* seen = (E_Int*) calloc(npts, sizeof(E_Int));
        E_Int stamp = 1;
    
        for (size_t i = 0; i < npts; ++i)
        {
          // Skip orphan points
          if (indir[i] == -1) { locIndir[i] = -1; continue; }

          prev = i;
          fi = locIndir[i];

          while (true)
          {
            // If orphan node, cut at last valid ancestor prev
            if (fi == -1 or indir[fi] == -1)
            {
              locIndir[i] = prev;
              break;
            }
            // If fi is a valid self-pointer, cut at fi
            else if (locIndir[fi] == fi)
            {
              locIndir[i] = fi;
              break;
            }
            // If fi is revisited in this traversal, break the cycle
            else if (seen[fi] == stamp)
            {
              locIndir[i] = prev;
              break;
            }

            seen[fi] = stamp;
            prev = fi;
            fi = locIndir[fi];
          }
          stamp++;
        }
        free(seen);
      }
      else
      {
        // Path compression step from a disjoint set data structure
        for (size_t i = 0; i < npts; ++i)
        {
          fi = locIndir[i];
          while (fi != locIndir[fi]) fi = locIndir[fi];
          locIndir[i] = fi;
        }
      }
    }
  }

  if(!rmOverlappingPts and rmOrphanPts)
  {
    // Orphans only as overlapping points were either not searched for or none
    // were found
    #pragma omp parallel for
    for (size_t i = 0; i < npts; ++i) locIndir[i] = indir[i];
  }

  // Vertex set to discard orphans and/or duplicates
  indir.clear(); indir.resize(npts);
  std::unordered_set<E_Int> vertexSet;

  // Loop over all vertices
  E_Int nuniquePts = 0;
  for (size_t i = 0; i < npts; ++i)
  {
    if (locIndir[i] == -1)
    {
      // Skip orphan points
      indir[i] = -1;
    }
    else if (vertexSet.insert(locIndir[i]).second)
    {
      // First time this vertex is encountered, added to the set
      indir[i] = nuniquePts; nuniquePts++;
    }
    else indir[i] = indir[locIndir[i]];
  }

  return vertexSet.size();
}

// Identify dirty faces topologically
E_Int K_CONNECT::V_identifyDirtyFacesNGon(
  E_Int dim, FldArrayF& f, FldArrayI& cn, std::vector<E_Int>& indirPG,
  E_Bool rmDegeneratedFaces
)
{
  E_Int *ngon = cn.getNGon(), *indPG = cn.getIndPG();
  return K_CONNECT::V_identifyDirtyFacesNGon(dim, cn, ngon, indPG, indirPG,
                                             rmDegeneratedFaces);
}

E_Int K_CONNECT::V_identifyDirtyFacesNGon(
  E_Int dim, FldArrayI& cn, E_Int* ngon, E_Int* indPG,
  std::vector<E_Int>& indirPG, E_Bool rmDegeneratedFaces
)
{
  E_Int nv, nuniqueFaces = 1;
  E_Int nfaces = cn.getNFaces();
  indirPG.clear(); indirPG.resize(nfaces);

  // Face map to eliminate duplicates and collapsed faces
  Topology F;
  std::unordered_map<Topology, E_Int, JenkinsHash<Topology> > faceMap;
  
  if (rmDegeneratedFaces && dim > 1)
  {
    // Loop over all faces of the NGON connectivity
    for (E_Int i = 0; i < nfaces; i++)
    {
      E_Int* face = cn.getFace(i, nv, ngon, indPG);
      F.set(face, nv, rmDegeneratedFaces);
      if (F.size_ < (size_t)dim) // collapsed
      {
        indirPG[i] = COLLAPSED;
      }
      else // degenerated or okay
      {
        // Use insert to ensure F is initially mapped to -1 if it doesn't exist
        auto res = faceMap.insert(std::make_pair(F, -1));
        // Check the value associated with F. If it is -1, then first time this
        // face is encountered
        if (res.first->second == -1)
        {
          res.first->second = nuniqueFaces; nuniqueFaces++;
        }
        indirPG[i] = res.first->second;
        // Flip sign of degenerated faces
        if (F.isDegen_) indirPG[i] *= -1;
      }
    }
  }
  else
  {
    for (E_Int i = 0; i < nfaces; i++)
    {
      E_Int* face = cn.getFace(i, nv, ngon, indPG);
      F.set(face, nv);
      // Use insert to ensure F is initially mapped to -1 if it doesn't exist
      auto res = faceMap.insert(std::make_pair(F, -1));
      // Check the value associated with F. If it is -1, then first time this
      // face is encountered
      if (res.first->second == -1)
      {
        res.first->second = nuniqueFaces; nuniqueFaces++;
      }
      indirPG[i] = res.first->second;
    }
  }
  return nuniqueFaces-1;
}

// Identify dirty elements topologically 
E_Int K_CONNECT::V_identifyDirtyElementsNGon(
  E_Int dim, FldArrayF& f, FldArrayI& cn, std::vector<E_Int>& indirPH,
  E_Bool rmDegeneratedElts
)
{
  E_Int *nface = cn.getNFace(), *indPH = cn.getIndPH();
  return K_CONNECT::V_identifyDirtyElementsNGon(dim, cn, nface, indPH, indirPH,
                                                rmDegeneratedElts);
}

E_Int K_CONNECT::V_identifyDirtyElementsNGon(
  E_Int dim, FldArrayI& cn, E_Int* nface, E_Int* indPH,
  std::vector<E_Int>& indirPH, E_Bool rmDegeneratedElts
)
{
  E_Int nf, nuniqueElts = 1;
  E_Int nelts = cn.getNElts();
  indirPH.clear(); indirPH.resize(nelts);

  // Element map to eliminate duplicates
  Topology E;
  std::unordered_map<Topology, E_Int, JenkinsHash<Topology> > eltMap;

  if (rmDegeneratedElts && dim > 0)
  {
    // Loop over all faces of the NGON connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int* elt = cn.getElt(i, nf, nface, indPH);
      E.set(elt, nf, true);
      if (E.size_ <= (size_t)dim) // collapsed
      {
        indirPH[i] = COLLAPSED;
      }
      else // degenerated or okay
      {
        // Use insert to ensure E is initially mapped to -1 if it doesn't exist
        auto res = eltMap.insert(std::make_pair(E, -1));
        // Check the value associated with E. If it is -1, then first time this
        // element is encountered
        if (res.first->second == -1) 
        {
          res.first->second = nuniqueElts; nuniqueElts++;
        }
        indirPH[i] = res.first->second;
        // Flip sign of degenerated elements
        if (E.isDegen_) indirPH[i] *= -1;
      }
    }
  }
  else
  {
    // Loop over all faces of the NGON connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int* elt = cn.getElt(i, nf, nface, indPH);
      E.set(elt, nf); 
      // Use insert to ensure E is initially mapped to -1 if it doesn't exist
      auto res = eltMap.insert(std::make_pair(E, -1));
      // Check the value associated with E. If it is -1, then first time this
      // element is encountered
      if (res.first->second == -1) 
      {
        res.first->second = nuniqueElts; nuniqueElts++;
      }
      indirPH[i] = res.first->second;
    }
  }
  return nuniqueElts-1;
}

// Nettoyage de la connectivite ME
PyObject* K_CONNECT::V_cleanConnectivityME(
  E_Int posx, E_Int posy, E_Int posz, const char* varString,
  FldArrayF& f, FldArrayI& cn, const char* eltType,
  E_Float tol, E_Bool rmOverlappingPts, E_Bool rmOrphanPts,
  E_Bool rmDuplicatedElts, E_Bool rmDegeneratedElts,
  E_Bool exportIndirPts
)
{
  E_Bool rmDirtyElts = (rmDuplicatedElts || rmDegeneratedElts);
  
  PyObject* tpl = NULL;
  E_Int nc = cn.getNConnect();
  E_Int vidx;
  E_Int nfld = f.getNfld(), npts = f.getSize(), api = f.getApi();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Get dimensionality
  E_Int dim = 3;
  if (strcmp(eltTypes[0], "NODE") == 0) dim = 0;
  else if (strcmp(eltTypes[0], "BAR") == 0) dim = 1;
  else if (strcmp(eltTypes[0], "TRI") == 0 ||
           strcmp(eltTypes[0], "QUAD") == 0) dim = 2;
  for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];
  
  //for (E_Int ic = 0; ic < nc; ic++) std::cout<<"eltType: " << eltTypes[ic] << std::endl;
  //std::cout<<"dim: " << dim << std::endl;

  // Compute total number of elements
  E_Int neltsTot = 0;
  std::vector<E_Int> nepc(nc); // initial number of elements per connectivity
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic] = nelts;
    neltsTot += nelts;
  }
  
  // --- 1. Points ---
  // 1a. Identify orphan points, ie, initialise indirection table used in 1b
  E_Int nuniquePts = 0;
  std::vector<E_Int> indir;
  if (dim == 0) rmOrphanPts = false;
  if (rmOrphanPts)
  {
    indir.resize(npts, -1);
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int nvpe = cm.getNfld();
      for (E_Int i = 0; i < nelts; i++)
      {
        for (E_Int j = 1; j <= nvpe; j++)
        {
          vidx = cm(i,j)-1;
          indir[vidx] = vidx;
        }
      }
    }
    
    rmOrphanPts = false;
    for (E_Int i = 0; i < npts; i++)
    {
      if (indir[i] == -1)
      {
        rmOrphanPts = true;
        break;
      }
    }
  }

  // 1b. Identify overlapping points geometrically
  E_Bool rmDirtyPts = (rmOverlappingPts || rmOrphanPts);
  if (rmDirtyPts)
  {
    nuniquePts = K_CONNECT::V_identifyDirtyPoints(posx, posy, posz, f, tol,
                                                  indir, rmOverlappingPts);
    if (nuniquePts < 0) return NULL;
    else if (nuniquePts == npts) rmDirtyPts = false;
  }
  //std::cout<<"npts: " << npts << std::endl;
  //std::cout<<"nuniquePts: " << nuniquePts << std::endl;
  
  // An update is necessary before topological operations in 2.
  if (rmDirtyPts)
  {
    // 1.c Reindex vertices (no change in size)
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int nvpe = cm.getNfld();
      for (E_Int i = 0; i < nelts; i++)
      {
        for (E_Int j = 1; j <= nvpe; j++)
        {
          vidx = cm(i,j)-1;
          if (indir[vidx] != vidx) cm(i,j) = indir[vidx]+1;
        }
      }
    }

    // 1.d Reindex and compress fields
    for (E_Int fld = 1; fld <= nfld; fld++)
    {
      E_Int j = 0, itrl = -1;
      E_Float* fp = f.begin(fld);
      for (E_Int i = 0; i < npts; i++)
      {
        // Write if point is neither a duplicate nor an orphan
        if (indir[i] < 0) continue;
        if (itrl == -1 || indir[i] - indir[itrl] == 1)
        {
          fp[j] = fp[i];
          j += 1; itrl = i;
        }
      }
    }
    
    if (!exportIndirPts) indir.clear();

    // 1.d Resize fields
    //f.reAllocMat(nuniquePts,nfld);
  }

  // --- 2. Identify dirty elements topologically ---
  E_Int nuniqueEltsTot = neltsTot;
  std::vector<E_Int> indirPH, nuniqueElts;
  if (rmDirtyElts)
  {
    nuniqueEltsTot = K_CONNECT::V_identifyDirtyElementsME(
        dim, cn, indirPH, nuniqueElts, neltsTot, rmDegeneratedElts);
    if (nuniqueEltsTot == neltsTot) rmDuplicatedElts = false;
  }
  else
  {
    nuniqueElts.resize(nc);
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      nuniqueElts[ic] = cm.getSize();
    }
  }
  //std::cout<<"neltsTot: " << neltsTot << std::endl;
  //std::cout<<"nuniqueEltsTot: " << nuniqueEltsTot << std::endl;

  // --- 3. Reindex & Compress connectivities ---
  if (rmDirtyElts)
  {
    E_Int k = 0; // write pointer
    E_Int itrl = -1;
    E_Int elOffset = 0;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      // Skip connectivity if none of its elements are duplicated/degenerated
      if (nuniqueElts[ic] == nepc[ic]) continue;
      
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int nvpe = cm.getNfld();
      for (E_Int i = 0; i < nelts; i++)
      {
        // Skip duplicated/degenerated elements
        if (indirPH[elOffset+i] == COLLAPSED ||
            (itrl != -1 && indirPH[elOffset+i] - indirPH[elOffset+itrl] != 1)) continue;
        
        for (E_Int j = 1; j <= nvpe; j++) cm(k,j) = cm(i,j);
        itrl = i; k += 1;
      }
      
      elOffset += nelts;
    }
  }

  // --- 4. Create resized connectivity ---
  if (rmOverlappingPts || rmDirtyElts)
  {  
    E_Bool center = false;
    tpl = K_ARRAY::buildArray3(nfld, varString, nuniquePts, nuniqueElts,
                               eltType, center, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);

    #pragma omp parallel
    {
      // Copy compressed fields to f2
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for
        for (E_Int i = 0; i < nuniquePts; i++) f2p[i] = fp[i];
      }

      // Copy compressed connectivity to cn2
      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cn.getConnect(ic));
        FldArrayI& cm2 = *(cn2->getConnect(ic));
        E_Int nvpe = cm.getNfld();
        #pragma omp for
        for (E_Int i = 0; i < nuniqueElts[ic]; i++)
          for (E_Int j = 1; j <= nvpe; j++) cm2(i,j) = cm(i,j);
      }
    }
    RELEASESHAREDU(tpl, f2, cn2);
  }
  
  if (exportIndirPts)
  {
    PyObject* vmap = K_NUMPY::buildNumpyArray(npts, 1, 1);
    E_Int* vmapp = K_NUMPY::getNumpyPtrI(vmap);
    #pragma omp parallel for
    for (E_Int i = 0; i < npts; i++) vmapp[i] = indir[i];
    return Py_BuildValue("(OO)", tpl, vmap);
  }
  return tpl;
}

E_Int K_CONNECT::V_identifyDirtyElementsME(
  E_Int dim, FldArrayI& cn, std::vector<E_Int>& indir,
  std::vector<E_Int>& nuniqueElts, E_Int neltsTot,
  E_Bool rmDegeneratedElts
)
{
  E_Int nc = cn.getNConnect();
  E_Int nuniqueEltsTot = 0;
  // Compute total number of elements if not provided
  if (neltsTot == 0)
  {
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      neltsTot += nelts;
    }
  }
  indir.clear(); indir.resize(neltsTot);
  nuniqueElts.clear(); nuniqueElts.resize(nc);

  // Element map to eliminate duplicates
  E_Int elOffset = 0;
  TopologyOpt E;
  std::unordered_map<TopologyOpt, E_Int, JenkinsHash<TopologyOpt> > eltMap;

  if (rmDegeneratedElts && dim > 0)
  {
    std::vector<E_Int> isTotDegen(4, -1);
    for (E_Int i = 0; i <= 3; i++) isTotDegen[i] = i;

    E_Int eltc[9];

    // Loop over ME connectivity
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      nuniqueElts[ic] = 1;
      E_Int nvpe = cm.getNfld();
      //std::vector<E_Int> elt(nvpe);
      for (E_Int i = 0; i < nelts; i++)
      {
        for (E_Int j = 1; j <= nvpe; j++) eltc[j-1] = cm(i,j);
        E.set(eltc, nvpe, true);
        if (E.isDegen_ and (E_Int)E.size_ <= isTotDegen[dim])
        {
          indir[elOffset+i] = COLLAPSED; continue;
        } // TODO no ME support yet: for ex. create TRI conn from QUAD conn
        // Use insert to ensure E is initially mapped to -1 if it doesn't exist
        auto res = eltMap.insert(std::make_pair(E, -1));
        // Check the value associated with E. If it is -1, then first time this
        // element is encountered
        if (res.first->second == -1)
        {
          res.first->second = nuniqueElts[ic]; nuniqueElts[ic]++;
        }
        indir[elOffset+i] = res.first->second;
      }
      elOffset += nelts;
      nuniqueElts[ic] -= 1;
      nuniqueEltsTot += nuniqueElts[ic];
    }
  }
  else
  {
    E_Int eltc[9];

    // Loop over ME connectivity
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      nuniqueElts[ic] = 1;
      E_Int nvpe = cm.getNfld();
      //std::vector<E_Int> elt(nvpe);
      for (E_Int i = 0; i < nelts; i++)
      {
        for (E_Int j = 1; j <= nvpe; j++) eltc[j-1] = cm(i,j);
        E.set(eltc, nvpe);
        // Use insert to ensure E is initially mapped to -1 if it doesn't exist
        auto res = eltMap.insert(std::make_pair(E, -1));
        // Check the value associated with E. If it is -1, then first time this
        // element is encountered
        if (res.first->second == -1)
        {
          res.first->second = nuniqueElts[ic]; nuniqueElts[ic]++;
        }
        indir[elOffset+i] = res.first->second;
      }
      elOffset += nelts;
      nuniqueElts[ic] -= 1;
      nuniqueEltsTot += nuniqueElts[ic];
    }
  }
  return nuniqueEltsTot;
}

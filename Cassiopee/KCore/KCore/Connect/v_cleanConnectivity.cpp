/*    
    Copyright 2013-2024 Onera.

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
#include "String/kstring.h"
#include "Connect/connect.h"

# include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/merge.h"

using namespace K_FLD;

const E_Int COLLAPSED = std::numeric_limits<E_Int>::min();

// Nettoyage de connectivites non-structurees
PyObject* K_CONNECT::V_cleanConnectivity(
  const char* varString, FldArrayF& f,
  FldArrayI& cn, const char* eltType,
  E_Float tol, E_Bool removeOverlappingPoints,
  E_Bool removeOrphanPoints,
  E_Bool removeDuplicatedFaces,
  E_Bool removeDuplicatedElements,
  E_Bool removeDegeneratedFaces,
  E_Bool removeDegeneratedElements
)
{
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "cleanConnectivity: coord must be present in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  PyObject* tpl = NULL;
  if (K_STRING::cmp(eltType, "NGON") == 0 || K_STRING::cmp(eltType, "NGON*") == 0)
  {
    tpl = V_cleanConnectivityNGon(posx, posy, posz, varString, f, cn,
                                  tol, removeOverlappingPoints,
                                  removeOrphanPoints,
                                  removeDuplicatedFaces,
                                  removeDuplicatedElements,
                                  removeDegeneratedFaces,
                                  removeDegeneratedElements);
  }
  else
  {
    tpl = V_cleanConnectivityME(posx, posy, posz, varString, f, cn,
                                eltType, tol, removeOverlappingPoints,
                                removeOrphanPoints,
                                removeDuplicatedElements,
                                removeDegeneratedElements);
  }
  return tpl;
}

// Nettoyage de la connectivite NGON
PyObject* K_CONNECT::V_cleanConnectivityNGon(
  E_Int posx, E_Int posy, E_Int posz, const char* varString,
  FldArrayF& f, FldArrayI& cn,
  E_Float tol, E_Bool removeOverlappingPoints,
  E_Bool removeOrphanPoints,
  E_Bool removeDuplicatedFaces,
  E_Bool removeDuplicatedElements,
  E_Bool removeDegeneratedFaces,
  E_Bool removeDegeneratedElements
)
{
  E_Bool removeDirtyFaces = (removeDuplicatedFaces || removeDegeneratedFaces);
  E_Bool removeDirtyElements = (removeDuplicatedElements || removeDegeneratedElements);
  
  PyObject* tpl = NULL;
  E_Int *ngon = cn.getNGon(), *indPG = cn.getIndPG();
  E_Int *nface = cn.getNFace(), *indPH = cn.getIndPH();
  E_Int sizeFN = cn.getSizeNGon(), sizeEF = cn.getSizeNFace();
  E_Int nfaces = cn.getNFaces(), nelts = cn.getNElts();
  E_Int nfld = f.getNfld(), npts = f.getSize(), api = f.getApi();
  E_Bool array23 = false;
  if (api == 2 || api == 3) array23 = true;
  if (api == 2) api = 3;
  E_Int shift = 1; if (api == 3) shift = 0;

  // Get dimensionality
  E_Int dim = 3;
  E_Int nvf0; cn.getFace(0, nvf0, ngon, indPG);
  if (nvf0 == 1) dim = 1;
  else if (nvf0 == 2) dim = 2;
  std::cout<<"dim: " << dim << std::endl;
  
  // --- 1. Points ---
  // 1a. Identify orphan points, ie, initialise indirection table for use in 1b
  E_Int nuniquePts = npts;
  std::vector<E_Int> indir;
  if (dim == 0) removeOrphanPoints = false;
  if (removeOrphanPoints)
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

    removeOrphanPoints = false;
    for (E_Int i = 0; i < npts; i++)
    {
      if (indir[i] == -1)
      {
        removeOrphanPoints = true;
        indir.clear();
        break;
      }
    }
  }

  // 1b. Identify overlapping points geometrically
  if (removeOverlappingPoints)
  {
    nuniquePts = K_CONNECT::V_identifyOverlappingPoints(posx, posy, posz, f, tol, indir);
    if (nuniquePts < 0) return NULL;
    else if (nuniquePts == npts) removeOverlappingPoints = false;
  }
  std::cout<<"npts: " << npts << std::endl;
  std::cout<<"nuniquePts: " << nuniquePts << std::endl;

  E_Bool removeDirtyPoints = (removeOverlappingPoints || removeOrphanPoints);

  // An update is necessary before topological operations in 2. & 3.
  if (removeDirtyPoints)
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

  // --- 2. Identify dirty elements topologically ---
  E_Int nuniqueElts = nelts;
  std::vector<E_Int> indirPH;
  if (removeDirtyElements)
  {
    nuniqueElts = K_CONNECT::V_identifyDirtyElementsNGon(dim, cn, nface, indPH, indirPH, removeDegeneratedElements);
    if (nuniqueElts != nelts) removeDuplicatedFaces = true;
    else removeDuplicatedElements = false;
  }
  std::cout<<"nelts: " << nelts << std::endl;
  std::cout<<"nuniqueElts: " << nuniqueElts << std::endl;

  // --- 3. Identify dirty faces topologically ---
  E_Int nuniqueFaces = nfaces;
  std::vector<E_Int> indirPG;
  if (removeDirtyFaces)
  {
    nuniqueFaces = K_CONNECT::V_identifyDirtyFacesNGon(dim, cn, ngon, indPG, indirPG, removeDegeneratedFaces);
    if (nuniqueFaces == nfaces) removeDuplicatedFaces = false;
  }
  std::cout<<"nfaces: " << nfaces << std::endl;
  std::cout<<"nuniqueFaces: " << nuniqueFaces << std::endl;

  // --- 4. Reindex & Compress connectivities ---
  E_Int j, k; // write pointers
  E_Int ind; // read pointer, same for i
  E_Int itrl, nv, nvins, vidx, nf, nfins, fidx;
  E_Int indirPHi, indirPHitrl;
  E_Int sizeFN2 = sizeFN, sizeEF2 = sizeEF;

  // 4.a Compress FN connectivity
  if (removeDirtyFaces)
  {
    j = 0; itrl = -1; k = 0; ind = 0;
    for (E_Int i = 0; i < nfaces; i++)
    {
      cn.getFace(i, nv, ngon, indPG);
      //std::cout<<"cmpt: " << i << " " << itrl << " " << k << " " << std::endl;
      
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
      if (array23) // array2 or array3
      {
        indPG[j] = k; j += 1;
      }
      
      // Shift read and write pointers
      ind += nv+shift; k += nvins+shift; itrl = i;
    }
    sizeFN2 = k;
  }
  
  // 4.b EF connectivity
  if (removeDirtyFaces || removeDirtyElements)
  {
    j = 0; itrl = -1; k = 0; ind = 0;
    for (E_Int i = 0; i < nelts; i++)
    {
      cn.getElt(i, nf, nface, indPH);
      if (removeDirtyElements)
      {
        indirPHi = indirPH[i]; indirPHitrl = indirPH[itrl];
      }
      else
      {
        indirPHi = i; indirPHitrl = itrl;
      }

      // Detect duplicated or collapsed elements
      if ((indirPHi == COLLAPSED) ||
          (removeDirtyElements && (itrl != -1 && indirPH[i] - indirPHitrl != 1)))
      {
        // Shift read pointer and skip
        ind += nf+shift;
        continue;
      }
      
      // Reindex and /or compress
      if (indirPHi >= 0 && !removeDegeneratedFaces) // non-degenerated face
      {
        nfins = nf;
        nface[k] = nface[ind];
        for (E_Int f = shift; f < nf+shift; f++)
        {
          fidx = nface[ind+f];
          if (removeDuplicatedFaces) fidx = indirPG[fidx-1];
          nface[k+f] = fidx;
        }
      }
      else
      {
        nfins = 0;
        for (E_Int f = shift; f < nf+shift; f++)
        {
          fidx = nface[ind+f];
          if (removeDuplicatedFaces) fidx = indirPG[fidx-1];
          if (fidx == COLLAPSED) continue;
          nface[k+shift+nfins] = abs(fidx); nfins++;
        }
        if (shift == 1) nface[k] = nfins;
      }
      
      if (removeDirtyElements && array23) // array2 or array3
      {
        indPH[j] = k; j += 1;
      }
      
      // Shift read and write pointers
      ind += nf+shift; k += nfins+shift; itrl = i;
    }
    sizeEF2 = k;

    std::cout<<"sizeFN2: " << sizeFN2 << std::endl;
    std::cout<<"sizeEF2: " << sizeEF2 << std::endl;
  }

  // --- 5. Create resized connectivity ---
  if (removeOverlappingPoints ||
      removeDuplicatedFaces   || removeDuplicatedElements)
  {  
    E_Int ngonType = 1; // CGNSv3 compact array1
    if (api == 2) ngonType = 2; // CGNSv3, array2
    else if (api == 3) ngonType = 3; // force CGNSv4, array3
    E_Boolean center = false;
    tpl = K_ARRAY::buildArray3(nfld, varString, nuniquePts, nuniqueElts,
                               nuniqueFaces, "NGON", sizeFN2, sizeEF2,
                               ngonType, center, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    E_Int *ngon2 = cn2->getNGon(), *nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (array23)
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
        #pragma omp for
        for (E_Int i = 0; i < nuniquePts; i++) f2p[i] = fp[i];
      }

      // Copy compressed connectivity to cn2
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeFN2; i++) ngon2[i] = ngon[i];
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeEF2; i++) nface2[i] = nface[i];

      if (array23) // array2 or array3
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < nuniqueFaces; i++) indPG2[i] = indPG[i];
        #pragma omp for nowait
        for (E_Int i = 0; i < nuniqueElts; i++) indPH2[i] = indPH[i];
      }
    }

    RELEASESHAREDU(tpl, f2, cn2)
  }
  return tpl;
}

// Nettoyage de la connectivite et du tableau des sommets.
E_Int K_CONNECT::V_identifyOverlappingPoints(
  E_Int posx, E_Int posy, E_Int posz, 
  FldArrayF& f, E_Float tol, std::vector<E_Int>& indir
)
{
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "removeOverlappingPoints: field must have coordinates.");
    return -1;
  }
  
  size_t npts = f.getSize();
  E_Bool removeOrphans = (indir.size() == npts);
  ArrayAccessor<FldArrayF> coordAcc(f, posx, posy, posz);
  // Local indirection table for the KdTree
  std::vector<E_Int> locIndir(npts);
  #pragma omp parallel for
  for (size_t i = 0; i < npts; ++i) locIndir[i] = i;

  // Build the KdTree
  K_SEARCH::KdTree<FldArrayF> moving_tree(coordAcc, locIndir, tol, true);
  
  std::vector<triplet_t> palma;
  std::vector<std::vector<triplet_t> > palma_thrd(__NUMTHREADS__);
  std::vector<std::vector<E_Int> > onodes_thrd(__NUMTHREADS__);
  std::vector<std::vector<E_Float> > dist2_thrd(__NUMTHREADS__);

  // Loop on target nodes and get all moving nodes in sphere of radius tol.
  #pragma omp parallel
  {
    size_t osz;
    E_Int id, fi, fj;
    E_Float d2;

    id = __CURRENT_THREAD__;

    #pragma omp for schedule(static)
    for (size_t i = 0; i < npts; ++i)
    {
      onodes_thrd[id].clear(); dist2_thrd[id].clear();

      fi = i; //utarget[i];
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
  E_Int nb_merges = 0;
  if (!palma.empty()) // duplicates found
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
      if ((locIndir[fi] == fi) && (locIndir[fj] == fj))
      {
        locIndir[fj] = fi;
        nb_merges++;
      }
    }

    // Remove circular references by instructing the pointers to point to the
    // leaves, and conditionally mark orphan points
    if (removeOrphans)
    {
      for (size_t i = 0; i < npts; ++i)
      {
        // Skip orphan points
        if (indir[i] == -1) { locIndir[i] = -1; continue; }
        fi = locIndir[i];
        while (fi != locIndir[fi]) fi = locIndir[fi];
        locIndir[i] = fi;
      }
    }
    else
    {
      for (size_t i = 0; i < npts; ++i)
      {
        fi = locIndir[i];
        while (fi != locIndir[fi]) fi = locIndir[fi];
        locIndir[i] = fi;
      }
    }
  }
  else if (removeOrphans) // orphans only
  {
    #pragma omp parallel for
    for (size_t i = 0; i < npts; ++i) locIndir[i] = indir[i];
  }
  else return npts; // no orphans nor duplicates

  // Vertex set to discard orphans and/or duplicates
  indir.clear(); indir.resize(npts);
  size_t nuniquePtsRef = npts - nb_merges;
  std::unordered_set<E_Int> vertexSet;

  // Loop over all vertices
  E_Int nuniquePts = 0;
  for (size_t i = 0; i < npts; ++i)
  {
    if (locIndir[i] == -1)
    {
      // Skip orphan points
      indir[i] = -1; continue;
    }
    else if (vertexSet.insert(locIndir[i]).second)
    {
      // First time this vertex is encountered, added to the set
      indir[i] = nuniquePts; nuniquePts++;
    }
    else indir[i] = indir[locIndir[i]];
  }

  std::cout<<"nuniquePtsRef: " << nuniquePtsRef << std::endl;
  std::cout<<"vertexSet.size(): " << vertexSet.size() << std::endl;
  //assert(vertexSet.size() == nuniquePtsRef);
  return vertexSet.size();
}

// Identify dirty faces topologically
E_Int K_CONNECT::V_identifyDirtyFacesNGon(
  E_Int dim, FldArrayF& f, FldArrayI& cn, std::vector<E_Int>& indirPG,
  E_Bool removeDegeneratedFaces
)
{
  E_Int *ngon = cn.getNGon(), *indPG = cn.getIndPG();
  return K_CONNECT::V_identifyDirtyFacesNGon(dim, cn, ngon, indPG, indirPG, removeDegeneratedFaces);
}

E_Int K_CONNECT::V_identifyDirtyFacesNGon(
  E_Int dim, FldArrayI& cn, E_Int* ngon, E_Int* indPG, std::vector<E_Int>& indirPG,
  E_Bool removeDegeneratedFaces
)
{
  E_Int nv, nuniqueFaces = 1;
  E_Int nfaces = cn.getNFaces();
  indirPG.clear(); indirPG.resize(nfaces);

  // Face map to eliminate duplicates and collapsed faces
  Topology F;
  std::unordered_map<Topology, E_Int, JenkinsHash<Topology> > faceMap;
  
  if (removeDegeneratedFaces && dim > 1)
  {
    // Loop over all faces of the NGON connectivity
    for (E_Int i = 0; i < nfaces; i++)
    {
      E_Int* face = cn.getFace(i, nv, ngon, indPG);
      F.set(face, nv, removeDegeneratedFaces);
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
  E_Bool removeDegeneratedElements
)
{
  E_Int *nface = cn.getNFace(), *indPH = cn.getIndPH();
  return K_CONNECT::V_identifyDirtyElementsNGon(dim, cn, nface, indPH, indirPH, removeDegeneratedElements);
}

E_Int K_CONNECT::V_identifyDirtyElementsNGon(
  E_Int dim, FldArrayI& cn, E_Int* nface, E_Int* indPH,
  std::vector<E_Int>& indirPH,
  E_Bool removeDegeneratedElements
)
{
  E_Int nf, nuniqueElts = 1;
  E_Int nelts = cn.getNElts();
  indirPH.clear(); indirPH.resize(nelts);

  // Element map to eliminate duplicates
  Topology E;
  std::unordered_map<Topology, E_Int, JenkinsHash<Topology> > eltMap;

  if (removeDegeneratedElements && dim > 0)
  {
    // Loop over all faces of the NGON connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int* elt = cn.getElt(i, nf, nface, indPH);
      E.set(elt, nf, removeDegeneratedElements);
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
  E_Float tol, E_Bool removeOverlappingPoints,
  E_Bool removeOrphanPoints,
  E_Bool removeDuplicatedElements,
  E_Bool removeDegeneratedElements
)
{
  E_Bool removeDirtyPoints = (removeOverlappingPoints || removeOrphanPoints);
  E_Bool removeDirtyElements = (removeDuplicatedElements || removeDegeneratedElements);
  
  PyObject* tpl = NULL;
  E_Int nc = cn.getNConnect();
  E_Int nfld = f.getNfld(), npts = f.getSize(), api = f.getApi();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Get dimensionality
  E_Int dim = 3;
  if (strcmp(eltTypes[0], "NODE") == 0) dim = 0;
  else if (strcmp(eltTypes[0], "BAR") == 0) dim = 1;
  else if (strcmp(eltTypes[0], "TRI") == 0 ||
           strcmp(eltTypes[0], "QUAD") == 0) dim = 2;
  for (E_Int ic = 0; ic < nc; ic++) {std::cout<<"eltType: " << eltTypes[ic] << std::endl; delete [] eltTypes[ic];}
  std::cout<<"dim: " << dim << std::endl;

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
  if (removeOrphanPoints and dim > 0)
  {
    E_Int vidx;
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
  }

  // 1b. Identify overlapping points geometrically
  if (removeOverlappingPoints)
  {
    nuniquePts = K_CONNECT::V_identifyOverlappingPoints(posx, posy, posz, f, tol, indir);
    if (nuniquePts < 0) return NULL;
    else if (nuniquePts == npts) removeDirtyPoints = false;
  }
  std::cout<<"npts: " << npts << std::endl;
  std::cout<<"nuniquePts: " << nuniquePts << std::endl;

  // An update is necessary before topological operations in 2.
  if (removeDirtyPoints)
  {
    E_Int vidx;
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
          if (indir[vidx] != vidx)
          {
            cm(i,j) = indir[vidx]+1;
          }
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
        if (itrl == -1 || abs(indir[i]) - abs(indir[itrl]) == 1)
        {
          fp[j] = fp[i];
          j += 1; itrl = i;
        }
      }
    }

    // 1.d Resize fields
    //f.reAllocMat(nuniquePts,nfld);
  }

  // --- 2. Identify duplicated elements topologically ---
  E_Int nuniqueEltsTot = neltsTot;
  std::vector<E_Int> indirPH, nuniqueElts;
  if (removeDirtyElements)
  {
    nuniqueEltsTot = K_CONNECT::V_identifyDuplicatedElementsME(cn, indirPH, nuniqueElts, neltsTot, removeDegeneratedElements);
    if (nuniqueEltsTot == neltsTot) removeDuplicatedElements = false;
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
  std::cout<<"neltsTot: " << neltsTot << std::endl;
  std::cout<<"nuniqueEltsTot: " << nuniqueEltsTot << std::endl;

  // --- 3. Reindex & Compress connectivities ---
  E_Int k = 0; // write pointer
  E_Int itrl = -1;
  if (removeDirtyElements)
  {
    for (E_Int ic = 0; ic < nc; ic++)
    {
      // Skip connectivity if none of its elements are duplicated
      if (nuniqueElts[nc] == nepc[nc]) continue;
      
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int nvpe = cm.getNfld();
      for (E_Int i = 0; i < nelts; i++)
      {
        // Skip duplicated elements
        if (itrl != -1 && indirPH[i] - indirPH[itrl] != 1) continue;

        for (E_Int j = 1; j <= nvpe; j++) cm(k,j) = cm(i,j);
        itrl = i; k += 1;
      }
    }
  }

  // --- 4. Create resized connectivity ---
  if (removeOverlappingPoints || removeDirtyElements)
  {  
    E_Boolean center = false;
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
    RELEASESHAREDU(tpl, f2, cn2)
  }
  return tpl;
}

E_Int K_CONNECT::V_identifyDuplicatedElementsME(
  FldArrayI& cn, std::vector<E_Int>& indir, std::vector<E_Int>& nuniqueElts,
  E_Int neltsTot, E_Bool removeDegeneratedElements
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
  
  // Loop over ME connectivity
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nuniqueElts[ic] = 1;
    E_Int nvpe = cm.getNfld();
    std::vector<E_Int> elt(nvpe);
    for (E_Int i = 0; i < nelts; i++)
    {
      for (E_Int j = 1; j <= nvpe; j++) elt[j-1] = cm(i,j);
      E.set(elt, nvpe, removeDegeneratedElements);
      //if (E.isDegen_) { indir[elOffset+i] = COLLAPSED; continue; } // TODO refine - no ME support yet
      // Use insert to ensure E is initially mapped to -1 if it doesn't exist
      auto res = eltMap.insert(std::make_pair(E, -1));
      // Check the value associated with E. If it is -1, then first time this
      // element is encountered
      if (res.first->second == -1) 
      {
        res.first->second = nuniqueElts[ic]; nuniqueElts[ic]++;
      }
      indir[elOffset+i] = res.first->second;
      // Flip sign of degenerated elements
      // TODO: not all element types can afford losing one or more vertices
      // check element type and E.size_
      //if (E.isDegen_) indir[elOffset+i] *= -1;
    }
    elOffset += nelts;
    nuniqueElts[ic] -= 1;
    nuniqueEltsTot += nuniqueElts[ic];
  }
  return nuniqueEltsTot;
}

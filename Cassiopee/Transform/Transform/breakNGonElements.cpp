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

# include "transform.h"
# include "Connect/connect.h"
#include <unordered_map>

using namespace K_FLD;
using namespace std;

//=============================================================================
PyObject* K_TRANSFORM::breakElements(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString,
                               f, ni, nj, nk, cnl, eltType);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "breakElements: array is invalid.");
    return NULL;
  }

  if (strcmp(eltType, "NGON") != 0 && strcmp(eltType, "MIXED") != 0)  // BE/ME
  {
    RELEASESHAREDU(array, f, cnl);
    return array;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;

  PyObject* l;
  if (strcmp(eltType, "NGON") == 0) l = breakNGonElements(*f, *cnl, varString);
  else l = breakMixedElements(*f, *cnl, varString);

  RELEASESHAREDU(array, f, cnl);
  return l;
}

//=============================================================================
PyObject* K_TRANSFORM::breakNGonElements(FldArrayF& field, FldArrayI& cNG,
                                         char* varString)
{
  E_Int nfld = field.getNfld();
  E_Int npts = field.getSize();
  E_Int api = field.getApi();

  E_Int* ngon = cNG.getNGon(); E_Int* nface = cNG.getNFace();
  E_Int* indPG = cNG.getIndPG(); E_Int* indPH = cNG.getIndPH();
  E_Int nelts = cNG.getNElts();
  E_Int dim = cNG.getDim();
  E_Int ngonType = cNG.getNGonType();
  E_Int shift = 1; if (ngonType == 3) shift = 0;

  // Connectivity Element->Vertex
  std::vector<vector<E_Int> > cEVNGon(nelts);
  K_CONNECT::connectNG2EV(cNG, cEVNGon);

  PyObject* l = PyList_New(0);
  PyObject* tpl;

  if (dim == 1)
  {
    // In 1D, all NGon elements become BARs
    tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts, "BAR", false, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);

    #pragma omp parallel
    {
      // Fields
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = field.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      }

      // Connectivity
      FldArrayI& cm2 = *(cn2->getConnect(0));
      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        const std::vector<E_Int>& vertices = cEVNGon[i];  // sommets associes a l'elt
        cm2(i, 1) = vertices[0];
        cm2(i, 2) = vertices[1];
      }
    }

    RELEASESHAREDU(tpl, f2, cn2);
    PyList_Append(l, tpl); Py_DECREF(tpl);
    return l;
  }

  // Input NGon is either 2D or 3D and can be broken into ME + NGon
  // Vertex masks for the resulting ME and NGon connectivities, where 1 means
  // that this vertex belongs to this connectivity
  std::vector<E_Int> vMaskME(npts, 0), vMaskNG(npts, 0);

  // Number of elements per connectivity of the output ME
  // NB1: 'tmp_' is uncompressed: all possible connectivities listed
  // NB2: output NGon listed in position 0
  const E_Int nbuckets = 7;
  std::vector<E_Int> tmp_nepc2(nbuckets, 0);

  // Thread-related arrays are prefixed with 't'.
  const E_Int nthreads = __NUMTHREADS__;
  // For each thread:
  //  - tindir: maps element indices from new to old ME. Note that
  //            tindir[tid][ic].size() is the number of elements found in the
  //            BE conn. of index 'ic' of the output ME for this thread
  std::vector<std::vector<std::vector<E_Int> > > tindir(nthreads);
  //  - toffset: cumulative number of elements found in each connectivity
  std::vector<std::vector<E_Int> > toffset(nthreads);
  // - tsizeEF2: size of the output NGon Element-Face connectivity
  std::vector<E_Int> tsizeEF2(nthreads);

  // Init. thread vars
  for (E_Int tid = 0; tid < nthreads; tid++)
  {
    tindir[tid].resize(nbuckets);
    toffset[tid].resize(nbuckets);

    tsizeEF2[tid] = 0;
    for (E_Int ic = 0; ic < nbuckets; ic++)
    {
      // tindir[tid][ic].reserve(X); TODO light first pass
      toffset[tid][ic] = 0;
    }
  }

  #pragma omp parallel
  {
    E_Int nv, nf, nvpe;
    std::vector<E_Int> verticesf;  // candidats aux sommets images

    // Local thread-related arrays are prefixed with 'loc_t'
    E_Int tid = __CURRENT_THREAD__;
    std::vector<std::vector<E_Int> >& loc_tindir = tindir[tid];
    
    #pragma omp for schedule(static)
    for (E_Int i = 0; i < nelts; i++)
    {
      std::vector<E_Int>& vertices = cEVNGon[i];  // sommets associes a l'elt
      nvpe = vertices.size();
      E_Int* elt = cNG.getElt(i, nf, nface, indPH);

      if (dim == 2)
      {
        if (nf == 3)  // TRI
        {
          // Tag vertices
          for (E_Int v = 0; v < nvpe; v++) vMaskME[vertices[v]-1] = 1;
          loc_tindir[1].push_back(i);
        }
        else if (nf == 4)  // QUAD
        {
          E_Int nv, vert0, vert1, vert2, vert3, vert4;
          E_Int fidx = elt[0];
          E_Int* face = cNG.getFace(fidx-1, nv, ngon, indPG);
          vert1 = face[0]; vert2 = face[1];

          // verification de la coherence de la numerotation des indices
          verticesf.clear();
          for (E_Int v = 0; v < nvpe; v++)
          {
            vert0 = vertices[v];
            if (vert0 != vert1 && vert0 != vert2) verticesf.push_back(vert0);
            vMaskME[vert0-1] = 1;  // Tag vertex
          }

          vert3 = K_CONNECT::image(vert2, fidx, i, verticesf, cNG,
                                   ngon, nface, indPG, indPH);
          if (vert3 == -1) goto ngonLabel;
          vert4 = K_CONNECT::image(vert1, fidx, i, verticesf, cNG,
                                   ngon, nface, indPG, indPH);
          if (vert4 == -1) goto ngonLabel;

          // Vertex order may have changed - update cEVNGon
          vertices[0] = vert1; vertices[1] = vert2; vertices[2] = vert3;
          vertices[3] = vert4;
          loc_tindir[2].push_back(i);
        }
        else goto ngonLabel;  // 2D polygon with more than 4 faces
      }
      else if (nf == 4)  // TETRA
      {
        // Recherche de la premiere face tri de l elt
        E_Int nv, vert0, vert1, vert2, vert3, vert4 = -1;
        E_Int* face = cNG.getFace(elt[0]-1, nv, ngon, indPG);
        vert1 = face[0]; vert2 = face[1]; vert3 = face[2];

        for (E_Int v = 0; v < nvpe; v++)
        {
          vert0 = vertices[v];
          if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3)
          {
            vert4 = vert0; break;
          }
        }

        // Vertex order may have changed - update cEVNGon
        for (E_Int v = 0; v < nvpe; v++) vMaskME[vertices[v]-1] = 1;
        vertices[0] = vert1; vertices[1] = vert2; vertices[2] = vert3;
        vertices[3] = vert4;
        loc_tindir[3].push_back(i);
      }
      else if (nf == 5)  // PENTA, PYRA or NGon
      {
        E_Int nv, vert0;
        E_Int vert1 = -1, vert2 = -1, vert3 = -1, vert4 = -1, vert5 = -1;
        E_Int nbnodes = 0;  // number of non-unique vertices composing the elt
        for (E_Int j = 0; j < nf; j++)
        {
          cNG.getFace(elt[j]-1, nv, ngon, indPG);
          nbnodes += nv;
        }

        if (nbnodes == 16)  // PYRA
        {
          vert5 = -1;
          for (E_Int j = 0; j < nf; j++)
          {
            E_Int* face = cNG.getFace(elt[j]-1, nv, ngon, indPG);
            if (nv == 4)  // quad base face found
            {
              vert1 = face[0]; vert2 = face[1]; vert3 = face[2]; vert4 = face[3];
              // verification de la coherence de la numerotation des indices
              for (E_Int v = 0; v < nvpe; v++)
              {
                vert0 = vertices[v];
                if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3 && vert0 != vert4)
                {
                  vert5 = vert0; break;
                }
              }
            }
            if (vert5 != -1) break;  // a PYRA elt was found
          }

          if (vert5 == -1) goto ngonLabel;

          // Vertex order may have changed - update cEVNGon
          for (E_Int v = 0; v < nvpe; v++) vMaskME[vertices[v]-1] = 1;
          vertices[0] = vert1; vertices[1] = vert2; vertices[2] = vert3;
          vertices[3] = vert4; vertices[4] = vert5;
          loc_tindir[4].push_back(i);
        }
        else if (nbnodes == 18) // PENTA
        {
          E_Int fidx, vert6 = -1;
          for (E_Int j = 0; j < nf; j++)
          {
            fidx = elt[j];
            E_Int* face = cNG.getFace(fidx-1, nv, ngon, indPG);
            if (nv == 3)  // tri base face found
            {
              vert1 = face[0]; vert2 = face[1]; vert3 = face[2];

              // verification de la coherence de la numerotation des indices
              verticesf.clear();
              for (E_Int v = 0; v < nvpe; v++)
              {
                vert0 = vertices[v];
                if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3)
                  verticesf.push_back(vert0);
              }

              vert4 = K_CONNECT::image(vert1, fidx, i, verticesf, cNG,
                                       ngon, nface, indPG, indPH);
              if (vert4 == -1) goto ngonLabel;
              vert5 = K_CONNECT::image(vert2, fidx, i, verticesf, cNG,
                                       ngon, nface, indPG, indPH);
              if (vert5 == -1) goto ngonLabel;
              vert6 = K_CONNECT::image(vert3, fidx, i, verticesf, cNG,
                                       ngon, nface, indPG, indPH);
              if (vert6 == -1) goto ngonLabel;
              break;
            }
          }

          // Vertex order may have changed - update cEVNGon
          for (E_Int v = 0; v < nvpe; v++) vMaskME[vertices[v]-1] = 1;
          vertices[0] = vert1; vertices[1] = vert2; vertices[2] = vert3;
          vertices[3] = vert4; vertices[4] = vert5; vertices[5] = vert6;
          loc_tindir[5].push_back(i);
        }
        else goto ngonLabel;  // 5-faced polygon that is not a PYRA nor a PENTA
      }
      else if (nf == 6) // HEXA
      {
        E_Int nv, vert0, vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8;
        E_Int fidx = elt[0];
        E_Int* face = cNG.getFace(fidx-1, nv, ngon, indPG);
        vert1 = face[0]; vert2 = face[1]; vert3 = face[2]; vert4 = face[3];

        // verification de la coherence de la numerotation des indices
        verticesf.clear();
        for (E_Int v = 0; v < nvpe; v++)
        {
          vert0 = vertices[v];
          if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3 && vert0 != vert4)
            verticesf.push_back(vert0);
        }

        vert5 = K_CONNECT::image(vert1, fidx, i, verticesf, cNG,
                                 ngon, nface, indPG, indPH);
        if (vert5 == -1) goto ngonLabel;
        vert6 = K_CONNECT::image(vert2, fidx, i, verticesf, cNG,
                                 ngon, nface, indPG, indPH);
        if (vert6 == -1) goto ngonLabel;
        vert7 = K_CONNECT::image(vert3, fidx, i, verticesf, cNG,
                                 ngon, nface, indPG, indPH);
        if (vert7 == -1) goto ngonLabel;
        vert8 = K_CONNECT::image(vert4, fidx, i, verticesf, cNG,
                                 ngon, nface, indPG, indPH);
        if (vert8 == -1) goto ngonLabel;
        
        // Vertex order may have changed - update cEVNGon
        for (E_Int v = 0; v < nvpe; v++) vMaskME[vertices[v]-1] = 1;
        vertices[0] = vert1; vertices[1] = vert2; vertices[2] = vert3;
        vertices[3] = vert4; vertices[4] = vert5; vertices[5] = vert6;
        vertices[6] = vert7; vertices[7] = vert8;
        loc_tindir[6].push_back(i);
      }
      else  // element remains an NGon
      {
        ngonLabel:;
        tsizeEF2[tid] += nf + shift;
        for (E_Int j = 0; j < nf; j++)
        {
          E_Int* face = cNG.getFace(elt[j]-1, nv, ngon, indPG);
          for (E_Int v = 0; v < nv; v++) vMaskNG[face[v]-1] = 1;
        }
        loc_tindir[0].push_back(i);
      }
    }
  }

  // Transform the vertex masks comprised of zeros and ones into vertex maps
  // from old to new connectivities, and get the number of vertices for each
  // output connectivity
  E_Int nptsME = K_CONNECT::prefixSum(vMaskME);
  E_Int nptsNG = K_CONNECT::prefixSum(vMaskNG);

  // Sum over all threads to fill tmp_nepc2 and update sizeEF2
  E_Int sizeEF2 = 0;
  for (E_Int tid = 0; tid < nthreads; tid++)
  {
    for (E_Int ic = 0; ic < nbuckets; ic++)
      tmp_nepc2[ic] += (E_Int)tindir[tid][ic].size();
    sizeEF2 += tsizeEF2[tid];
  }

  // Compute thread element offsets in the output ME for each connectivity
  // toffset is a cumulative thread-tmp_nepc2 over all conns
  // Used to build cm2 using multiple threads
  for (E_Int tid = 1; tid < nthreads; tid++)
    for (E_Int ic = 0; ic < nbuckets; ic++)
      toffset[tid][ic] = toffset[tid-1][ic] + (E_Int)tindir[tid-1][ic].size();

  // Build output ME
  if (nptsME > 0)
  {
    // Build new eltType from connectivities that have at least one element
    E_Int nc2 = 0;
    char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
    eltType2[0] = '\0';
    if (tmp_nepc2[1] > 0)
    {
      strcat(eltType2, "TRI");
      nc2++;
    }
    if (tmp_nepc2[2] > 0)
    {
      if (nc2 > 0) strcat(eltType2, ",");
      strcat(eltType2, "QUAD");
      nc2++;
    }
    if (tmp_nepc2[3] > 0)
    {
      if (nc2 > 0) strcat(eltType2, ",");
      strcat(eltType2, "TETRA");
      nc2++;
    }
    if (tmp_nepc2[4] > 0)
    {
      if (nc2 > 0) strcat(eltType2, ",");
      strcat(eltType2, "PYRA");
      nc2++;
    }
    if (tmp_nepc2[5] > 0)
    {
      if (nc2 > 0) strcat(eltType2, ",");
      strcat(eltType2, "PENTA");
      nc2++;
    }
    if (tmp_nepc2[6] > 0)
    {
      if (nc2 > 0) strcat(eltType2, ",");
      strcat(eltType2, "HEXA");
      nc2++;
    }

    std::cout << "api = " << api << std::endl;
    std::cout << "eltType2 = " << eltType2 << std::endl;
    std::cout << "nc2 = " << nc2 << std::endl;
    std::cout << "nptsME = " << nptsME << std::endl;
    for (E_Int ic = 0; ic < nbuckets; ic++)
      std::cout << "tmp_nepc2["<<ic<<"] = " << tmp_nepc2[ic] << std::endl;

    // Compress the number of elements per connectivity of the output ME, ie,
    // drop connectivities containing no elements
    std::vector<E_Int> nepc2(nc2);
    nc2 = 0;
    for (E_Int ic = 1; ic < nbuckets; ic++)  // from TRI (1) to HEXA (6)
    {
      if (tmp_nepc2[ic] > 0) { nepc2[nc2] = tmp_nepc2[ic]; nc2++; }
    }

    // Build new ME connectivity
    // if (nc2 > 1) api = 3;
    tpl = K_ARRAY::buildArray3(nfld, varString, nptsME,
                               nepc2, eltType2, false, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);

    #pragma omp parallel
    {
      E_Int nvpe, loc_nelts, ind, offW, eidx;
      E_Int tid = __CURRENT_THREAD__;
      std::vector<std::vector<E_Int> >& loc_tindir = tindir[tid];
      std::vector<E_Int>& loc_toffset = toffset[tid];

      // Copy fields
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = field.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < npts; i++)
        {
          ind = vMaskME[i];
          if (ind > 0) f2p[ind-1] = fp[i];
        }
      }

      // Copy connectivities
      E_Int ic2 = 0;
      for (E_Int ic = 1; ic < nbuckets; ic++)  // from TRI (1) to HEXA (6)
      {
        if (tmp_nepc2[ic] == 0) continue;  // no elements in this conn., skip
        FldArrayI& cm2 = *(cn2->getConnect(ic2));
        nvpe = cm2.getNfld();
        std::vector<E_Int>& loc_tindirIc = loc_tindir[ic];
        loc_nelts = loc_tindirIc.size();
        offW = loc_toffset[ic];

        for (E_Int i = 0; i < loc_nelts; i++)
        {
          eidx = loc_tindirIc[i];  // global element index
          const std::vector<E_Int>& vertices = cEVNGon[eidx];
          for (E_Int j = 0; j < nvpe; j++)
          {
            ind = vertices[j] - 1;
            cm2(offW+i, j+1) = vMaskME[ind];
          }
        }
        ic2++;
      }
    }

    RELEASESHAREDU(tpl, f2, cn2);
    delete [] eltType2;
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  vMaskME.clear(); vMaskME.shrink_to_fit();

  // Build output NGon
  if (nptsNG > 0)
  {
    // Hash NGon faces to get their unique count
    E_Int nv, nf, fidx, eidx, nfaces2 = 0, sizeFN2 = 0;
    std::unordered_map<E_Int, E_Int> fMapNG;
  
    // Loop over all input NGon elements which remain NGon
    for (E_Int tid = 0; tid < nthreads; tid++)
    {
      std::vector<E_Int>& loc_tindirNG = tindir[tid][0];
      for (size_t e = 0; e < loc_tindirNG.size(); e++)
      {
        eidx = loc_tindirNG[e];  // global element index
        E_Int* elt = cNG.getElt(eidx, nf, nface, indPH);
        for (E_Int j = 0; j < nf; j++)
        {
          fidx = elt[j] - 1;
          auto resF = fMapNG.insert({fidx, nfaces2});
          if (resF.first->second == nfaces2)  // first time this face is encountered
          {
            cNG.getFace(fidx, nv, ngon, indPG);
            sizeFN2 += nv + shift;
            nfaces2++;
          }
        }
      }
    }

    // Build new NGon connectivity
    E_Int nelts2 = tmp_nepc2[0];
    tpl = K_ARRAY::buildArray3(nfld, varString, nptsNG, nelts2,
                               nfaces2, "NGON", sizeFN2, sizeEF2,
                               ngonType, false, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    E_Int *ngon2 = cn2->getNGon(), *nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (ngonType == 2 || ngonType == 3)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }

    #pragma omp parallel
    {
      E_Int ind;
      // TODO: multithread copy of connectivity if possible
      // E_Int eidx;
      // E_Int tid = __CURRENT_THREAD__;
      // std::vector<E_Int>& loc_tindirNG = tindir[tid][0];
      // E_Int loc_nelts = loc_tindirNG.size();
      // E_Int offW = toffset[tid][0];

      // Copy fields
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = field.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < npts; i++)
        {
          ind = vMaskNG[i];
          if (ind > 0) f2p[ind-1] = fp[i];
        }
      }
    }

    // Copy connectivity
    E_Int ind, eidx2 = 0, ind3 = 0, ind4 = 0, ind5 = 1;

    if (ngonType == 2 || ngonType == 3) { indPG2[0] = 0; indPH2[0] = 0; }
    for (E_Int tid = 0; tid < nthreads; tid++)
    {
      std::vector<E_Int>& loc_tindirNG = tindir[tid][0];
      for (size_t e = 0; e < loc_tindirNG.size(); e++)
      {
        eidx = loc_tindirNG[e];  // global element index of the input NGon
        E_Int* elt = cNG.getElt(eidx, nf, nface, indPH);
        nface2[ind4] = nf;
        if ((ngonType == 2 || ngonType == 3) && (eidx2+1 < nelts2)) indPH2[eidx2+1] = nf;

        for (E_Int j = 0; j < nf; j++)
        {
          fidx = elt[j] - 1;
          E_Int* face = cNG.getFace(fidx, nv, ngon, indPG);
          ngon2[ind3] = nv;
          if ((ngonType == 2 || ngonType == 3) && (eidx2+1 < nfaces2))
          {
            indPG2[ind5] = nv; ind5++;
          }
          for (E_Int v = 0; v < nv; v++)
          {
            ind = face[v] - 1;
            ngon2[ind3+v+shift] = vMaskNG[ind];
          }
          ind3 += nv + shift;
          nface2[ind4+j+shift] = fMapNG[fidx] + 1;
        }
        ind4 += nf + shift;
        eidx2++;  // global element index of the output NGon
      }
    }

    RELEASESHAREDU(tpl, f2, cn2);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  
  return l;
}

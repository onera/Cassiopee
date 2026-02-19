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

// selectExteriorVertices
#include <unordered_map>
#include "post.h"

//=============================================================================
/* Find exterior vertices and return an unstructured NODE connectivity (BE) */
// ============================================================================
PyObject* K_POST::selectExteriorVertices(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* indices;
  if (!PYPARSETUPLE_(args, OO_, &array, &indices)) return NULL;
  
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);

  PyObject* tpl = NULL;
  if (res == 1)
  {
    tpl = exteriorVerticesStructured(varString, *f, ni, nj, nk, indices);
    RELEASESHAREDS(array, f);
  }
  else if (res == 2)
  {
    if (K_STRING::cmp(eltType, "NODE") == 0)
    {
      RELEASESHAREDU(array, f, cn);
      return array;
    }
    else if (strcmp(eltType, "NGON") == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                    "exteriorVertices: NGON array not supported yet.");
      // E_Int dim = cn->getDim();
      // if (dim == 2) tpl = selectExteriorVerticesNGon2D(varString, *f, *cn, indices);
      // else tpl = selectExteriorVerticesNGon3D(varString, *f, *cn, indices);
    }
    else tpl = selectExteriorVerticesME(varString, *f, *cn, eltType, indices);
    RELEASESHAREDU(array, f, cn);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorVertices: invalid array.");
    RELEASESHAREDB(res, array, f, cn); 
  } 
  return tpl;
}
//=============================================================================
PyObject* K_POST::exteriorVerticesStructured(char* varString, FldArrayF& f, 
                                             E_Int ni, E_Int nj, E_Int nk, 
                                             PyObject* indices)
{
  E_Int api = f.getApi();
  E_Int nfld = f.getNfld();
  E_Int npts = f.getSize();

  PyObject* tpl = NULL;
  FldArrayF* f2;
  E_Bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;
  PyObject* indir = NULL;
  E_Int* indirp = NULL;

  // 1D arrays
  if ((ni == 1 && nj == 1) || (nj == 1 && nk == 1) || (ni == 1 && nk == 1))
  {
    E_Int nptsExt = 2;
    tpl = K_ARRAY::buildArray3(nfld, varString, nptsExt, 0, "NODE", false, api);
    K_ARRAY::getFromArray3(tpl, f2);

    for (E_Int n = 1; n <= nfld; n++)
    {
      (*f2)(0,n) = f(0,n);
      (*f2)(1,n) = f(npts-1,n);
    }

    if (boolIndir)
    {
      indir = K_NUMPY::buildNumpyArray(nptsExt, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);
      indirp[0] = 1; indirp[1] = npts;
    }
  }
  else if (ni == 1 || nj == 1 || nk == 1) // arrays 2D
  {
    E_Int p = 0, q = 0;
    if (ni == 1) { p = nj; q = nk; }
    if (nj == 1) { p = ni; q = nk; }
    if (nk == 1) { p = ni; q = nj; }
    E_Int ninti = 2*(q - 1); E_Int nintj = 2*(p - 1);
    E_Int nfacesExt = ninti + nintj;
    E_Int nptsExt = 2*p + 2*q - 4;

    tpl = K_ARRAY::buildArray3(nfld, varString, nptsExt, 0, "NODE", false, api);
    K_ARRAY::getFromArray3(tpl, f2);

    if (boolIndir)
    {
      indir = K_NUMPY::buildNumpyArray(nptsExt, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);
    }

    #pragma omp parallel
    {
      E_Int ind;
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        // Border j=1, i=1..p-1
        #pragma omp for nowait
        for (E_Int i = 0; i < p-1; i++) f2p[i] = fp[i];
        // Border i=p, j=0..q-1
        #pragma omp for nowait
        for (E_Int j = 0; j < q-1; j++)
        { 
          ind = p-1+j;
          f2p[ind] = fp[p-1+j*p];
        }
        // Border j=jmax, i=p,0
        #pragma omp for nowait
        for (E_Int i = p-1; i > 0; i--)
        { 
          ind = nintj+q-1-i;
          f2p[ind] = fp[i+(q-1)*p];
        }
        // Border i=1, j=q,1
        #pragma omp for nowait
        for (E_Int j = q-1; j > 0; j--)
        { 
          ind = nfacesExt-j;
          f2p[ind] = fp[j*ni];
        }
      }
    }
  }
  else
  {
    E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
    E_Int nptsExt = 2*(ni*nj + ni*nk + nj*nk) - 4*(ni + nj + nk) + 8;

    tpl = K_ARRAY::buildArray3(nfld, varString, nptsExt, 0, "NODE", false, api);
    K_ARRAY::getFromArray3(tpl, f2);

    if (boolIndir)
    {
      indir  = K_NUMPY::buildNumpyArray(nptsExt, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);       
    }

    #pragma omp parallel
    {
      E_Int nj2 = nj-2; E_Int nk2 = nk-2;
      E_Int ind, ind2, ns;
      
      // face k = 0
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < ni*nj; i++) f2p[i] = fp[i];
      }

      // face k = nk-1
      ns = ni*nj;
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for collapse(2) nowait
        for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        { 
          ind = ns + j*ni + i;
          ind2 = i + j*ni + nk1*ni*nj;
          f2p[ind] = fp[ind2];
        }
      }

      // face j = 0
      ns += ni*nj;
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for collapse(2) nowait
        for (E_Int k = 1; k < nk1; k++)
        for (E_Int i = 0; i < ni; i++)
        { 
          ind = ns + (k-1)*ni + i;
          ind2 = i + k*ni*nj;
          f2p[ind] = fp[ind2];
        }
      }
      
      // face j = nj-1
      ns += ni*(nk-2);
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for collapse(2) nowait
        for (E_Int k = 1; k < nk1; k++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = ns + (k-1)*ni + i;
          ind2 = i + nj1*ni + k*ni*nj;
          f2p[ind] = fp[ind2];
        }
      }

      // face i = 0
      ns += ni*nk2;
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for collapse(2) nowait
        for (E_Int k = 1; k < nk1; k++)
        for (E_Int j = 1; j < nj1; j++)
        { 
          ind = ns + (k-1)*nj2 + (j-1);
          ind2 = j*ni + k*ni*nj;
          f2p[ind] = fp[ind2];
        }
      }

      // face i = ni1
      ns += (nj-2)*(nk-2);
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for collapse(2) nowait
        for (E_Int k = 1; k < nk1; k++)
        for (E_Int j = 1; j < nj1; j++)
        { 
          ind = ns + (k-1)*nj2 + (j-1);
          ind2 = ni1 + j*ni + k*ni*nj;
          f2p[ind] = fp[ind2];
        }
      }
    }
  }

  if (boolIndir)
  {
    PyList_Append(indices, indir); Py_DECREF(indir);
  }

  RELEASESHAREDS(tpl, f2);
  return tpl;
}

//=============================================================================
// IN: varString: varString de f
// IN: f: les champs
// IN: cn: connectivite NGon
// IN: boolIndir: si true, sort un tableau d'indirection des vertices en plus
// OUT: array des vertices exterieurs
// CAS NGON 3D
//==============================================================================
PyObject* K_POST::selectExteriorVerticesNGon3D(char* varString, FldArrayF& f, 
                                               FldArrayI& cn, PyObject* indices)
{
  // E_Bool boolIndir = false;
  // if (indices != Py_None) boolIndir = true;
  // PyObject* indir = NULL;
  // E_Int* indirp = NULL;

  // E_Int npts = f.getSize(), nfld = f.getNfld(), api = f.getApi();
  // E_Int* ngon = cn.getNGon(); E_Int* indPG = cn.getIndPG();
  // E_Int nfaces = cn.getNFaces();

  // // cFE: connectivite face/elements.
  // // Si une face n'a pas d'element gauche ou droit, retourne 0 pour 
  // // cet element.
  // FldArrayI cFE; K_CONNECT::connectNG2FE(cn, cFE);
  // E_Int* cFE1 = cFE.begin(1);
  // E_Int* cFE2 = cFE.begin(2);
  // E_Int fidx, e1, e2, nbnodes, nptsExt;
  // E_Int indvertn = 1, indedgen = 1;

  // // Calcul du nombre de points uniques, aretes uniques, et faces dans la
  // // nouvelle connectivite 2D
  // E_Int nedgesExt = 0, nfacesExt = 0;
  // vector<E_Int> indirVertices(npts, -1);

  // // Les vertices sont mappes pour calculer leur table d'indirection sachant 
  // // que leur nombre n'est pas connu a priori et qu'ils ne sont pas parcourus
  // // dans l'ordre.
  // std::unordered_map<E_Int, E_Int> vertexMap; vertexMap.reserve(npts);
  // // Les aretes sont hashees pour determiner le nombre unique d'aretes et
  // // construire la nouvelle connectivite 2D
  // vector<E_Int> edge(2);
  // std::pair<E_Int, E_Bool> initEdge(-1, false);
  // TopologyOpt E;
  // std::unordered_map<TopologyOpt, std::pair<E_Int, E_Bool>, BernsteinHash<TopologyOpt> > edgeMap;

  // for (E_Int i = 0; i < nfaces; i++)
  // {
  //   e1 = cFE1[i]; // element voisin 1
  //   e2 = cFE2[i]; // element voisin 2
  //   if ((e1 == 0 && e2 != 0) || (e2 == 0 && e1 != 0)) // exterior element
  //   {
  //     E_Int* face = cn.getFace(i, nbnodes, ngon, indPG);
      
  //     for (E_Int p = 0; p < nbnodes; p++)
  //     {
  //       edge[0] = face[p]-1;
  //       edge[1] = face[(p+1)%nbnodes]-1;
  //       E.set(edge.data(), 2);
  //       // Ensure V and E are initially mapped to an initial value if either
  //       // doesn't exist
  //       auto resV = vertexMap.insert(std::make_pair(edge[0], 0));
  //       auto resE = edgeMap.insert(std::make_pair(E, initEdge));
  //       // Increment the value associated with V. If it is 1, then first
  //       // time this vertex is met, set indirVertices
  //       if (++resV.first->second == 1)
  //       {
  //         indirVertices[edge[0]] = indvertn;
  //         indvertn++;
  //       }
  //       // If the value associated with E is -1, then first time this edge
  //       // is encountered, set to current unique edge count
  //       if (resE.first->second.first == -1)
  //       {
  //         resE.first->second.first = indedgen;
  //         indedgen++;
  //       }
  //     }
  //   }
  // }

  // nptsExt = vertexMap.size();
  // cFE.malloc(0);

  // // Build new NODE connectivity
  // PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, nptsExt,
  //                                      0, "NODE", false, api);
  // FldArrayF* f2;
  // K_ARRAY::getFromArray3(tpl, f2);

  // #pragma omp parallel
  // {
  //   E_Int ind;
  //   // Copy fields
  //   for(E_Int n = 1; n <= nfld; n++)
  //   {
  //     E_Float* fp = f.begin(n);
  //     E_Float* f2p = f2->begin(n);
  //     #pragma omp for nowait
  //     for (E_Int i = 0; i < npts; i++)
  //     {
  //       ind = indirVertices[i]-1;
  //       if (ind > -1) f2p[ind] = fp[i];
  //     }
  //   }
  // }
  
  // if (boolIndir)
  // {
  //   indir = K_NUMPY::buildNumpyArray(nptsExt, 1, 1, 0);
  //   indirp = K_NUMPY::getNumpyPtrI(indir);

  //   for (E_Int i = 0; i < npts; i++)
  //   {
  //     E_Int indv = indirVertices[i];
  //     if (indv > 0) indirp[indv-1] = indv;
  //   }
    
  //   PyList_Append(indices, indir); Py_DECREF(indir);
  // }
  // RELEASESHAREDS(tpl, f2);
  // return tpl;
  return NULL;
}

//=============================================================================
// IN: varString: varString de f
// IN: f: le champ
// IN: cn: connectivite NGon
// IN: boolIndir: si true, sort un tableau d'indirection des vertices en plus
// OUT: array des vertices exterieurs
// CAS NGON 2D
//==============================================================================
PyObject* K_POST::selectExteriorVerticesNGon2D(char* varString, FldArrayF& f, 
                                               FldArrayI& cn, PyObject* indices)
{
  // // CAS 2D
  // E_Bool boolIndir = false;
  // if (indices != Py_None) boolIndir = true;

  // // Acces non universel sur le ptrs
  // E_Int* ngon = cn.getNGon();
  // E_Int* indPG = cn.getIndPG();
  // E_Int nfaces = cn.getNFaces();
  // E_Int npts = f.getSize(), nfld = f.getNfld(), api = f.getApi();
  // E_Int ngonType = cn.getNGonType();
  // E_Int shift = 1; if (ngonType == 3) shift = 0;

  // // cFE: connectivite face/elements.
  // // Si une face n'a pas d'element gauche ou droit, retourne 0 pour 
  // // cet element. 
  // FldArrayI cFE; K_CONNECT::connectNG2FE(cn, cFE);
  // E_Int* cFE1 = cFE.begin(1);
  // E_Int* cFE2 = cFE.begin(2);
  // E_Int e1, e2, dummy, nptsExt;
  // E_Int indvertn = 1;

  // // Calcul du nombre de points uniques et aretes uniques dans la
  // // nouvelle connectivite 1D
  // vector<E_Int> indirVertices(npts, -1);

  // // Les vertices sont mappes pour calculer leur table d'indirection sachant 
  // // que leur nombre n'est pas connu a priori et qu'ils ne sont pas parcourus
  // // dans l'ordre.
  // std::unordered_map<E_Int, E_Int > vertexMap; vertexMap.reserve(npts);
  // // Les aretes sont hashees pour determiner le nombre unique d'aretes et
  // // ainsi construire les "elements" de la nouvelle connectivite 1D
  // E_Int nedgesExt = 0;
  // vector<E_Int> edge(2);
  // vector<E_Int> exteriorEdges;
  // TopologyOpt E;
  // std::unordered_map<TopologyOpt, E_Int, BernsteinHash<TopologyOpt> > edgeMap;

  // for (E_Int i = 0; i < nfaces; i++)
  // {
  //   e1 = cFE1[i]; // element voisin 1
  //   e2 = cFE2[i]; // element voisin 2
  //   E_Int* face = cn.getFace(i, dummy, ngon, indPG);
  //   if ((e1 == 0 && e2 != 0) || (e2 == 0 && e1 != 0))
  //   {
  //     // Increment the value associated with V/E. If it is 1, then first
  //     // time this vertex/edge is encountered
  //     edge[0] = face[0]-1;
  //     auto resV = vertexMap.insert(std::make_pair(edge[0], 0));
  //     if (++resV.first->second == 1)
  //     {
  //       indirVertices[edge[0]] = indvertn;
  //       indvertn++;
  //     }

  //     edge[1] = face[1]-1;
  //     resV = vertexMap.insert(std::make_pair(edge[1], 0));
  //     if (++resV.first->second == 1)
  //     {
  //       indirVertices[edge[1]] = indvertn;
  //       indvertn++;
  //     }

  //     E.set(edge.data(), 2);
  //     //E.set(edge);
      
  //     auto resE = edgeMap.insert(std::make_pair(E, 0));
  //     if (++resE.first->second == 1)
  //     {
  //       exteriorEdges.push_back(i+1);
  //     }
  //   }
  // }

  // nptsExt = vertexMap.size();
  // sizeFN2 = (1+shift)*nptsExt;
  // nedgesExt = edgeMap.size();
  // sizeEF2 = (2+shift)*nedgesExt;

  // PyObject* indir = NULL;
  // E_Int* indirp = NULL;
  // if (boolIndir)
  // {
  //   indir = K_NUMPY::buildNumpyArray(nedgesExt, 1, 1, 0);
  //   indirp = K_NUMPY::getNumpyPtrI(indir);
  // }
  // cFE.malloc(0);
  // // Calcul des nouvelles connectivites Elmt/Faces et Face/Noeuds
  // E_Bool center = false;
  // PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, nptsExt, nedgesExt,
  //                                      nptsExt, "NGON", sizeFN2, sizeEF2,
  //                                      ngonType, center, api);
  // FldArrayF* f2; FldArrayI* cn2;
  // K_ARRAY::getFromArray3(tpl, f2, cn2);
  // E_Int* ngon2 = cn2->getNGon();
  // E_Int* nface2 = cn2->getNFace();
  // E_Int *indPG2 = NULL, *indPH2 = NULL;
  // if (ngonType == 2 || ngonType == 3) // set offsets
  // {
  //   indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
  // }
  
  // #pragma omp parallel
  // {
  //   E_Int ind, fidx, v1, v2;

  //   #pragma omp for nowait
  //   for(E_Int i = 0; i < nptsExt; i++)
  //   {
  //     ind = i*(1+shift);
  //     ngon2[ind] = 1;
  //     ngon2[ind+shift] = i+1;
  //   }
  //   #pragma omp for nowait
  //   for(E_Int i = 0; i < nedgesExt; i++)
  //   {
  //     fidx = exteriorEdges[i];
  //     E_Int* face = cn.getFace(fidx-1, dummy, ngon, indPG);
  //     if (boolIndir) indirp[i] = fidx;
  //     v1 = face[0]-1; v2 = face[1]-1;

  //     ind = i*(2+shift);
  //     nface2[ind] = 2;
  //     nface2[ind+shift] = indirVertices[v1];
  //     nface2[ind+1+shift] = indirVertices[v2];
  //   }

  //   if (ngonType == 2 || ngonType == 3)
  //   {
  //     #pragma omp for nowait
  //     for(E_Int i = 0; i < nptsExt; i++) indPG2[i] = 1;
  //     #pragma omp for nowait
  //     for(E_Int i = 0; i < nedgesExt; i++) indPH2[i] = 2;
  //   }
  
  //   for(E_Int n = 1; n <= nfld; n++)
  //   {
  //     E_Float* fp = f.begin(n);
  //     E_Float* f2p = f2->begin(n);
  //     #pragma omp for nowait
  //     for (E_Int i = 0; i < npts; i++)
  //     {
  //       ind = indirVertices[i]-1;
  //       if (ind > -1) f2p[ind] = fp[i];
  //     }
  //   }
  // }

  // RELEASESHAREDU(tpl, f2, cn2);
  // if (boolIndir) 
  // {
  //   PyList_Append(indices, indir);  Py_DECREF(indir);
  // }
  // return tpl;
  return NULL;
}

//=============================================================================
// Recherche topologique des vertices exterieurs utilisant le hashing des faces
//=============================================================================
PyObject* K_POST::selectExteriorVerticesME(char* varString, FldArrayF& f,
                                           FldArrayI& cn, char* eltType,
                                           PyObject* indices)
{
  // Numpy of indices of exterior vertices
  E_Bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;
  PyObject* indicesVertices = NULL;
  E_Int* indicesv = NULL;

  E_Int nc = cn.getNConnect();
  E_Int nfld = f.getNfld();
  E_Int api = f.getApi();
  E_Int npts = f.getSize();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Compute number of faces per element, nfpe
  std::vector<E_Int> nfpe;
  std::vector<E_Int> cumnfpc(nc+1); cumnfpc[0] = 0;  // cumulative number of faces per conn.
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, true);
  if (ierr != 0) return NULL;
  
  // Compute total number of faces across all connectivities
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    cumnfpc[ic+1] = cumnfpc[ic] + nelts*nfpe[ic];
  }
  E_Int ntotFaces = cumnfpc[nc];

  // Hash faces: if a face is only visited once it is an external face
  std::vector<E_Int> faceMask(ntotFaces);  // 0: interior, 1: exterior
  TopologyOpt F;
  std::unordered_map<TopologyOpt, E_Int, BernsteinHash<TopologyOpt> > faceMap;
  faceMap.reserve(ntotFaces);
  E_Int face[5];

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int nelts, nvpf, fidx;
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    nelts = cm.getSize();
    std::vector<std::vector<E_Int> > facets;
    K_CONNECT::getEVFacets(facets, eltTypes[ic], false);
  
    for (E_Int i = 0; i < nelts; i++)
    {
      // Loop over each facet of this element
      for (E_Int f = 0; f < nfpe[ic]; f++)
      {
        fidx = cumnfpc[ic] + i*nfpe[ic] + f;  // global face index
        nvpf = facets[f].size();  // number of vertices per facet
        // Fill face and insert in map
        for (E_Int j = 0; j < nvpf; j++) face[j] = cm(i, facets[f][j]);
        F.set(face, nvpf);
        auto res = faceMap.insert(std::make_pair(F, fidx));
        if (res.second) faceMask[fidx] = 1;  // first occurence of that face: tag as exterior
        else
        {
          // duplicate: this face and the one found in the map are interior faces
          faceMask[fidx] = 0;
          faceMask[res.first->second] = 0;
        }
      }
    }
  }

  // Free memory
  faceMap.clear(); faceMap.rehash(0);

  // In a first pass, tag vertex indices that belong to exterior faces
  std::vector<E_Int> vindir(npts, 0);

  #pragma omp parallel
  {
    E_Int nelts, nvpf, fidx, indv;
    std::vector<std::vector<E_Int> > facets;

    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      K_CONNECT::getEVFacets(facets, eltTypes[ic], false);

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        // Loop over each facet of this element
        for (E_Int f = 0; f < nfpe[ic]; f++)
        {
          fidx = cumnfpc[ic] + i*nfpe[ic] + f;  // global face index
          if (faceMask[fidx] == 1)  // exterior face found
          {
            nvpf = facets[f].size();  // number of vertices per face
            // Tag vertices of that faces as exterior vertices
            for (E_Int j = 0; j < nvpf; j++)
            {
              indv = cm(i, facets[f][j]) - 1;
              vindir[indv] = 1;
            }
          }
        }
      }
    }
  }

  // Transform the exterior vertex mask of zeros and ones into a vertex map
  // from old to new connectivities, and get the number of unique exterior
  // vertices, npts2
  E_Int npts2 = K_CONNECT::prefixSum(vindir);

  // Build new NODE connectivity
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts2,
                                       0, "NODE", false, api);
  FldArrayF* f2;
  K_ARRAY::getFromArray3(tpl, f2);

  #pragma omp parallel
  {
    E_Int indv;
    // Copy fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f.begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < npts; i++)
      {
        indv = vindir[i];
        if (indv > 0) f2p[indv-1] = fp[i];
      }
    }
  }

  if (boolIndir)
  {
    indicesVertices = K_NUMPY::buildNumpyArray(npts2, 1, 1, 0);
    indicesv = K_NUMPY::getNumpyPtrI(indicesVertices);
    
    #pragma omp parallel for
    for (E_Int i = 0; i < npts; i++)
    {
      E_Int indv = vindir[i];
      if (indv > 0) indicesv[indv-1] = indv;
    }

    PyList_Append(indices, indicesVertices); Py_DECREF(indicesVertices);
  }

  // Free memory
  RELEASESHAREDS(tpl, f2);
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  return tpl;
}

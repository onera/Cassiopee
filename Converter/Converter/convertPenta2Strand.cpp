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
#include "converter.h"
#include <unordered_map>

// Find top vertex from vertex
E_Int findTopVertex(E_Int vertex, FldArrayI& cm, std::vector< std::vector<E_Int> >& cVE)
{
  std::vector<E_Int>& elts = cVE[vertex-1];
  E_Int elt;
  for (size_t i = 0; i < elts.size(); i++)
  {
    elt = elts[i];
    if (cm(elt, 1) == vertex) return cm(elt,4);
    if (cm(elt, 3) == vertex) return cm(elt,6);
    if (cm(elt, 2) == vertex) return cm(elt,5);
  }
  return -1;
}

E_Int findBottomVertex(E_Int vertex, FldArrayI& cm, std::vector< std::vector<E_Int> >& cVE)
{
  std::vector<E_Int>& elts = cVE[vertex-1];
  E_Int elt;
  for (size_t i = 0; i < elts.size(); i++)
  {
    elt = elts[i];
    if (cm(elt, 4) == vertex) return cm(elt,1);
    if (cm(elt, 6) == vertex) return cm(elt,3);
    if (cm(elt, 5) == vertex) return cm(elt,2);
  }
  return -1;
}

E_Int findShellBottomVertex(E_Int vertex, FldArrayI& cm, std::vector< std::vector<E_Int> >& cVE, E_Int& k)
{
  E_Int vertexp;
  k = 0;
  while (vertex != -1) 
  {
    vertexp = vertex;
    vertex = findBottomVertex(vertexp, cm, cVE);
    k++;
  }
  return vertexp;
}

// ============================================================================
/* Convert a strand array to a penta mesh */
// ============================================================================
PyObject* K_CONVERTER::convertPenta2Strand(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, f,
                               nil, njl, nkl, cnl, eltType);

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertPenta2Strand: array must be unstructured PENTA.");
    return NULL; 
  }

  if (res == 2)
  {
    if (strcmp(eltType, "PENTA") != 0)
    {
      RELEASESHAREDU(array, f, cnl);
      PyErr_SetString(PyExc_TypeError, 
                      "convertPenta2Strand: array must be a single PENTA connectivity.");
      return NULL; 
    }
  }

  E_Int api = f->getApi(); if (api == 2) api = 3;

  // the strand grid is a TRI array with extra vertices
  FldArrayI& cm = *(cnl->getConnect(0));
  E_Int nelts = cm.getSize();
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();

  std::vector< std::vector<E_Int> > cVE(npts);
  K_CONNECT::connectEV2VE(cm, cVE);

  // Find nk
  E_Int nk = 0;
  E_Int vertex = 1; // premier vertex
  
  while (vertex != -1)
  {
    vertex = findTopVertex(vertex, cm, cVE);
    nk += 1;
  }
  vertex = 1;
  while (vertex != -1)
  {
    vertex = findBottomVertex(vertex, cm, cVE);
    nk += 1;
  }
  nk = nk-1;
  //printf("I detected %d layers\n", nk);
  
  E_Int nptsShell = npts / nk;
  std::unordered_map<E_Int, E_Int> mapv; // tag of bottom shell vertex -> 1
  std::unordered_map<E_Int, E_Int> mape; // tag of bottom shell triangles -> 1
  
  E_Int* kk = new E_Int [npts]; // k of global vertex
  E_Int* shellVert = new E_Int [npts]; // shell vertex of global vertex

  // Find bottom shell vertices
  E_Int k;
  for (E_Int i = 0; i < npts; i++)
  {
    vertex = findShellBottomVertex(i+1, cm, cVE, k);
    mapv[vertex] = 1;
    kk[i] = k;
    shellVert[i] = vertex;
    std::vector<E_Int>& elts = cVE[vertex-1];
    for (E_Int e = 0; e < elts.size(); e++) mape[elts[e]] = 1;
  }


  //printf("size of mapv = %ld, size of shell=%d\n", mapv.size(), nptsShell);
  //printf("size of mape = %ld\n", mape.size());
  E_Int nt = mape.size(); // nbre de triangles de la bottom shell

  // numerote le shell localement suivant la map
  E_Int* no = new E_Int [nptsShell]; // shell vertex -> local shell vertex
  E_Int c = 0;
  for (auto it=mapv.begin(); it != mapv.end(); it++) 
  {
    no[it->first-1] = c;
    //printf("no: %d %d\n", it->first-1, c);
    c++;
  }
  
  // Build output
  // maillage strand : 
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nt,
                                       "TRI", false, api);
  FldArrayF* f2; FldArrayI* cnl2;
  K_ARRAY::getFromArray3(tpl, f2, cnl2);
  FldArrayI& cm2 = *(cnl2->getConnect(0));

  // copie du field avec renumerotation
  E_Int ind;
  for (E_Int n = 1; n <= nfld; n++)
  {
    E_Float* fp = f->begin(n);
    E_Float* f2p = f2->begin(n);
    for (E_Int i = 0; i < npts; i++)
    {
      //printf("nodes: %d: %d %d\n", i, shellVert[i], kk[i]);
      ind = no[shellVert[i]-1]+(kk[i]-1)*nptsShell;
      f2p[ind] = fp[i];
    }
  }

  // copie de la connectivite
  E_Int i = 0;
  E_Int elt, v1, v2, v3, ind1, ind2, ind3;
  for (auto it = mape.begin(); it != mape.end(); it++)
  {
    elt = it->first; // shell element
    //printf("connect: %d: %d\n", i, elt);
    v1 = cm(elt, 1)-1; v2 = cm(elt, 2)-1; v3 = cm(elt, 3)-1;
    //printf("bottom vertices: %d %d %d\n", v1, v2, v3);
    
    ind1 = no[shellVert[v1]-1]+(kk[v1]-1)*nptsShell+1;
    ind2 = no[shellVert[v2]-1]+(kk[v2]-1)*nptsShell+1;
    ind3 = no[shellVert[v3]-1]+(kk[v3]-1)*nptsShell+1;
    
    //printf("bottom new vertices: %d %d %d\n", ind1, ind2, ind3);

    cm2(i, 1) = ind1;
    cm2(i, 2) = ind2;
    cm2(i, 3) = ind3;
    i++;
  }

  delete [] kk;
  delete [] shellVert;
  delete [] no;

  return tpl;
}
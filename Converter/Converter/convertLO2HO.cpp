/*    
    Copyright 2013-2019 Onera.

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
# include <unordered_map>

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

#define ETK E_LONG
#define EDGEINDEX(n1,n2,ind) \
    if (n1 < n2) k = (ETK)(n1-1)*nvertex+(ETK)(n2-1); \
    else k = (ETK)(n2-1)+nvertex*(ETK)(n1-1); \
    ind = map[k];
#define ADDEDGE(n1,n2) \
    if (n1 < n2) k = (ETK)(n1-1)*nvertex+(ETK)(n2-1); \
    else k = (ETK)(n2-1)+nvertex*(ETK)(n1-1); \
    it = map.find(k); \
    if (it == map.end()) { map[k] = compt; compt++; } \
    else map[k] = it->second;

// ============================================================================
/* Convert LO mesh to HO mesh */
// ============================================================================
PyObject* K_CONVERTER::convertLO2HO(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int mode;
  if (!PYPARSETUPLEI(args, "Ol", "Oi", &array, &mode)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray2(array, varString, 
                               f, ni, nj, nk, cn, eltType);

  if (res != 1 && res != 2)
  {
     PyErr_SetString(PyExc_TypeError, 
                     "convertLO2HO: array is invalid.");
     return NULL;
  }
  if (res == 1)
  {   
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertLO2HO: array must be unstructured.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, 4, "NGON") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "convertLO2HO: array must not be NGON.");
    return NULL;
  }
  
  // Caracteristiques de l'array input
  E_Int nelts = cn->getSize();
  E_Int nfld = f->getNfld();
  E_Int nvertex = f->getSize();
  E_Int api = f->getApi();
  PyObject* o = NULL;

  // BAR -> BAR_3
  if (K_STRING::cmp(eltType, 3, "BAR") == 0 && mode == 0)
  {
    E_Int nvertexHO = nvertex + nelts;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "BAR_3", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1, p2;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++)
        (*fo)(i,n) = (*f)(i,n);
      // ajout des sommets milieux a la fin
      for (E_Int i = 0; i < nelts; i++)
      {
        p1 = (*cn)(i,1)-1; p2 = (*cn)(i,2)-1;
        (*fo)(i+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
    }
    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      (*co)(i,1) = (*cn)(i,1);
      (*co)(i,2) = (*cn)(i,2);
      (*co)(i,3) = nvertex+i+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // TRI -> TRI_6
  else if (K_STRING::cmp(eltType, 3, "TRI") == 0 && mode == 0)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    E_Int compt = 0; std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3);
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n1,n3);
    }
    E_Int nedges = map.size();

    E_Int nvertexHO = nvertex + nedges;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "TRI_6", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1, p2,ind,n4,n5,n6;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<ETK,E_Int>& elt : map)
      {
        k = elt.first;
        ind = elt.second;
        p1 = k/nvertex;
        p2 = k - p1*nvertex;
        (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
    }
    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3);
      EDGEINDEX(n1,n2,n4);
      EDGEINDEX(n2,n3,n5);
      EDGEINDEX(n1,n3,n6);
    
      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4+nvertex+1;
      (*co)(i,5) = n5+nvertex+1; 
      (*co)(i,6) = n6+nvertex+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // QUAD -> QUAD_8
  else if (K_STRING::cmp(eltType, 4, "QUAD") == 0 && mode == 0)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    E_Int compt = 0; std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3,n4; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n3,n4);
      ADDEDGE(n1,n4);
    }
    E_Int nedges = map.size();

    E_Int nvertexHO = nvertex + nedges;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_8", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,ind,n5,n6,n7,n8;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<ETK,E_Int>& elt : map)
      {
        k = elt.first;
        ind = elt.second;
        p1 = k/nvertex;
        p2 = k - p1*nvertex;
        (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
    }
    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      EDGEINDEX(n1,n2,n5);
      EDGEINDEX(n2,n3,n6);
      EDGEINDEX(n3,n4,n7);
      EDGEINDEX(n1,n4,n8);
      
      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;
      (*co)(i,5) = n5+nvertex+1;
      (*co)(i,6) = n6+nvertex+1;
      (*co)(i,7) = n7+nvertex+1;
      (*co)(i,8) = n8+nvertex+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // QUAD -> QUAD_9
  else if (K_STRING::cmp(eltType, 4, "QUAD") == 0 && mode == 1)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    E_Int compt = 0; std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3,n4; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n3,n4);
      ADDEDGE(n1,n4);
    }
    E_Int nedges = map.size();

    E_Int nvertexHO = nvertex + nedges + nelts;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_9", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,ind,n5,n6,n7,n8;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<ETK,E_Int>& elt : map)
      {
        k = elt.first;
        ind = elt.second;
        p1 = k/nvertex;
        p2 = k - p1*nvertex;
        (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
      // ajout pour les centres
      for (E_Int i = 0; i < nelts; i++)
      {
        n1 = (*cn)(i,1)-1; n2 = (*cn)(i,2)-1; n3 = (*cn)(i,3)-1; n4 = (*cn)(i,4)-1;
        (*fo)(i+nvertex+nedges,n) = 0.25*( (*f)(n1,n)+ (*f)(n2,n) + (*f)(n3,n) + (*f)(n4,n) );
      }
    }
    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      EDGEINDEX(n1,n2,n5);
      EDGEINDEX(n2,n3,n6);
      EDGEINDEX(n3,n4,n7);
      EDGEINDEX(n1,n4,n8);
      
      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;
      (*co)(i,5) = n5+nvertex+1;
      (*co)(i,6) = n6+nvertex+1;
      (*co)(i,7) = n7+nvertex+1;
      (*co)(i,8) = n8+nvertex+1;
      (*co)(i,9) = i+nvertex+nedges+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // TETRA -> TETRA_10
  else if (K_STRING::cmp(eltType, 5, "TETRA") == 0 && mode == 0)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    E_Int compt = 0; std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3,n4; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n1,n3);
      ADDEDGE(n1,n4);
      ADDEDGE(n2,n4);
      ADDEDGE(n3,n4);
    }
    E_Int nedges = map.size();

    E_Int nvertexHO = nvertex + nedges;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "TETRA_10", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,ind,n5,n6,n7,n8,n9,n10;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<ETK,E_Int>& elt : map)
      {
        k = elt.first;
        ind = elt.second;
        p1 = k/nvertex;
        p2 = k - p1*nvertex;
        (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
    }
    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      EDGEINDEX(n1,n2,n5);
      EDGEINDEX(n2,n3,n6);
      EDGEINDEX(n1,n3,n7);
      EDGEINDEX(n1,n4,n8);
      EDGEINDEX(n2,n4,n9);
      EDGEINDEX(n3,n4,n10);
      
      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;
      (*co)(i,5) = n5+nvertex+1;
      (*co)(i,6) = n6+nvertex+1;
      (*co)(i,7) = n7+nvertex+1;
      (*co)(i,8) = n8+nvertex+1;
      (*co)(i,9) = n9+nvertex+1;
      (*co)(i,10) = n10+nvertex+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  RELEASESHAREDB(res, array, f, cn);
  
  return o;
}

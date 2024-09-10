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
# include <unordered_map>

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

#define ETK E_LONG
#define EDGEINDEX(n1,n2,ind) \
    if (n1 < n2) { k = (ETK)(n1-1)*nvertex+(ETK)(n2-1); }\
    else { k = (ETK)(n1-1) + nvertex*(ETK)(n2-1); } \
    ind = map[k];
#define ADDEDGE(n1,n2) \
    if (n1 < n2)   k = (ETK)(n1-1)*nvertex+(ETK)(n2-1); \
    else k = (ETK)(n1-1)+nvertex*(ETK)(n2-1); \
    it = map.find(k); \
    if (it == map.end()) { map[k] = compt; compt++; } 
    
#define ADDFACE(n1,n2,n3,n4) \
    if (n1 < n2) { l1 = n1; h1 = n2; } \
    else { l1 = n2; h1 = n1; } \
    if (n3 < n4) { l2 = n3; h2 = n4; } \
    else { l2 = n4; h2 = n3; } \
    if (l1 < l2) { ls = l1; md1 = l2; } \
    else { ls = l2; md1 = l1; } \
    if (h1 > h2) { hs = h1; md2 = h2; } \
    else { hs = h2; md2 = h1; } \
    if (md1 < md2) k = (ETK)(ls-1)+(ETK)(md1-1)*nvertex+(ETK)(md2-1)*nvertex*nvertex+(ETK)(hs-1)*nvertex*nvertex*nvertex; \
    else k = (ETK)(ls-1)+(ETK)(md2-1)*nvertex+(ETK)(md1-1)*nvertex*nvertex+(ETK)(hs-1)*nvertex*nvertex*nvertex; \
    it = map2.find(k); \
    if (it == map2.end()) { map2[k] = compt2; compt2++; } \
    else map2[k] = it->second;
#define FACEINDEX(n1,n2,n3,n4,ind) \
    if (n1 < n2) { l1 = n1; h1 = n2; } \
    else { l1 = n2; h1 = n1; } \
    if (n3 < n4) { l2 = n3; h2 = n4; } \
    else { l2 = n4; h2 = n3; } \
    if (l1 < l2) { ls = l1; md1 = l2; } \
    else { ls = l2; md1 = l1; } \
    if (h1 > h2) { hs = h1; md2 = h2; } \
    else { hs = h2; md2 = h1; } \
    if (md1 < md2) k = (ETK)(ls-1)+(ETK)(md1-1)*nvertex+(ETK)(md2-1)*nvertex*nvertex+(ETK)(hs-1)*nvertex*nvertex*nvertex; \
    else k = (ETK)(ls-1)+(ETK)(md2-1)*nvertex+(ETK)(md1-1)*nvertex*nvertex+(ETK)(hs-1)*nvertex*nvertex*nvertex; \
    ind = map2[k];

// ============================================================================
/* Convert LO mesh to HO mesh */
// ============================================================================
PyObject* K_CONVERTER::convertLO2HO(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int mode, order;
  order = 2;
  if (!PYPARSETUPLE_(args, O_ I_ "|" I_, &array, &mode, &order)) return NULL;

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
  if ( (order < 2) || (order > 4) )
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "convertLO2HO: order must be 2, 3 or 4.");
    return NULL; 
  }

  constexpr const double onethird = 1./3.;
  constexpr const double twothird = 2./3.;
  constexpr const double onequart = 1./4.;
  constexpr const double onehalf  = 1./2.;
  constexpr const double threequart = 3./4.;

  // Caracteristiques de l'array input
  E_Int nelts = cn->getSize();
  E_Int nfld = f->getNfld();
  E_Int nvertex = f->getSize();
  E_Int api = f->getApi();
  PyObject* o = NULL;

  E_Int ls,md1,md2,hs,l1,l2,h1,h2;

  // BAR -> BAR_3
  if (K_STRING::cmp(eltType, 3, "BAR") == 0 && mode == 0)
  {

    E_Int nvertexHO = nvertex + (order-1)*nelts;
    E_Int neltsHO = nelts;
    if (order == 2 )
    {
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "BAR_3", false, 0, 0, 0, api);
    }
    else if (order == 3)
    {
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "BAR_4", false, 0, 0, 0, api);      
    }
    else // order == 4
    {
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "BAR_5", false, 0, 0, 0, api);
    }
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
      if (order == 2)
      {
        for (E_Int i = 0; i < nelts; i++)
        {
          p1 = (*cn)(i,1)-1; p2 = (*cn)(i,2)-1;
          (*fo)(i+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
        }
      }
      else if (order == 3)
      {
        for (E_Int i = 0; i < nelts; i++)
        {
          p1 = (*cn)(i,1)-1; p2 = (*cn)(i,2)-1;
          (*fo)(2*i+0+nvertex,n) = twothird*(*f)(p1,n) + onethird*(*f)(p2,n);
          (*fo)(2*i+1+nvertex,n) = onethird*(*f)(p1,n) + twothird*(*f)(p2,n);
        }
      }
      else //if (order == 4)
      {
        for (E_Int i = 0; i < nelts; i++)
        {
          p1 = (*cn)(i,1)-1; p2 = (*cn)(i,2)-1;
          (*fo)(3*i+0+nvertex,n) = threequart*(*f)(p1,n) + onequart * (*f)(p2,n);
          (*fo)(3*i+1+nvertex,n) = onehalf   *( (*f)(p1,n)+ (*f)(p2,n) );
          (*fo)(3*i+2+nvertex,n) = onequart  *(*f)(p1,n) + threequart*(*f)(p2,n);
        }
      }
    }
    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      (*co)(i,1) = (*cn)(i,1);
      (*co)(i,2) = (*cn)(i,2);
      for ( E_Int ord = 1; ord < order; ++ord )
        (*co)(i,2+ord) = nvertex+(order-1)*i+ord;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // TRI -> TRI_6 ou TRI_9 ou TRI_12
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

    E_Int nvertexHO = nvertex + (order-1)*nedges;
    E_Int neltsHO = nelts;
    if (order == 2)
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "TRI_6", false, 0, 0, 0, api);
    else if (order == 3)
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "TRI_9", false, 0, 0, 0, api);
    else // order == 4
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "TRI_12", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1, p2,ind,n4,n5,n6;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      if (order == 2)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
        }
      }
      else if (order == 3)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(2*ind+0+nvertex,n) = twothird * (*f)(p1,n) + onethird * (*f)(p2,n);
          (*fo)(2*ind+1+nvertex,n) = onethird * (*f)(p1,n) + twothird * (*f)(p2,n);
        }
      }
      else // if (order == 4)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(3*ind+0+nvertex,n) = threequart*(*f)(p1,n) + onequart * (*f)(p2,n);
          (*fo)(3*ind+1+nvertex,n) = onehalf   *( (*f)(p1,n)+ (*f)(p2,n) );
          (*fo)(3*ind+2+nvertex,n) = onequart  *(*f)(p1,n) + threequart*(*f)(p2,n);
        }
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

      if (n1<n2) // Si n1 < n2, alors arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+ord) = (order-1)*n4 + nvertex + ord;
      else // Arete a l'envers par rapport a la convention cgns :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+ord) = (order-1)*n4 + nvertex + (order-ord);

      if (n2 < n3) // Si n2 < n3, arete dans le bon sens :      
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+(order-1)+ord) = (order-1)*n5 + nvertex + ord;
      else
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+(order-1)+ord) = (order-1)*n5 + nvertex + (order-ord);
      if (n3 < n1) // Si n3 < n1, arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+2*(order-1)+ord) = (order-1)*n6 + nvertex + ord;
      else // Arete pas dans le bon sens
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+2*(order-1)+ord) = (order-1)*n6 + nvertex + (order-ord);
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // QUAD -> QUAD_8 ou QUAD_12 ( et QUAD_P4_16 mais pour sa visu, on verra plus tard, car ce couillon a le meme nombre de sommets que QUAD_16 => pas cool pour le shader ! )
  // Pour Quad_P4_16 il faudra donc faire un shader a part...
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

    E_Int nvertexHO = nvertex + (order-1)*nedges;
    E_Int neltsHO = nelts;
    if (order == 2)
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_8", false, 0, 0, 0, api);
    else if (order == 3)
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_12", false, 0, 0, 0, api);
    else // order == 4
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_P4_16", false, 0, 0, 0, api);

    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,ind,n5,n6,n7,n8;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      if (order == 2)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
           p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
        }
      }
      else if (order == 3)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(2*ind+0+nvertex,n) = twothird * (*f)(p1,n) + onethird * (*f)(p2,n);
          (*fo)(2*ind+1+nvertex,n) = onethird * (*f)(p1,n) + twothird * (*f)(p2,n);
        }
      }
      else // order == 4
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(3*ind+0+nvertex,n) = threequart*(*f)(p1,n) + onequart * (*f)(p2,n);
          (*fo)(3*ind+1+nvertex,n) = onehalf * ((*f)(p1,n) + (*f)(p2,n) );
          (*fo)(3*ind+2+nvertex,n) = onequart  *(*f)(p1,n) + threequart*(*f)(p2,n);
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

      if (n1<n2) // Si n1 < n2, alors arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+ord) = (order-1)*n5 + nvertex + ord;
      else // Arete a l'envers par rapport a la convention cgns :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+ord) = (order-1)*n5 + nvertex + (order-ord);

      if (n2 < n3) // Si n2 < n3, arete dans le bon sens :      
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+(order-1)+ord) = (order-1)*n6 + nvertex + ord;
      else
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+(order-1)+ord) = (order-1)*n6 + nvertex + (order-ord);

      if (n3 < n4) // Si n3 < n4, arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+2*(order-1)+ord) = (order-1)*n7 + nvertex + ord;
      else // Arete pas dans le bon sens
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+2*(order-1)+ord) = (order-1)*n7 + nvertex + (order-ord);

      if (n4 < n1) // Si n4 < n1, arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+3*(order-1)+ord) = (order-1)*n8 + nvertex + ord;
      else // Arete pas dans le bon sens
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+3*(order-1)+ord) = (order-1)*n8 + nvertex + (order-ord);

    }
    RELEASESHAREDU(o, fo, co); 
  }
  // TRI -> TRI_10 ou TRI_15
  else if (K_STRING::cmp(eltType, 4, "TRI") == 0 && mode == 1)
  {
    if (order < 3)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "convertLO2HO: order must be 3 or 4 for triangle with mode 1.");
      return NULL; 
    }
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    E_Int compt = 0; 
    std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3);
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n1,n3);
    }
    E_Int nedges = map.size();
    E_Int nvertexHO = nvertex + (order-1)*nedges + (order == 3 ? 1 : 3) * nelts;
    E_Int neltsHO = nelts;
    if (order == 3)
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "TRI_10", false, 0, 0, 0, api);
    else // order == 4
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "TRI_15", false, 0, 0, 0, api);

    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,ind,n4,n5,n6;//,n7,n8;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      if (order == 3)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(2*ind+0+nvertex,n) = twothird * (*f)(p1,n) + onethird * (*f)(p2,n);
          (*fo)(2*ind+1+nvertex,n) = onethird * (*f)(p1,n) + twothird * (*f)(p2,n);
        }
      }
      else // if (order == 4)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(3*ind+0+nvertex,n) = threequart*(*f)(p1,n) + onequart * (*f)(p2,n);
          (*fo)(3*ind+1+nvertex,n) = onehalf   *( (*f)(p1,n)+ (*f)(p2,n) );
          (*fo)(3*ind+2+nvertex,n) = onequart  *(*f)(p1,n) + threequart*(*f)(p2,n);
        }
      }
      // ajout pour les centres
      if (order == 3)
      {
        E_Int beg_centers = 2*nedges;
        for (E_Int i = 0; i < nelts; i++)
        {
          n1 = (*cn)(i,1)-1; n2 = (*cn)(i,2)-1; n3 = (*cn)(i,3)-1;
          (*fo)(i+nvertex+beg_centers,n) = (onethird)*( (*f)(n1,n)+ (*f)(n2,n) + (*f)(n3,n) );
        }
      }
      else // order == 4
      {
        E_Int beg_centers = 3*nedges;
        for (E_Int i = 0; i < nelts; i++)
        {
          n1 = (*cn)(i,1)-1; n2 = (*cn)(i,2)-1; n3 = (*cn)(i,3)-1;
          (*fo)(3*i+0+nvertex+beg_centers,n) = onehalf*(*f)(n1,n) + onequart*((*f)(n2,n) + (*f)(n3,n) );
          (*fo)(3*i+1+nvertex+beg_centers,n) = onehalf*(*f)(n2,n) + onequart*((*f)(n3,n) + (*f)(n1,n) );
          (*fo)(3*i+2+nvertex+beg_centers,n) = onehalf*(*f)(n3,n) + onequart*((*f)(n1,n) + (*f)(n2,n) );
        }
      }
    }
    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3);

      EDGEINDEX(n1,n2,n4);
      EDGEINDEX(n2,n3,n5);
      EDGEINDEX(n3,n1,n6);

      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;

      if (n1<n2) // Si n1 < n2, alors arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+ord) = (order-1)*n4 + nvertex + ord;
      else // Arete a l'envers par rapport a la convention cgns :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+ord) = (order-1)*n4 + nvertex + (order-ord);

      if (n2 < n3) // Si n2 < n3, arete dans le bon sens :      
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+(order-1)+ord) = (order-1)*n5 + nvertex + ord;
      else
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+(order-1)+ord) = (order-1)*n5 + nvertex + (order-ord);

      if (n3 < n1) // Si n3 < n1, arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+2*(order-1)+ord) = (order-1)*n6 + nvertex + ord;
      else // Arete pas dans le bon sens
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,3+2*(order-1)+ord) = (order-1)*n6 + nvertex + (order-ord);

      if (order == 3)
        (*co)(i,10) = i+nvertex+2*nedges+1;
      else // order == 4
      {
        (*co)(i,10) = 3*i+0+nvertex+2*nedges+1;
        (*co)(i,11) = 3*i+1+nvertex+2*nedges+1;
        (*co)(i,12) = 3*i+2+nvertex+2*nedges+1;
      }
    }
    RELEASESHAREDU(o, fo, co);
  }
  // QUAD -> QUAD_9 ou QUAD_16 ou QUAD_25
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

    E_Int nvertexHO = nvertex + (order-1)*nedges + (order-1)*(order-1)*nelts;
    E_Int neltsHO = nelts;
    if (order == 2)
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_9", false, 0, 0, 0, api);
    else if (order == 3)
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_16", false, 0, 0, 0, api);
    else // order == 4
      o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "QUAD_25", false, 0, 0, 0, api);

    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,ind,n5,n6,n7,n8;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      if (order == 2)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
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
          (*fo)(i+nvertex+nedges,n) = onequart*( (*f)(n1,n)+ (*f)(n2,n) + (*f)(n3,n) + (*f)(n4,n) );
        }
      }
      else if (order == 3)
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(2*ind+0+nvertex,n) = twothird * (*f)(p1,n) + onethird * (*f)(p2,n);
          (*fo)(2*ind+1+nvertex,n) = onethird * (*f)(p1,n) + twothird * (*f)(p2,n);
        }
        // ajout pour les centres
        for (E_Int i = 0; i < nelts; i++)
        {
          n1 = (*cn)(i,1)-1; n2 = (*cn)(i,2)-1; n3 = (*cn)(i,3)-1; n4 = (*cn)(i,4)-1;
          (*fo)(4*i+0+nvertex+2*nedges,n) = (4./9.)*(*f)(n1,n) + (2./9.) * ( (*f)(n2,n) + (*f)(n4,n) ) + (1./9.)*(*f)(n3,n);
          (*fo)(4*i+1+nvertex+2*nedges,n) = (4./9.)*(*f)(n2,n) + (2./9.) * ( (*f)(n3,n) + (*f)(n1,n) ) + (1./9.)*(*f)(n4,n);
          (*fo)(4*i+2+nvertex+2*nedges,n) = (4./9.)*(*f)(n3,n) + (2./9.) * ( (*f)(n4,n) + (*f)(n2,n) ) + (1./9.)*(*f)(n1,n);
          (*fo)(4*i+3+nvertex+2*nedges,n) = (4./9.)*(*f)(n4,n) + (2./9.) * ( (*f)(n1,n) + (*f)(n3,n) ) + (1./9.)*(*f)(n2,n);
        }
      }
      else // order == 4
      {
        for (const std::pair<const ETK,E_Int>& elt : map)
        {
          k = elt.first;
          ind = elt.second;
          p1 = k/nvertex;
          p2 = k - p1*nvertex;
          (*fo)(3*ind+0+nvertex,n) = threequart*(*f)(p1,n) + onequart * (*f)(p2,n);
          (*fo)(3*ind+1+nvertex,n) = onehalf * ((*f)(p1,n) + (*f)(p2,n) );
          (*fo)(3*ind+2+nvertex,n) = onequart  *(*f)(p1,n) + threequart*(*f)(p2,n);
        }
        // ajout pour les centres
        E_Int beg_centers = nvertex + 3*nedges;
        for (E_Int i = 0; i < nelts; i++)
        {
          n1 = (*cn)(i,1)-1; n2 = (*cn)(i,2)-1; n3 = (*cn)(i,3)-1; n4 = (*cn)(i,4)-1;
          (*fo)(9*i+0+beg_centers,n) = (9./16.)*(*f)(n1,n) + (3./16.) * ( (*f)(n2,n) + (*f)(n4,n) ) + (1./16.)*(*f)(n3,n);
          (*fo)(9*i+1+beg_centers,n) = (3./8.0)*( (*f)(n1,n) + (*f)(n2,n) ) + (1./8.0) * ( (*f)(n4,n) + (*f)(n3,n) );
          (*fo)(9*i+2+beg_centers,n) = (9./16.)*(*f)(n2,n) + (3./16.) * ( (*f)(n3,n) + (*f)(n1,n) ) + (1./16.)*(*f)(n4,n);
          (*fo)(9*i+3+beg_centers,n) = (3./8.0)*( (*f)(n2,n) + (*f)(n3,n) ) + (1./8.0) * ( (*f)(n1,n) + (*f)(n4,n) );
          (*fo)(9*i+4+beg_centers,n) = (9./16.)*(*f)(n3,n) + (3./16.) * ( (*f)(n4,n) + (*f)(n2,n) ) + (1./16.)*(*f)(n1,n);
          (*fo)(9*i+5+beg_centers,n) = (3./8.0)*( (*f)(n3,n) + (*f)(n4,n) ) + (1./8.0) * ( (*f)(n2,n) + (*f)(n1,n) );
          (*fo)(9*i+6+beg_centers,n) = (9./16.)*(*f)(n4,n) + (3./16.) * ( (*f)(n1,n) + (*f)(n3,n) ) + (1./16.)*(*f)(n2,n);
          (*fo)(9*i+7+beg_centers,n) = (3./8.0)*( (*f)(n4,n) + (*f)(n1,n) ) + (1./8.0) * ( (*f)(n3,n) + (*f)(n2,n) );
          (*fo)(9*i+8+beg_centers,n) = (1./4.0)*( (*f)(n1,n) + (*f)(n2,n) + (*f)(n3,n) + (*f)(n4,n) );
        }
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

      E_Int order2 = (order-1)*(order-1);
      
      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;

      if (n1<n2) // Si n1 < n2, alors arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+ord) = (order-1)*n5 + nvertex + ord;
      else // Arete a l'envers par rapport a la convention cgns :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+ord) = (order-1)*n5 + nvertex + (order-ord);

      if (n2 < n3) // Si n2 < n3, arete dans le bon sens :      
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+(order-1)+ord) = (order-1)*n6 + nvertex + ord;
      else
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+(order-1)+ord) = (order-1)*n6 + nvertex + (order-ord);

      if (n3 < n4) // Si n3 < n4, arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+2*(order-1)+ord) = (order-1)*n7 + nvertex + ord;
      else // Arete pas dans le bon sens
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+2*(order-1)+ord) = (order-1)*n7 + nvertex + (order-ord);

      if (n4 < n1) // Si n4 < n1, arete dans le bon sens :
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+3*(order-1)+ord) = (order-1)*n8 + nvertex + ord;
      else // Arete pas dans le bon sens
        for (E_Int ord = 1; ord < order; ++ord )
          (*co)(i,4+3*(order-1)+ord) = (order-1)*n8 + nvertex + (order-ord);


      for (E_Int ord = 1; ord <= order2; ++ord)
        (*co)(i,4+4*(order-1)+ord) = order2*i + nvertex + (order-1)*nedges + ord;
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
      for (const std::pair<const ETK,E_Int>& elt : map)
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
  // HEXA -> HEXA_20
  else if (K_STRING::cmp(eltType, 4, "HEXA") == 0 && mode == 0)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    E_Int compt = 0; std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3,n4,n5,n6,n7,n8; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5); n6 = (*cn)(i,6); n7 = (*cn)(i,7); n8 = (*cn)(i,8);
      
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n3,n4);
      ADDEDGE(n1,n4);
      ADDEDGE(n1,n5);
      ADDEDGE(n2,n6);
      ADDEDGE(n3,n7);
      ADDEDGE(n4,n8);
      ADDEDGE(n5,n6);
      ADDEDGE(n6,n7);
      ADDEDGE(n7,n8);
      ADDEDGE(n5,n8);
    }
    E_Int nedges = map.size();

    E_Int nvertexHO = nvertex + nedges;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "HEXA_20", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,ind,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<const ETK,E_Int>& elt : map)
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
      n5 = (*cn)(i,5); n6 = (*cn)(i,6); n7 = (*cn)(i,7); n8 = (*cn)(i,8);
      
      EDGEINDEX(n1,n2,n9);
      EDGEINDEX(n2,n3,n10);
      EDGEINDEX(n3,n4,n11);
      EDGEINDEX(n1,n4,n12);
      EDGEINDEX(n1,n5,n13);
      EDGEINDEX(n2,n6,n14);
      EDGEINDEX(n3,n7,n15);
      EDGEINDEX(n4,n8,n16);
      EDGEINDEX(n5,n6,n17);
      EDGEINDEX(n6,n7,n18);
      EDGEINDEX(n7,n8,n19);
      EDGEINDEX(n5,n8,n20);
      
      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;
      (*co)(i,5) = n5;
      (*co)(i,6) = n6;
      (*co)(i,7) = n7;
      (*co)(i,8) = n8;
      (*co)(i,9) = n9+nvertex+1;
      (*co)(i,10) = n10+nvertex+1;
      (*co)(i,11) = n11+nvertex+1;
      (*co)(i,12) = n12+nvertex+1;
      (*co)(i,13) = n13+nvertex+1;
      (*co)(i,14) = n14+nvertex+1;
      (*co)(i,15) = n15+nvertex+1;
      (*co)(i,16) = n16+nvertex+1;
      (*co)(i,17) = n17+nvertex+1;
      (*co)(i,18) = n18+nvertex+1;
      (*co)(i,19) = n19+nvertex+1;
      (*co)(i,20) = n20+nvertex+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // HEXA -> HEXA_27
  else if (K_STRING::cmp(eltType, 4, "HEXA") == 0 && mode == 1)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    std::unordered_map<ETK, E_Int> map2;
    E_Int compt = 0; E_Int compt2 = 0; 
    std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3,n4,n5,n6,n7,n8; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5); n6 = (*cn)(i,6); n7 = (*cn)(i,7); n8 = (*cn)(i,8);
      
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n3,n4);
      ADDEDGE(n1,n4);
      ADDEDGE(n1,n5);
      ADDEDGE(n2,n6);
      ADDEDGE(n3,n7);
      ADDEDGE(n4,n8);
      ADDEDGE(n5,n6);
      ADDEDGE(n6,n7);
      ADDEDGE(n7,n8);
      ADDEDGE(n5,n8);
    }
    E_Int nedges = map.size();
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5); n6 = (*cn)(i,6); n7 = (*cn)(i,7); n8 = (*cn)(i,8);
      ADDFACE(n1,n2,n3,n4);
      ADDFACE(n1,n2,n6,n5);
      ADDFACE(n2,n3,n7,n6);
      ADDFACE(n3,n4,n8,n7);
      ADDFACE(n1,n4,n8,n5);
      ADDFACE(n5,n6,n7,n8);
    }
    E_Int nfaces = map2.size();

    E_Int nvertexHO = nvertex + nedges + nfaces + nelts;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "HEXA_27", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,p3,p4,ind,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<const ETK,E_Int>& elt : map)
      {
        k = elt.first;
        ind = elt.second;
        p1 = k/nvertex;
        p2 = k - p1*nvertex;
        (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
    }

    // faces + pt interieur
    for (E_Int n = 1; n <= nfld; n++)
    {
      // ajout pour chaque face
      for (const std::pair<const ETK,E_Int>& elt : map2)
      {
        k = elt.first;
        ind = elt.second;
        p4 = k/(nvertex*nvertex*nvertex);
        p3 = k - p4*nvertex*nvertex*nvertex;
        p3 = p3/(nvertex*nvertex);
        p2 = k - p3*nvertex*nvertex - p4*nvertex*nvertex*nvertex;
        p2 = p2/nvertex;
        p1 = k - p2*nvertex - p3*nvertex*nvertex - p4*nvertex*nvertex*nvertex;
        (*fo)(ind+nvertex+nedges,n) = 0.25*( (*f)(p1,n)+ (*f)(p2,n) + (*f)(p3,n)+ (*f)(p4,n) );
      }

      // pt interieur
      for (E_Int i = 0; i < nelts; i++)
      {
        n1 = (*cn)(i,1)-1; n2 = (*cn)(i,2)-1; n3 = (*cn)(i,3)-1; n4 = (*cn)(i,4)-1;
        n5 = (*cn)(i,5)-1; n6 = (*cn)(i,6)-1; n7 = (*cn)(i,7)-1; n8 = (*cn)(i,8)-1;
        (*fo)(i+nvertex+nedges+nfaces,n) = 0.125*( (*f)(n1,n)+ (*f)(n2,n) + (*f)(n3,n)+ (*f)(n4,n) + 
                                                  (*f)(n5,n)+ (*f)(n6,n) + (*f)(n7,n)+ (*f)(n8,n)  );
      }
    }

    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5); n6 = (*cn)(i,6); n7 = (*cn)(i,7); n8 = (*cn)(i,8);
      
      EDGEINDEX(n1,n2,n9);
      EDGEINDEX(n2,n3,n10);
      EDGEINDEX(n3,n4,n11);
      EDGEINDEX(n1,n4,n12);
      EDGEINDEX(n1,n5,n13);
      EDGEINDEX(n2,n6,n14);
      EDGEINDEX(n3,n7,n15);
      EDGEINDEX(n4,n8,n16);
      EDGEINDEX(n5,n6,n17);
      EDGEINDEX(n6,n7,n18);
      EDGEINDEX(n7,n8,n19);
      EDGEINDEX(n5,n8,n20);
      FACEINDEX(n1,n2,n3,n4,n21);
      FACEINDEX(n1,n2,n6,n5,n22);
      FACEINDEX(n2,n3,n7,n6,n23);
      FACEINDEX(n3,n4,n8,n7,n24);
      FACEINDEX(n1,n5,n8,n4,n25);
      FACEINDEX(n5,n6,n7,n8,n26);

      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;
      (*co)(i,5) = n5;
      (*co)(i,6) = n6;
      (*co)(i,7) = n7;
      (*co)(i,8) = n8;
      (*co)(i,9) = n9+nvertex+1;
      (*co)(i,10) = n10+nvertex+1;
      (*co)(i,11) = n11+nvertex+1;
      (*co)(i,12) = n12+nvertex+1;
      (*co)(i,13) = n13+nvertex+1;
      (*co)(i,14) = n14+nvertex+1;
      (*co)(i,15) = n15+nvertex+1;
      (*co)(i,16) = n16+nvertex+1;
      (*co)(i,17) = n17+nvertex+1;
      (*co)(i,18) = n18+nvertex+1;
      (*co)(i,19) = n19+nvertex+1;
      (*co)(i,20) = n20+nvertex+1;
      (*co)(i,21) = n21+nvertex+nedges+1;
      (*co)(i,22) = n22+nvertex+nedges+1;
      (*co)(i,23) = n23+nvertex+nedges+1;
      (*co)(i,24) = n24+nvertex+nedges+1;
      (*co)(i,25) = n25+nvertex+nedges+1;
      (*co)(i,26) = n26+nvertex+nedges+1;
      (*co)(i,27) = i+nvertex+nedges+nfaces+1;
    }
  
    RELEASESHAREDU(o, fo, co); 
  }
  // PENTA -> PENTA_18
  else if (K_STRING::cmp(eltType, 5, "PENTA") == 0 && mode == 1)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    std::unordered_map<ETK, E_Int> map2;
    E_Int compt = 0; E_Int compt2 = 0; 
    std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3,n4,n5,n6; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5); n6 = (*cn)(i,6);
      
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n1,n3);
      ADDEDGE(n1,n4);
      ADDEDGE(n2,n5);
      ADDEDGE(n3,n6);
      ADDEDGE(n4,n5);
      ADDEDGE(n5,n6);
      ADDEDGE(n4,n6);
    }
    E_Int nedges = map.size();
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5); n6 = (*cn)(i,6); 
      ADDFACE(n1,n2,n5,n4);
      ADDFACE(n2,n3,n6,n5);
      ADDFACE(n1,n3,n6,n4);
    }
    E_Int nfaces = map2.size();

    E_Int nvertexHO = nvertex + nedges + nfaces;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "PENTA_18", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,p3,p4,ind,n7,n8,n9,n10,n11,n12,n13,n14,n15,n16,n17,n18;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<const ETK,E_Int>& elt : map)
      {
        k = elt.first;
        ind = elt.second;
        p1 = k/nvertex;
        p2 = k - p1*nvertex;
        (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
    }

    // faces 
    for (E_Int n = 1; n <= nfld; n++)
    {
      // ajout pour chaque face
      for (const std::pair<const ETK,E_Int>& elt : map2)
      {
        k = elt.first;
        ind = elt.second;
        p4 = k/(nvertex*nvertex*nvertex);
        p3 = k - p4*nvertex*nvertex*nvertex;
        p3 = p3/(nvertex*nvertex);
        p2 = k - p3*nvertex*nvertex - p4*nvertex*nvertex*nvertex;
        p2 = p2/nvertex;
        p1 = k - p2*nvertex - p3*nvertex*nvertex - p4*nvertex*nvertex*nvertex;
        (*fo)(ind+nvertex+nedges,n) = 0.25*( (*f)(p1,n)+ (*f)(p2,n) + (*f)(p3,n)+ (*f)(p4,n) );
      }
    }

    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5); n6 = (*cn)(i,6); 
      
      EDGEINDEX(n1,n2,n7);
      EDGEINDEX(n2,n3,n8);
      EDGEINDEX(n1,n3,n9);
      EDGEINDEX(n1,n4,n10);
      EDGEINDEX(n2,n5,n11);
      EDGEINDEX(n3,n6,n12);
      EDGEINDEX(n4,n5,n13);
      EDGEINDEX(n5,n6,n14);
      EDGEINDEX(n4,n6,n15);
      
      FACEINDEX(n1,n2,n5,n4,n16);
      FACEINDEX(n3,n2,n5,n6,n17);
      FACEINDEX(n1,n3,n6,n4,n18);

      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;
      (*co)(i,5) = n5;
      (*co)(i,6) = n6;
      (*co)(i,7) = n7+nvertex+1;
      (*co)(i,8) = n8+nvertex+1;
      (*co)(i,9) = n9+nvertex+1;
      (*co)(i,10) = n10+nvertex+1;
      (*co)(i,11) = n11+nvertex+1;
      (*co)(i,12) = n12+nvertex+1;
      (*co)(i,13) = n13+nvertex+1;
      (*co)(i,14) = n14+nvertex+1;
      (*co)(i,15) = n15+nvertex+1;

      (*co)(i,16) = n16+nvertex+nedges+1;
      (*co)(i,17) = n17+nvertex+nedges+1;
      (*co)(i,18) = n18+nvertex+nedges+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }
  // PYRA -> PYRA_14
  else if (K_STRING::cmp(eltType, 4, "PYRA") == 0 && mode == 1)
  {
    // compte les edges, cree la map des cles
    std::unordered_map<ETK, E_Int> map;
    std::unordered_map<ETK, E_Int> map2;
    E_Int compt = 0; E_Int compt2 = 0; 
    std::unordered_map<ETK,E_Int>::iterator it;
    E_Int n1,n2,n3,n4,n5; ETK k;
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5);
      
      ADDEDGE(n1,n2);
      ADDEDGE(n2,n3);
      ADDEDGE(n3,n4);
      ADDEDGE(n4,n1);
      ADDEDGE(n1,n5);
      ADDEDGE(n2,n5);
      ADDEDGE(n3,n5);
      ADDEDGE(n4,n5);
    }
    E_Int nedges = map.size();
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5);  
      ADDFACE(n1,n2,n3,n4);
    }
    E_Int nfaces = map2.size();

    E_Int nvertexHO = nvertex + nedges + nfaces;
    E_Int neltsHO = nelts;
    o = K_ARRAY::buildArray2(nfld, varString, nvertexHO, neltsHO, -1, "PYRA_14", false, 0, 0, 0, api);
    FldArrayF* fo; FldArrayI* co;
    K_ARRAY::getFromArray2(o, fo, co);
    E_Int p1,p2,p3,p4,ind,n6,n7,n8,n9,n10,n11,n12,n13,n14;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      // reprise des sommets LO
      for (E_Int i = 0; i < nvertex; i++) (*fo)(i,n) = (*f)(i,n);
      // ajout pour chaque edge
      for (const std::pair<const ETK,E_Int>& elt : map)
      {
        k = elt.first;
        ind = elt.second;
        p1 = k/nvertex;
        p2 = k - p1*nvertex;
        (*fo)(ind+nvertex,n) = 0.5*( (*f)(p1,n)+ (*f)(p2,n) );
      }
    }

    // faces 
    for (E_Int n = 1; n <= nfld; n++)
    {
      // ajout pour chaque face
      for (const std::pair<const ETK,E_Int>& elt : map2)
      {
        k = elt.first;
        ind = elt.second;
        p4 = k/(nvertex*nvertex*nvertex);
        p3 = k - p4*nvertex*nvertex*nvertex;
        p3 = p3/(nvertex*nvertex);
        p2 = k - p3*nvertex*nvertex - p4*nvertex*nvertex*nvertex;
        p2 = p2/nvertex;
        p1 = k - p2*nvertex - p3*nvertex*nvertex - p4*nvertex*nvertex*nvertex;
        (*fo)(ind+nvertex+nedges,n) = 0.25*( (*f)(p1,n)+ (*f)(p2,n) + (*f)(p3,n)+ (*f)(p4,n) );
      }
    }

    // Connectivity
    for (E_Int i = 0; i < nelts; i++)
    {
      n1 = (*cn)(i,1); n2 = (*cn)(i,2); n3 = (*cn)(i,3); n4 = (*cn)(i,4);
      n5 = (*cn)(i,5);
      
      EDGEINDEX(n1,n2,n6);
      EDGEINDEX(n2,n3,n7);
      EDGEINDEX(n3,n4,n8);
      EDGEINDEX(n4,n1,n9);
      EDGEINDEX(n1,n5,n10);
      EDGEINDEX(n2,n5,n11);
      EDGEINDEX(n3,n5,n12);
      EDGEINDEX(n4,n5,n13);
      FACEINDEX(n1,n2,n3,n4,n14);

      (*co)(i,1) = n1;
      (*co)(i,2) = n2;
      (*co)(i,3) = n3;
      (*co)(i,4) = n4;
      (*co)(i,5) = n5;
      (*co)(i,6) = n6+nvertex+1;
      (*co)(i,7) = n7+nvertex+1;
      (*co)(i,8) = n8+nvertex+1;
      (*co)(i,9) = n9+nvertex+1;
      (*co)(i,10) = n10+nvertex+1;
      (*co)(i,11) = n11+nvertex+1;
      (*co)(i,12) = n12+nvertex+1;
      (*co)(i,13) = n13+nvertex+1;
      (*co)(i,14) = n14+nvertex+nedges+1;
    }
    RELEASESHAREDU(o, fo, co); 
  }

  RELEASESHAREDB(res, array, f, cn);
  
  return o;
}

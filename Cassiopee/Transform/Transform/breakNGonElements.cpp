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

# include "transform.h"
# include "Connect/connect.h"
#include <unordered_map>

using namespace K_FLD;
using namespace std;

//=============================================================================
PyObject* K_TRANSFORM::breakElements(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, "O", &array))
  {
    return NULL;
  }
  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltTypes;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, ni, nj, nk, cnl, eltTypes);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "breakElements: array is invalid.");
    return NULL;
  }
  if (strcmp(eltTypes, "TRI")   == 0 || strcmp(eltTypes, "QUAD") == 0 ||
      strcmp(eltTypes, "TETRA") == 0 || strcmp(eltTypes, "HEXA") == 0 ||
      strcmp(eltTypes, "PENTA") == 0 || strcmp(eltTypes, "BAR")  == 0 || 
      strcmp(eltTypes, "PYRA")  == 0 || strcmp(eltTypes, "NODE") == 0)
  { RELEASESHAREDU(array, f, cnl); return array; }
  if (strcmp(eltTypes, "NGON") != 0 && strcmp(eltTypes, "MIXED") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "breakElements: elt type must be NGON or MIXED.");
    RELEASESHAREDU(array, f, cnl); return NULL;    
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  
  vector<E_Int> eltTypev; vector<FldArrayI*> cEV; vector<FldArrayF*> fields;
  if (strcmp(eltTypes, "NGON") == 0)
    breakNGonElements(*f, *cnl, cEV, fields, eltTypev, varString);
  else breakMixedElements(*f, *cnl, cEV, fields, eltTypev);

  PyObject* tpl;
  PyObject* l = PyList_New(0);
  char eltType[10]; strcpy(eltType, "BAR");

  for (size_t v = 0; v < cEV.size(); v++)
  {
    if (fields[v]->getSize() != 0) 
    {
      if (eltTypev[v] == 1) strcpy(eltType, "BAR");
      else if (eltTypev[v] == 2) strcpy(eltType, "TRI");
      else if (eltTypev[v] == 3) strcpy(eltType, "QUAD");
      else if (eltTypev[v] == 4) strcpy(eltType, "TETRA");
      else if (eltTypev[v] == 7) strcpy(eltType, "HEXA");
      else if (eltTypev[v] == 6) strcpy(eltType, "PENTA");
      else if (eltTypev[v] == 5) strcpy(eltType, "PYRA");
      else if (eltTypev[v] == 8) strcpy(eltType, "NGON");

      if (posx != 0 && posy != 0 && posz != 0)
        K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, 
                                     *fields[v], *cEV[v]);   
      tpl = K_ARRAY::buildArray3(*fields[v], varString, *cEV[v], eltType);
      PyList_Append(l, tpl); Py_DECREF(tpl);
    }
    delete fields[v]; delete cEV[v];
  }
  fields.clear(); cEV.clear(); eltTypev.clear();
  RELEASESHAREDU(array, f, cnl);
  return l;
}

//=============================================================================
void K_TRANSFORM::breakNGonElements(
  FldArrayF& field, FldArrayI& cNG, vector<FldArrayI*>& cEV,
  vector<FldArrayF*>& fields, vector<E_Int>& eltType, char* varString)
{ 
  E_Int nfld = field.getNfld();
  E_Int api = field.getApi(); if (api == 2) api = 3;
  E_Int shift = 1; if (api == 3) shift = 0;
  
  E_Int* ngon = cNG.getNGon(); E_Int* nface = cNG.getNFace();
  E_Int* indPG = cNG.getIndPG(); E_Int* indPH = cNG.getIndPH();
  E_Int ncells = cNG.getNElts();

  FldArrayI dimElt(ncells);
  K_CONNECT::getDimElts(cNG, dimElt);
  vector<vector<E_Int> > cEVNGon(ncells);
  K_CONNECT::connectNG2EV(cNG, cEVNGon);

  E_Int dim, nfacesl, sizeFN2 = 0, sizeEF2 = 0, nfacesngon = 0;
  vector<E_Int> verticesf; // sommets candidats image de la face

  E_Int nptsbar = 0, nptstri = 0, nptsquad = 0, nptstetra = 0, nptshexa = 0,
    nptspenta = 0, nptspyra = 0, nptsngon = 0;
  std::unordered_map<E_Int, E_Int> vMapBar, vMapTri, vMapQuad, vMapTetra,
    vMapPenta, vMapPyra, vMapHexa, vMapNGon, fMapNGon;
  vector<E_Int> etListBar, etListTri, etListQuad, etListTetra, etListPenta,
    etListPyra, etListHexa, etListNGon;

  for (E_Int et = 0; et < ncells; et++)
  {
    dim = dimElt[et];
    vector<E_Int>& vertices = cEVNGon[et]; // sommets associes a l'elt
    E_Int* elem = cNG.getElt(et, nfacesl, nface, indPH);

    if (nfacesl == 2 && dim == 1) // BAR
    {
      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        auto res = vMapBar.insert(std::make_pair(vertices[nov]-1, nptsbar));
        if (res.first->second == nptsbar) nptsbar++; // first time this vertex is encountered
      }
      etListBar.push_back(et);
    }
    else if (nfacesl == 3 && dim == 2) // TRI
    {
      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        auto res = vMapTri.insert(std::make_pair(vertices[nov]-1, nptstri));
        if (res.first->second == nptstri) nptstri++;
      }
      etListTri.push_back(et);
    }
    else if (nfacesl == 4 && dim == 2) // QUAD
    {
      // verification de la coherence de la numerotation des indices
      E_Int nv, vert0 = -1, fidx = elem[0];
      E_Int* face = cNG.getFace(fidx-1, nv, ngon, indPG);
      E_Int vert1 = face[0], vert2 = face[1];

      verticesf.clear();
      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        vert0 = vertices[nov];
        if (vert0 != vert1 && vert0 != vert2) verticesf.push_back(vert0);
      }
      E_Int vert4 = K_CONNECT::image(vert1, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);
      E_Int vert3 = K_CONNECT::image(vert2, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);

      if (vert3 == -1 || vert4 == -1) goto ngonLabel;
      else
      {
        vertices = {vert1, vert2, vert3, vert4};
        for (size_t nov = 0; nov < vertices.size(); nov++)
        {
          auto res = vMapQuad.insert(std::make_pair(vertices[nov]-1, nptsquad));
          if (res.first->second == nptsquad) nptsquad++;
        }
        etListQuad.push_back(et);
      }
    }
    else if (nfacesl == 4 && dim == 3) // TETRA
    {
      // recherche de la premiere face tri de l elt
      E_Int dummy, vert0 = -1, vert4 = -1;
      E_Int* face = cNG.getFace(elem[0]-1, dummy, ngon, indPG);
      E_Int vert1 = face[0], vert2 = face[1], vert3 = face[2];
      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        vert0 = vertices[nov];
        if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3) 
        {
          vert4 = vert0; break;
        }
      }
      vertices = {vert1, vert2, vert3, vert4};
      
      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        auto res = vMapTetra.insert(std::make_pair(vertices[nov]-1, nptstetra));
        if (res.first->second == nptstetra) nptstetra++;
      }
      etListTetra.push_back(et);
    }
    else if (nfacesl == 5 && dim == 3) // PENTA / PYRA
    {
      E_Int nv, nbnodes = 0;
      E_Int vert0 = -1, vert1 = -1, vert2 = -1, vert3 = -1, vert4 = -1, vert5 = -1;
      for (E_Int nf = 0; nf < nfacesl; nf++)
      {
        cNG.getFace(elem[nf]-1, nv, ngon, indPG);
        nbnodes += nv;
      }
      if (nbnodes == 16) // PYRA
      {
        // verification de la coherence de la numerotation des indices
        for (E_Int nf = 0; nf < nfacesl; nf++)
        {
          E_Int* face = cNG.getFace(elem[nf]-1, nv, ngon, indPG);
          if (nv == 4) // face = base quad
          {       
            vert1 = face[0]; vert2 = face[1]; vert3 = face[2]; vert4 = face[3];
            for (size_t nov = 0; nov < vertices.size(); nov++)
            {
              vert0 = vertices[nov];
              if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3 && vert0 != vert4)
              {
                vert5 = vert0; break;
              }
            }
          }
          if (vert5 != -1) break;
        }

        if (vert5 == -1) { etListNGon.push_back(et); } // TODO NGon
        else
        {
          vertices = {vert1, vert2, vert3, vert4, vert5};
          for (size_t nov = 0; nov < vertices.size(); nov++)
          {
            auto res = vMapPyra.insert(std::make_pair(vertices[nov]-1, nptspyra));
            if (res.first->second == nptspyra) nptspyra++;
          }
          etListPyra.push_back(et);
        }
      }
      else // PENTA
      {
        // verification de la coherence de la numerotation des indices
        E_Int fidx, vert6 = -1;
        for (E_Int nf = 0; nf < nfacesl; nf++)
        {
          fidx = elem[nf];
          E_Int* face = cNG.getFace(fidx-1, nv, ngon, indPG);
          if (nv == 3) // face = base tri
          {
            vert1 = face[0]; vert2 = face[1]; vert3 = face[2];
            verticesf.clear();
            for (size_t nov = 0; nov < vertices.size(); nov++)
            {
              vert0 = vertices[nov];
              if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3)
                verticesf.push_back(vert0);
            }

            vert4 = K_CONNECT::image(vert1, fidx, et, verticesf, cNG,
                                    ngon, nface, indPG, indPH);
            vert5 = K_CONNECT::image(vert2, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);
            vert6 = K_CONNECT::image(vert3, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);      
            break;        
          }
        }

        if ((vert4 != -1) && (vert5 != -1) && (vert6 != -1))
        {
          vertices = {vert1, vert2, vert3, vert4, vert5, vert6};
          for (size_t nov = 0; nov < vertices.size(); nov++)
          {
            auto res = vMapPenta.insert(std::make_pair(vertices[nov]-1, nptspenta));
            if (res.first->second == nptspenta) nptspenta++;
          }
          etListPenta.push_back(et);
        }
        else goto ngonLabel;
      }
    }
    else if (nfacesl == 6 && dim == 3) // HEXA
    {
      // verification de la coherence de la numerotation des indices
      E_Int nv, fidx = elem[0];
      E_Int* face = cNG.getFace(fidx-1, nv, ngon, indPG);
      E_Int vert1 = face[0], vert2 = face[1], vert3 = face[2], vert4 = face[3];

      verticesf.clear(); 
      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        E_Int vert0 = vertices[nov];
        if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3 && vert0 != vert4)
          verticesf.push_back(vert0);
      }
      E_Int vert5 = K_CONNECT::image(vert1, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);
      E_Int vert6 = K_CONNECT::image(vert2, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);
      E_Int vert7 = K_CONNECT::image(vert3, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);
      E_Int vert8 = K_CONNECT::image(vert4, fidx, et, verticesf, cNG,
                                     ngon, nface, indPG, indPH);

      if ((vert5 == -1) || (vert6 == -1) || (vert7 == -1) || (vert8 == -1)) goto ngonLabel;
      else
      {
        vertices = {vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8};
        for (size_t nov = 0; nov < vertices.size(); nov++)
        {
          auto res = vMapHexa.insert(std::make_pair(vertices[nov]-1, nptshexa));
          if (res.first->second == nptshexa) nptshexa++;
        }
        etListHexa.push_back(et);
      }
    }
    else
    {
      ngonLabel:;
      E_Int nv, fidx;
      for (E_Int nf = 0; nf < nfacesl; nf++)
      {
        fidx = elem[nf]-1;
        E_Int* face = cNG.getFace(fidx, nv, ngon, indPG);
        auto resF = fMapNGon.insert(std::make_pair(fidx, nfacesngon));
        if (resF.first->second == nfacesngon) // first time this face is encountered
        {
          for (E_Int nov = 0; nov < nv; nov++)
          {
            auto resV = vMapNGon.insert(std::make_pair(face[nov]-1, nptsngon));
            if (resV.first->second == nptsngon) nptsngon++; // first time this vertex is encountered
          }
          sizeFN2 += nv+shift;
          nfacesngon++;
        }
      }
      sizeEF2 += nfacesl+shift;
      etListNGon.push_back(et);
    }
  }

  E_Int netbar = etListBar.size();
  E_Int nettri = etListTri.size();
  E_Int netquad = etListQuad.size();
  E_Int nettetra = etListTetra.size();
  E_Int netpenta = etListPenta.size();
  E_Int netpyra = etListPyra.size();
  E_Int nethexa = etListHexa.size();
  E_Int netngon = etListNGon.size();

  FldArrayI* cEVbarp = new FldArrayI(netbar,2);
  FldArrayF* fbarp = new FldArrayF(nptsbar,nfld);
  FldArrayF& fbar = *fbarp; FldArrayI& cEVbar = *cEVbarp;

  FldArrayI* cEVtrip = new FldArrayI(nettri,3);
  FldArrayF* ftrip = new FldArrayF(nptstri,nfld);
  FldArrayF& ftri = *ftrip; FldArrayI& cEVtri = *cEVtrip;

  FldArrayI* cEVquadp = new FldArrayI(netquad,4);
  FldArrayF* fquadp = new FldArrayF(nptsquad,nfld);
  FldArrayF& fquad = *fquadp; FldArrayI& cEVquad = *cEVquadp;

  FldArrayI* cEVtetrap = new FldArrayI(nettetra,4);
  FldArrayF* ftetrap = new FldArrayF(nptstetra,nfld);
  FldArrayF& ftetra = *ftetrap; FldArrayI& cEVtetra = *cEVtetrap;

  FldArrayI* cEVpyrap = new FldArrayI(netpyra,5); 
  FldArrayF* fpyrap = new FldArrayF(nptspyra,nfld);
  FldArrayF& fpyra = *fpyrap; FldArrayI& cEVpyra = *cEVpyrap;
  
  FldArrayI* cEVpentap = new FldArrayI(netpenta,6);
  FldArrayF* fpentap = new FldArrayF(nptspenta, nfld);
  FldArrayF& fpenta = *fpentap; FldArrayI& cEVpenta = *cEVpentap;

  FldArrayI* cEVhexap = new FldArrayI(nethexa,8); 
  FldArrayF* fhexap = new FldArrayF(nptshexa,nfld);
  FldArrayF& fhexa = *fhexap; FldArrayI& cEVhexa = *cEVhexap;

  FldArrayI* cn2; FldArrayF* f2;
  E_Int* ngon2 = NULL; E_Int* nface2 = NULL;
  E_Int *indPG2 = NULL; E_Int* indPH2 = NULL;
  if (netngon)
  {
    E_Int ngonType = 1; // CGNSv3 compact array1
    if (api == 2) ngonType = 2; // CGNSv3, array2
    else if (api == 3) ngonType = 3; // force CGNSv4, array3
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, nptsngon, netngon,
                                         nfacesngon, "NGON", sizeFN2, sizeEF2,
                                         ngonType, false, api);
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    ngon2 = cn2->getNGon();
    nface2 = cn2->getNFace();
    if (api == 2 || api == 3)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }
  }

  #pragma omp parallel
  {
    E_Int et, ind1, ind2;
    
    #pragma omp for
    for (E_Int i = 0; i < netbar; i++) // BAR
    {
      et = etListBar[i];
      const vector<E_Int>& vertices = cEVNGon[et];

      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        ind1 = vertices[nov]-1;
        ind2 = vMapBar[ind1];
        cEVbar(i,nov+1) = ind2+1;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fbar(ind2,eq) = field(ind1,eq);
        }
      }
    }

    #pragma omp for
    for (E_Int i = 0; i < nettri; i++) // TRI
    {
      et = etListTri[i];
      const vector<E_Int>& vertices = cEVNGon[et];

      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        ind1 = vertices[nov]-1;
        ind2 = vMapTri[ind1];
        cEVtri(i,nov+1) = ind2+1;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          ftri(ind2,eq) = field(ind1,eq);
        }
      }
    }

    #pragma omp for
    for (E_Int i = 0; i < netquad; i++) // QUAD
    {
      et = etListQuad[i];
      const vector<E_Int>& vertices = cEVNGon[et];

      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        ind1 = vertices[nov]-1;
        ind2 = vMapQuad[ind1];
        cEVquad(i,nov+1) = ind2+1;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fquad(ind2,eq) = field(ind1,eq);
        }
      }
    }

    #pragma omp for
    for (E_Int i = 0; i < nettetra; i++) // TETRA
    {
      et = etListTetra[i];
      const vector<E_Int>& vertices = cEVNGon[et];

      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        ind1 = vertices[nov]-1;
        ind2 = vMapTetra[ind1];
        cEVtetra(i,nov+1) = ind2+1;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          ftetra(ind2,eq) = field(ind1,eq);
        }
      }
    }
    
    #pragma omp for
    for (E_Int i = 0; i < netpyra; i++) // PYRA
    {
      et = etListPyra[i];
      const vector<E_Int>& vertices = cEVNGon[et];

      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        ind1 = vertices[nov]-1;
        ind2 = vMapPyra[ind1];
        cEVpyra(i,nov+1) = ind2+1;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fpyra(ind2,eq) = field(ind1,eq);
        }
      }
    }
    
    #pragma omp for
    for (E_Int i = 0; i < netpenta; i++) // PENTA
    {
      et = etListPenta[i];
      const vector<E_Int>& vertices = cEVNGon[et];

      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        ind1 = vertices[nov]-1;
        ind2 = vMapPenta[ind1];
        cEVpenta(i,nov+1) = ind2+1;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fpenta(ind2,eq) = field(ind1,eq);
        }
      }
    }

    #pragma omp for
    for (E_Int i = 0; i < nethexa; i++) // HEXA
    {
      et = etListHexa[i];
      const vector<E_Int>& vertices = cEVNGon[et];

      for (size_t nov = 0; nov < vertices.size(); nov++)
      {
        ind1 = vertices[nov]-1;
        ind2 = vMapHexa[ind1];
        cEVhexa(i,nov+1) = ind2+1;
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          fhexa(ind2,eq) = field(ind1,eq);
        }
      }
    }
  }

  E_Int nv, fidx, et, ind1, ind2, ind3 = 0, ind4 = 0, ind5 = 1;
  if ((api == 2 || api == 3) && netngon) { indPG2[0] = 0; indPH2[0] = 0; }

  for (E_Int i = 0; i < netngon; i++) // NGON
  {
    et = etListNGon[i];
    E_Int* elem = cNG.getElt(et, nfacesl, nface, indPH);
    nface2[ind4] = nfacesl;
    if ((api == 2 || api == 3) && (i+1 < netngon)) indPH2[i+1] = nfacesl;

    for (E_Int nf = 0; nf < nfacesl; nf++)
    {
      fidx = elem[nf]-1;
      E_Int* face = cNG.getFace(fidx, nv, ngon, indPG);
      ngon2[ind3] = nv;
      if ((api == 2 || api == 3) && (i+1 < nfacesngon))
      {
        indPG2[ind5] = nv; ind5++;
      }
      for (E_Int nov = 0; nov < nv; nov++)
      {
        ind1 = face[nov]-1;
        ind2 = vMapNGon[ind1];
        for (E_Int eq = 1; eq <= nfld; eq++)
        {
          E_Float* f2p = f2->begin(eq);
          f2p[ind2] = field(ind1,eq);
        }
        ngon2[ind3+nov+shift] = ind2+1;
      }
      ind3 += nv+shift;
      
      nface2[ind4+nf+shift] = fMapNGon[fidx]+1;
    }
    ind4 += nfacesl+shift;
  }

  // Append non-void connectivities & fields
  if (netbar) {cEV.push_back(cEVbarp); fields.push_back(fbarp); eltType.push_back(1);}
  else { delete cEVbarp; delete fbarp; }
  if (nettri) {cEV.push_back(cEVtrip); fields.push_back(ftrip); eltType.push_back(2);}
  else { delete cEVtrip; delete ftrip; }
  if (netquad) {cEV.push_back(cEVquadp); fields.push_back(fquadp); eltType.push_back(3);}
  else { delete cEVquadp; delete fquadp; }
  if (nettetra) {cEV.push_back(cEVtetrap); fields.push_back(ftetrap); eltType.push_back(4);}
  else { delete cEVtetrap; delete ftetrap; }
  if (netpyra) {cEV.push_back(cEVpyrap); fields.push_back(fpyrap); eltType.push_back(5);}
  else { delete cEVpyrap; delete fpyrap; }
  if (netpenta) {cEV.push_back(cEVpentap); fields.push_back(fpentap); eltType.push_back(6);}
  else { delete cEVpentap; delete fpentap; }
  if (nethexa) {cEV.push_back(cEVhexap); fields.push_back(fhexap); eltType.push_back(7);}
  else { delete cEVhexap; delete fhexap; }
  if (netngon) {cEV.push_back(cn2); fields.push_back(f2); eltType.push_back(8);}
}

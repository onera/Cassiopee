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

# include "transform.h"
# include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
PyObject* K_TRANSFORM::breakElements(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array))
  {
      return NULL;
  }
  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varStrings; char* eltTypes;
  res = K_ARRAY::getFromArray(array, varStrings, 
                              f, ni, nj, nk, cnl, eltTypes, true);

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
  if (strcmp(eltTypes, "NGON") != 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "breakElements: elt type must be NGON.");
    RELEASESHAREDU(array, f, cnl); return NULL;    
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varStrings); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varStrings); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varStrings); posz++;
  
  PyObject* tpl;
  PyObject* l = PyList_New(0);

  vector<E_Int> eltTypev; vector<FldArrayI*> cEV; vector<FldArrayF*> fields;
  breakNGonElements(*f, *cnl, cEV, fields, eltTypev);
  char eltType[10];
  strcpy(eltType, "BAR");
  E_Int cEVsize = cEV.size();

  for (int v = 0; v < cEVsize; v++)
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
      tpl = K_ARRAY::buildArray(*fields[v], varStrings, *cEV[v], eltTypev[v]);
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
  FldArrayF& field, FldArrayI& cNG, 
  vector<FldArrayI*>& cEV, vector<FldArrayF*>& fields, vector<E_Int>& eltType)
{ 
  E_Int npts = field.getSize(); E_Int nfld = field.getNfld();
  E_Int* cnp = cNG.begin();
  E_Int nfaces = cnp[0];
  E_Int sizeFN = cnp[1];
  E_Int ncells = cnp[sizeFN+2]; 
  E_Int sizeEF = cnp[sizeFN+3];
  //E_Int* cEF = cnp+4+sizeFN;// debut connectivite EF
  //E_Int* cFN = cnp+2;// debut connectivite FN
  FldArrayI posFace(nfaces);
  K_CONNECT::getPosFaces(cNG, posFace);
  FldArrayI posElt(ncells);
  K_CONNECT::getPosElts(cNG, posElt);
  FldArrayI dimElt(ncells);
  K_CONNECT::getDimElts(cNG, posFace, dimElt);

  vector< vector<E_Int> > cEVNGon(ncells);
  K_CONNECT::connectNG2EV(cNG, cEVNGon);

  E_Int netbar = 0; E_Int nptsbar = 0;
  E_Int nettri = 0; E_Int nptstri = 0;
  E_Int netquad = 0; E_Int nptsquad = 0;
  E_Int nettetra = 0; E_Int nptstetra = 0;
  E_Int nethexa = 0; E_Int nptshexa = 0;
  E_Int netpenta = 0; E_Int nptspenta = 0;
  E_Int netpyra = 0; E_Int nptspyra = 0;
  //E_Int netngon = 0; E_Int nptsngon = 0;

  FldArrayI* cEVbarp = new FldArrayI(ncells, 2); 
  FldArrayF* fbarp = new FldArrayF(npts, nfld);
  FldArrayF& fbar = *fbarp; FldArrayI& cEVbar = *cEVbarp;
  FldArrayI indirbF(npts); indirbF.setAllValuesAt(-1);
  E_Int* indirb = indirbF.begin();

  FldArrayI* cEVtrip = new FldArrayI(ncells,3); 
  FldArrayF* ftrip = new FldArrayF(npts,nfld);
  FldArrayF& ftri = *ftrip; FldArrayI& cEVtri = *cEVtrip;
  FldArrayI indirtF(npts); indirtF.setAllValuesAt(-1);
  E_Int* indirt = indirtF.begin();

  FldArrayI* cEVquadp = new FldArrayI(ncells,4); 
  FldArrayF* fquadp = new FldArrayF(npts,nfld);
  FldArrayF& fquad = *fquadp; FldArrayI& cEVquad = *cEVquadp;
  FldArrayI indirqF(npts); indirqF.setAllValuesAt(-1);
  E_Int* indirq = indirqF.begin();
  
  FldArrayI* cEVtetrap = new FldArrayI(ncells,4); 
  FldArrayF* ftetrap = new FldArrayF(npts,nfld);
  FldArrayF& ftetra = *ftetrap; FldArrayI& cEVtetra = *cEVtetrap;
  FldArrayI indirttF(npts); indirttF.setAllValuesAt(-1);
  E_Int* indirtt = indirttF.begin();

  FldArrayI* cEVpentap = new FldArrayI(ncells, 6); 
  FldArrayF* fpentap = new FldArrayF(npts, nfld);
  FldArrayF& fpenta = *fpentap; FldArrayI& cEVpenta = *cEVpentap;
  FldArrayI indirpF(npts); indirpF.setAllValuesAt(-1);
  E_Int* indirp = indirpF.begin();

  FldArrayI* cEVpyrap = new FldArrayI(ncells, 5); 
  FldArrayF* fpyrap = new FldArrayF(npts,nfld);
  FldArrayF& fpyra = *fpyrap; FldArrayI& cEVpyra = *cEVpyrap;
  FldArrayI indiryF(npts); indiryF.setAllValuesAt(-1);
  E_Int* indiry = indiryF.begin();

  FldArrayI* cEVhexap = new FldArrayI(ncells, 8); 
  FldArrayF* fhexap = new FldArrayF(npts,nfld);
  FldArrayF& fhexa = *fhexap; FldArrayI& cEVhexa = *cEVhexap;
  FldArrayI indirhF(npts); indirhF.setAllValuesAt(-1);
  E_Int* indirh = indirhF.begin();

  FldArrayF* fngonp = new FldArrayF(npts, nfld); FldArrayF& fngon = *fngonp;
  FldArrayI cEF2(sizeEF); FldArrayI cFN2(sizeFN);
  FldArrayI indirnF(npts); indirnF.setAllValuesAt(-1);
  E_Int* indirn = indirnF.begin();

  // nfaces tot ?
  E_Int* ptref = cNG.begin(); E_Int* ptrfn = cNG.begin();
  FldArrayI indirFF(nfaces); indirFF.setAllValuesAt(-1); // pour les faces NGON
  E_Int* indirF = indirFF.begin();
  E_Int vert1, vert2, vert3, vert4, vert5, vert6, vert7, vert8;
  E_Int dim;
  E_Int sizeEF2=0; E_Int sizeFN2=0; E_Int npts2=0; 
  E_Int nelts2=0; E_Int nfaces2=0;
  vector<E_Int> verticesf;// sommets candidats image de la face

  for (E_Int et = 0; et < ncells; et++)
  {
    dim = dimElt[et];
    ptref = cnp+posElt[et];
    E_Int nfacesl = ptref[0];//nb faces pour l element et
    vector<E_Int>& vertices = cEVNGon[et];// sommets associes a l'elt
    if (nfacesl == 2 && dim == 1) //BAR
    {
      for (E_Int nov = 0; nov < 2; nov++)
      {
        if (indirb[vertices[nov]-1] == -1)  
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            fbar(nptsbar,eq) = field(vertices[nov]-1,eq);    
          nptsbar++;
          indirb[vertices[nov]-1] = nptsbar;
        }
      }
      for (E_Int nov = 0; nov < 2; nov++)
      {  
        cEVbar(netbar, nov+1) = indirb[vertices[nov]-1];
      }
      netbar++;
    }
    else if ((nfacesl == 3)&&(dim == 2)) // TRI
    {
      for (E_Int nov = 0; nov < 3; nov++)
      {
        if (indirt[vertices[nov]-1] == -1)  
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            ftri(nptstri,eq) = field(vertices[nov]-1,eq);    
          nptstri++;
          indirt[vertices[nov]-1] = nptstri;
        }
      }
      for (E_Int nov = 0; nov < 3; nov++)
      {  
        cEVtri(nettri,nov+1) = indirt[vertices[nov]-1];
      }
      nettri++;     
    }
    else if (nfacesl == 4 && dim == 2) // QUAD
    {
      E_Int nbnodes = 0;
      for (E_Int nf = 1; nf <= nfacesl; nf++)
      {
        E_Int face = ptref[nf]-1;
        ptrfn = cnp+posFace[face];
        nbnodes += ptrfn[0];
      }
      vert1 = -1; vert2 = -1; vert3 = 0; vert4 = 0;
      // recherche des faces de l'elt
      for (E_Int nf = 1; nf <= nfacesl; nf++)
      {
        E_Int face = ptref[nf];
        ptrfn = cnp+posFace[face-1];// pointeur sur la face
        vert1 = ptrfn[1]; vert2 = ptrfn[2]; 
        for (unsigned int novert = 0; novert < vertices.size(); novert++)
        {
          E_Int vert0 = vertices[novert];
          if (vert0 != vert1 && vert0 != vert2) verticesf.push_back(vert0);
        }
        // verification de la coherence de la numerotation des indices
        vert4 = K_CONNECT::image(vert1, face, et, verticesf, posFace, posElt, cNG);
        vert3 = K_CONNECT::image(vert2, face, et, verticesf, posFace, posElt, cNG);
        if (vert3 == -1 || vert4 == -1) goto ngon;
        else goto quad0;
      }
      quad0:;
      //construction de l'elt quad
      if (indirq[vert1-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fquad(nptsquad,eq) = field(vert1-1,eq); 
        nptsquad++; indirq[vert1-1] = nptsquad;
      }
      if (indirq[vert2-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fquad(nptsquad,eq) = field(vert2-1,eq); 
        nptsquad++; indirq[vert2-1] = nptsquad;
      }
      if (indirq[vert3-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fquad(nptsquad,eq) = field(vert3-1,eq); 
        nptsquad++; indirq[vert3-1] = nptsquad;
      }
      if (indirq[vert4-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fquad(nptsquad,eq) = field(vert4-1,eq); 
        nptsquad++; indirq[vert4-1] = nptsquad;
      }
      cEVquad(netquad,1) = indirq[vert1-1];
      cEVquad(netquad,2) = indirq[vert2-1];
      cEVquad(netquad,3) = indirq[vert3-1];
      cEVquad(netquad,4) = indirq[vert4-1];
      netquad++;
    }//fin quad
    else if ((nfacesl == 4)&&(dim == 3)) // TETRA
    {
      E_Int nbnodes = 0;
      for (E_Int nf = 1; nf <= nfacesl; nf++)
      {
        E_Int face = ptref[nf]-1;
        ptrfn = cnp+posFace[face];
        nbnodes += ptrfn[0];
      }
      vert1 = -1; vert2 = -1; vert3 = -1; vert4 = -1; 
      //recherche de la premiere face tri de l elt
      E_Int face = ptref[1];
      ptrfn = cnp+posFace[face-1];// pointeur sur la face    
      vert1 = ptrfn[1]; vert2 = ptrfn[2]; vert3 = ptrfn[3];          
      for (unsigned int novert = 0; novert < vertices.size(); novert++)
      {
        E_Int vert0 = vertices[novert];
        if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3) 
        {vert4 = vert0; goto tetra0; }
      }
      tetra0:;
      //construction de l'elt tetra
      if (indirtt[vert1-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          ftetra(nptstetra,eq) = field(vert1-1,eq); 
        nptstetra++; indirtt[vert1-1] = nptstetra;
      }
      if (indirtt[vert2-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          ftetra(nptstetra,eq) = field(vert2-1,eq); 
        nptstetra++; indirtt[vert2-1] = nptstetra;
      }
      if (indirtt[vert3-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          ftetra(nptstetra,eq) = field(vert3-1,eq); 
        nptstetra++; indirtt[vert3-1] = nptstetra;
      }
      if (indirtt[vert4-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          ftetra(nptstetra,eq) = field(vert4-1,eq); 
        nptstetra++; indirtt[vert4-1] = nptstetra;
      }
      cEVtetra(nettetra,1) = indirtt[vert1-1];
      cEVtetra(nettetra,2) = indirtt[vert2-1];
      cEVtetra(nettetra,3) = indirtt[vert3-1];
      cEVtetra(nettetra,4) = indirtt[vert4-1];
      nettetra++;
    }//fin tetra
    else if ((nfacesl == 5)&&(dim == 3)) // PENTA/PYRA
    {
      E_Int nbnodes = 0;
      for (E_Int nf = 1; nf <= nfacesl; nf++)
      {
        E_Int face = ptref[nf]-1;
        ptrfn = cnp+posFace[face];
        nbnodes += ptrfn[0];
      }
      if (nbnodes == 16) // PYRA
      {
        vert1 = -1; vert2 = -1; vert3 = -1; vert4 = -1; vert5 = -1;
        //recherche de la face quad de l elt
        for (E_Int nf = 1; nf <= nfacesl; nf++)
        {
          E_Int face = ptref[nf];
          ptrfn = cnp+posFace[face-1];// pointeur sur la face
          if (ptrfn[0] == 4) // face = base quad
          {       
            vert1 = ptrfn[1]; vert2 = ptrfn[2]; vert3 = ptrfn[3]; vert4 = ptrfn[4];         
            for (unsigned int novert = 0; novert < vertices.size(); novert++)
            {
              E_Int vert0 = vertices[novert];
              if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3 && vert0 != vert4) 
              { vert5 = vert0; goto pyra0; }
            }
            goto ngon;
          }
        }
        pyra0:;
        //construction de l'elt pyra
        for (E_Int eq = 1; eq <= nfld; eq++) 
        {
          if (indiry[vert1-1] == -1) 
          {
            for (E_Int eq = 1; eq <= nfld; eq++) {fpyra(nptspyra,eq) = field(vert1-1,eq);}
            nptspyra++; indiry[vert1-1] = nptspyra;
          }
          if (indiry[vert2-1] == -1) 
          {
            for (E_Int eq = 1; eq <= nfld; eq++) {fpyra(nptspyra,eq) = field(vert2-1,eq);}
            nptspyra++; indiry[vert2-1] = nptspyra;
          }
          if (indiry[vert3-1] == -1) 
          {
            for (E_Int eq = 1; eq <= nfld; eq++) {fpyra(nptspyra,eq) = field(vert3-1,eq);}
            nptspyra++; indiry[vert3-1] = nptspyra;
          }
          if (indiry[vert4-1] == -1) 
          {
            for (E_Int eq = 1; eq <= nfld; eq++) {fpyra(nptspyra,eq) = field(vert4-1,eq);}
            nptspyra++; indiry[vert4-1] = nptspyra;
          }
          if (indiry[vert5-1] == -1) 
          {
            for (E_Int eq = 1; eq <= nfld; eq++) {fpyra(nptspyra,eq) = field(vert5-1,eq);}
            nptspyra++; indiry[vert5-1] = nptspyra;
          }
        }
        cEVpyra(netpyra,1) = indiry[vert1-1];
        cEVpyra(netpyra,2) = indiry[vert2-1];
        cEVpyra(netpyra,3) = indiry[vert3-1];
        cEVpyra(netpyra,4) = indiry[vert4-1];
        cEVpyra(netpyra,5) = indiry[vert5-1];
        netpyra++; 
      }//fin PYRA
      else // PENTA
      {
        vert1 = -1; vert2 = -1; vert3 = -1; vert4 = -1; vert5 = -1; vert6 = -1;
        // recherche des faces tri de l'elt
        for (E_Int nf = 1; nf <= nfacesl; nf++)
        {
          E_Int face = ptref[nf];
          ptrfn = cnp+posFace[face-1]; // pointeur sur la face
          if (ptrfn[0] == 3) // face = base tri
          {
            vert1 = ptrfn[1]; vert2 = ptrfn[2]; vert3 = ptrfn[3];
            verticesf.clear(); 
            for (unsigned int novert = 0; novert < vertices.size(); novert++)
            {
              E_Int vert0 = vertices[novert];
              if (vert0 != vert1 && vert0 != vert2 &&  vert0 != vert3) verticesf.push_back(vert0);
            }
            // verification de la coherence de la numerotation des indices
            vert4 = K_CONNECT::image(vert1, face, et, verticesf, posFace, posElt, cNG);
            vert5 = K_CONNECT::image(vert2, face, et, verticesf, posFace, posElt, cNG);
            vert6 = K_CONNECT::image(vert3, face, et, verticesf, posFace, posElt, cNG);      
            if ((vert4 == -1)||(vert5 == -1)||(vert6 == -1)) goto ngon;
            else goto penta0;            
          }
        }
        goto ngon;
        penta0:;
        //construction de l'elt penta
        if (indirp[vert1-1] == -1) 
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            fpenta(nptspenta,eq) = field(vert1-1,eq); 
          nptspenta++; indirp[vert1-1] = nptspenta;
        }
        if (indirp[vert2-1] == -1) 
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            fpenta(nptspenta,eq) = field(vert2-1,eq); 
          nptspenta++; indirp[vert2-1] = nptspenta;
        }
        if (indirp[vert3-1] == -1) 
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            fpenta(nptspenta,eq) = field(vert3-1,eq); 
          nptspenta++; indirp[vert3-1] = nptspenta;
        }
        if (indirp[vert4-1] == -1) 
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            fpenta(nptspenta,eq) = field(vert4-1,eq); 
          nptspenta++; indirp[vert4-1] = nptspenta;
        }
        if (indirp[vert5-1] == -1) 
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            fpenta(nptspenta,eq) = field(vert5-1,eq); 
          nptspenta++; indirp[vert5-1] = nptspenta;
        }
        if (indirp[vert6-1] == -1) 
        {
          for (E_Int eq = 1; eq <= nfld; eq++) 
            fpenta(nptspenta,eq) = field(vert6-1,eq); 
          nptspenta++; indirp[vert6-1] = nptspenta;
        }
        
        cEVpenta(netpenta,1) = indirp[vert1-1];
        cEVpenta(netpenta,2) = indirp[vert2-1];
        cEVpenta(netpenta,3) = indirp[vert3-1];
        cEVpenta(netpenta,4) = indirp[vert4-1];
        cEVpenta(netpenta,5) = indirp[vert5-1];
        cEVpenta(netpenta,6) = indirp[vert6-1];
        netpenta++; 

      }//fin PENTA
    }// FIN PYRA/PENTA
    else if ((nfacesl == 6)&&(dim == 3)) // HEXA
    {
      vert1 = -1; vert2 = -1; vert3 = -1; vert4 = -1; 
      vert5 = -1; vert6 = -1; vert7 = -1; vert8 = -1;
      //recherche de la face quad de l elt
      for (E_Int nf = 1; nf <= nfacesl; nf++)
      {
        E_Int face = ptref[nf];
        ptrfn = cnp+posFace[face-1];// pointeur sur la face
        vert1 = ptrfn[1]; vert2 = ptrfn[2]; vert3 = ptrfn[3]; vert4 = ptrfn[4];        
        verticesf.clear(); 
        for (unsigned int novert = 0; novert < vertices.size(); novert++)
        {
          E_Int vert0 = vertices[novert];
          if (vert0 != vert1 && vert0 != vert2 && vert0 != vert3 && vert0 != vert4) verticesf.push_back(vert0);
        }
        // verification de la coherence de la numerotation des indices
        vert5 = K_CONNECT::image(vert1, face, et, verticesf, posFace, posElt, cNG);
        vert6 = K_CONNECT::image(vert2, face, et, verticesf, posFace, posElt, cNG);
        vert7 = K_CONNECT::image(vert3, face, et, verticesf, posFace, posElt, cNG);   
        vert8 = K_CONNECT::image(vert4, face, et, verticesf, posFace, posElt, cNG);   
        if ((vert5 == -1)||(vert6 == -1)||(vert7 == -1)||(vert8 == -1)) goto ngon;
        else goto hexa0;
      }
      hexa0:;
      //construction de l'elt Hexa
      if (indirh[vert1-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert1-1,eq); 
        nptshexa++; indirh[vert1-1] = nptshexa;
      }
      if (indirh[vert2-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert2-1,eq); 
        nptshexa++; indirh[vert2-1] = nptshexa;
      }
      if (indirh[vert3-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert3-1,eq); 
        nptshexa++; indirh[vert3-1] = nptshexa;
      }
      if (indirh[vert4-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert4-1,eq); 
        nptshexa++; indirh[vert4-1] = nptshexa;
      }
      if (indirh[vert5-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert5-1,eq); 
        nptshexa++; indirh[vert5-1] = nptshexa;
      }
      if (indirh[vert6-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert6-1,eq); 
        nptshexa++; indirh[vert6-1] = nptshexa;
      }
      if (indirh[vert7-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert7-1,eq); 
        nptshexa++; indirh[vert7-1] = nptshexa;
      }
      if (indirh[vert8-1] == -1) 
      {
        for (E_Int eq = 1; eq <= nfld; eq++) 
          fhexa(nptshexa,eq) = field(vert8-1,eq); 
        nptshexa++; indirh[vert8-1] = nptshexa;
      }
      cEVhexa(nethexa,1) = indirh[vert1-1];
      cEVhexa(nethexa,2) = indirh[vert2-1];
      cEVhexa(nethexa,3) = indirh[vert3-1];
      cEVhexa(nethexa,4) = indirh[vert4-1];
      cEVhexa(nethexa,5) = indirh[vert5-1];
      cEVhexa(nethexa,6) = indirh[vert6-1];
      cEVhexa(nethexa,7) = indirh[vert7-1];
      cEVhexa(nethexa,8) = indirh[vert8-1];
      nethexa++; 
    }
    else // NGON
    {
      ngon:;
      //connectivite elt/faces
      cEF2[sizeEF2] = nfacesl; 
      for (E_Int nf = 1; nf <= nfacesl; nf++)
      {
        E_Int face = ptref[nf]-1;
        if (indirF[face] == -1) 
        {
          ptrfn = cnp+posFace[face];
          cFN2[sizeFN2] = ptrfn[0];          
          for (E_Int nv = 1; nv <= ptrfn[0]; nv++)//parcours des noeuds
          {
            E_Int indvert = ptrfn[nv]-1;
            if ( indirn[indvert] == -1 ) 
            {
              for (E_Int eq = 1; eq <= nfld; eq++)
                fngon(npts2,eq) = field(indvert,eq);
              npts2++;
              indirn[indvert] = npts2; 
            }
            cFN2[sizeFN2+nv] = indirn[indvert];
          }
          sizeFN2+= ptrfn[0]+1;
          nfaces2++; 
          indirF[face] = nfaces2;//demarre a 1
        }
        cEF2[sizeEF2+nf] = indirF[face];
      }// fin boucle sur les faces 
      sizeEF2 += nfacesl+1;
      nelts2++;
    }//NGON
  }
  // BAR
  fbarp->reAllocMat(nptsbar,nfld); cEVbarp->reAllocMat(netbar,2);
  cEV.push_back(cEVbarp); fields.push_back(fbarp); eltType.push_back(1);
  //TRI
  ftrip->reAllocMat(nptstri,nfld); cEVtrip->reAllocMat(nettri,3);
  cEV.push_back(cEVtrip); fields.push_back(ftrip); eltType.push_back(2);
  //QUAD
  fquadp->reAllocMat(nptsquad,nfld); cEVquadp->reAllocMat(netquad,4);
  cEV.push_back(cEVquadp); fields.push_back(fquadp); eltType.push_back(3);
  //TETRA
  ftetrap->reAllocMat(nptstetra,nfld); cEVtetrap->reAllocMat(nettetra,4);
  cEV.push_back(cEVtetrap); fields.push_back(ftetrap); eltType.push_back(4);
  //PYRA
  fpyrap->reAllocMat(nptspyra,nfld); cEVpyrap->reAllocMat(netpyra,5);
  cEV.push_back(cEVpyrap); fields.push_back(fpyrap); eltType.push_back(5);
  //PENTA
  fpentap->reAllocMat(nptspenta,nfld); cEVpentap->reAllocMat(netpenta,6);
  cEV.push_back(cEVpentap); fields.push_back(fpentap); eltType.push_back(6);
  //HEXA
  fhexap->reAllocMat(nptshexa,nfld); cEVhexap->reAllocMat(nethexa,8);
  cEV.push_back(cEVhexap); fields.push_back(fhexap); eltType.push_back(7);
  // NGON
  fngonp->reAllocMat(npts2,nfld); cFN2.resize(sizeFN2); cEF2.resize(sizeEF2);
  FldArrayI* cnNGON = new FldArrayI(4+sizeEF2+sizeFN2);
  E_Int* cnNGONp = cnNGON->begin();
  cnNGONp[0] = nfaces2; cnNGONp[1] = sizeFN2;
  for (E_Int i = 0; i < sizeFN2; i++) cnNGONp[2+i] = cFN2[i];
  cnNGONp[sizeFN2+2] = nelts2;
  cnNGONp[sizeFN2+3] = sizeEF2;
  for (E_Int i = 0; i < sizeEF2; i++) cnNGONp[4+sizeFN2+i] = cEF2[i];
  cEV.push_back(cnNGON); fields.push_back(fngonp); eltType.push_back(8);
}

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
# include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Rm ghost cells to a NGON array */
//=============================================================================
PyObject* K_CONVERTER::rmGhostCellsNGonNodes(PyObject* self, PyObject* args)
{
  PyObject *arrayN;
  E_Int depth;
  if (!PYPARSETUPLE_(args, O_ I_, &arrayN, &depth)) return NULL;

  if (depth < 1) 
  {
    PyErr_SetString(PyExc_ValueError, 
                    "rmGhostCellsNGon: depth must be > 0.");
    return NULL;
  }
  E_Int ni, nj, nk;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(arrayN, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(arrayN, f);
    PyErr_SetString(PyExc_TypeError, 
                    "rmGhostCells: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType,"NGON") != 0 && strcmp(eltType,"NGON*") != 0 ) 
  {
    RELEASESHAREDU(arrayN,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "rmGhostCells: array must be a NGON.");
    return NULL;
  }

  E_Int api = f->getApi();
  FldArrayI cFE; K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int* ptr = cn->begin();
  E_Int sizeFN = ptr[1]; ptr += sizeFN+2;
  E_Int ne = ptr[0];
  E_Int sizeEF = ptr[1]; ptr += 2;
  
  E_Int nfe, e1, e2, face;
  FldArrayI co(sizeEF); // elements interieurs
  E_Int nint=0; // nbre d'elts interieurs
  E_Int sizeco = 0; // taille de la connectivite d'elements interieurs
  E_Int* ptro = co.begin();
  E_Int sizengbrs, indv;
  //-----
  // Compute tag=1 : elts to be removed
  // 1st layer: external face  - Upper layers: neighbours of an elt tagged 1
  FldArrayI tag(ne); tag.setAllValuesAtNull();
  for (E_Int i = 0; i < ne; i++)
  {
    nfe = ptr[0];
    for (E_Int n = 1; n <= nfe; n++)
    {
      face = ptr[n]-1; 
      e1 = cFE1[face]; e2 = cFE2[face];
      if (e1 == 0 && e2 > 0 )  tag[i] = 1;
      else if (e1 > 0 && e2 == 0) tag[i] = 1;
    }
    ptr += nfe+1;
  }


  if (depth > 1)
  {
    vector < vector<E_Int> > cEEN(ne);
    K_CONNECT::connectFE2EENbrs(cFE, cEEN);  
    for (E_Int d = 2; d <= depth; d++)
    {
      ptr = cn->begin()+sizeFN+4;
      for (E_Int i = 0; i < ne; i++)
      {
        nfe = ptr[0];
        vector<E_Int>& ngbrs = cEEN[i];
        sizengbrs = ngbrs.size();
        if (tag[i] == d-1)
        {
          for (E_Int ng = 0; ng < sizengbrs; ng++)
          {
            indv = ngbrs[ng];
            if ( tag[indv] == 0) tag[indv] = d;
          }
        }
        ptr += nfe+1;
      }
    }
  }
 
  //-----
  // Select interior elements
  ptr = cn->begin()+sizeFN+4;
  for (E_Int i = 0; i < ne; i++)
  {
    nfe = ptr[0];
    
    if (tag[i] == 0) // recopie l'element 
    {
      ptro[0] = nfe; 
      
      for (E_Int n = 1; n <= nfe; n++) 
        ptro[n] = ptr[n];  
      
      ptro += nfe+1; sizeco += nfe+1; nint++;
    }
    ptr += nfe+1;
  }
  co.reAlloc(sizeco); ptro = co.begin();
  // Cree la nouvelle connectivite complete
  FldArrayI* cout = new FldArrayI(sizeFN+4+sizeco);
  E_Int* coutp = cout->begin();
  ptr = cn->begin();
  for (E_Int i = 0; i < sizeFN+2; i++) coutp[i] = ptr[i];
  coutp += sizeFN+2;
  coutp[0] = nint;
  coutp[1] = sizeco; coutp += 2;
  for (E_Int i = 0; i < sizeco; i++) coutp[i] = ptro[i];

  // Sortie  
  K_CONNECT::cleanUnreferencedFacesNGon(*cout);

  PyObject* tpl = K_ARRAY::buildArray3(*f, varString, *cout, eltType, api);
  RELEASESHAREDU(arrayN, f, cn); 
  delete cout;
  return tpl;
}  

//=============================================================================
/* Rm ghost cells to a NGON array */
//=============================================================================
PyObject* K_CONVERTER::rmGhostCellsNGonCenters(PyObject* self, PyObject* args)
{
  PyObject *arrayC;
  E_Int depth;
  if (!PYPARSETUPLE_(args, O_ I_, &arrayC, &depth)) return NULL;

  if (depth < 1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "rmGhostCells: depth must be > 0.");
    return NULL;
  }
 
  E_Int nic, njc, nkc;
  char* varStringC; char* eltTypec;
  FldArrayF* fc; FldArrayI* cn;
  E_Int resc = K_ARRAY::getFromArray3(arrayC, varStringC, fc, nic, njc, nkc, 
                                      cn, eltTypec);
  if (resc != 2)
  {
    if ( resc == 1) RELEASESHAREDS(arrayC,fc);
    PyErr_SetString(PyExc_TypeError, 
                    "rmGhostCells: array is invalid.");
    return NULL;
  }
  if (strcmp(eltTypec,"NGON*") != 0 ) 
  {
    
    PyErr_SetString(PyExc_TypeError,
                    "rmGhostCells: array must be NGON*.");
    RELEASESHAREDU(arrayC,fc,cn); return NULL;
  }

  E_Int api = fc->getApi();
  FldArrayI cFE; K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int* ptr = cn->begin();
  E_Int sizeFN = ptr[1]; ptr += sizeFN+2;
  E_Int ne = ptr[0];
  E_Int sizeEF = ptr[1]; ptr += 2;
  
  E_Int nfe, e1, e2, face;
  FldArrayI co(sizeEF); // elements interieurs
  E_Int nint=0; // nbre d'elts interieurs
  E_Int sizeco = 0; // taille de la connectivite d'elements interieurs
  E_Int* ptro = co.begin();
  E_Int sizengbrs, indv;
  //-----
  // Compute tag=1 : elts to be removed
  // 1st layer: external face  - Upper layers: neighbours of an elt tagged 1
  FldArrayI tag(ne); tag.setAllValuesAtNull();
  for (E_Int i = 0; i < ne; i++)
  {
    nfe = ptr[0];
    for (E_Int n = 1; n <= nfe; n++)
    {
      face = ptr[n]-1; 
      e1 = cFE1[face]; e2 = cFE2[face];      
      if (e1 == 0 && e2 > 0 )  tag[i] = 1;
      else if (e1 > 0 && e2 == 0) tag[i] = 1;
    }
    ptr += nfe+1;
  }

  if ( depth > 1)
  {
    vector < vector<E_Int> > cEEN(ne);
    K_CONNECT::connectFE2EENbrs(cFE, cEEN);  
    for (E_Int d = 2; d <= depth; d++)
    {
      ptr = cn->begin()+sizeFN+4;
      for (E_Int i = 0; i < ne; i++)
      {
        nfe = ptr[0];
        vector<E_Int>& ngbrs = cEEN[i];
        sizengbrs = ngbrs.size();
        if (tag[i] == d-1)
        {
          for (E_Int ng = 0; ng < sizengbrs; ng++)
          {
            indv = ngbrs[ng];
            if ( tag[indv] == 0) tag[indv] = d;
          }
        }
        ptr += nfe+1;
      }
    }
  }
 
  //-----
  // Select interior elements
  ptr = cn->begin()+sizeFN+4;
  for (E_Int i = 0; i < ne; i++)
  {
    nfe = ptr[0];
    if (tag[i] == 0) // recopie l'element 
    {
      ptro[0] = nfe; 
      
      for (E_Int n = 1; n <= nfe; n++) 
        ptro[n] = ptr[n];  
      
      ptro += nfe+1; sizeco += nfe+1; nint++;
    }
    ptr += nfe+1;
  }
  co.reAlloc(sizeco); ptro = co.begin();

  // Les points ajoutes sont des points crees a la fin par addGhostCells
  // si on a conserve l ordre de addGhostCells, pas besoin de remodifier la connectivite Faces/Noeuds
  E_Int nfld = fc->getNfld();
  FldArrayF* fcout = new FldArrayF(nint,nfld);
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fpc = fcout->begin(eq);
    E_Float* fp = fc->begin(eq);
    for (E_Int i = 0; i < nint; i++)
    {
      fpc[i] = fp[i];
    }
  }
  // Cree la nouvelle connectivite complete
  FldArrayI* cout = new FldArrayI(sizeFN+4+sizeco);
  E_Int* coutp = cout->begin();
  ptr = cn->begin();
  for (E_Int i = 0; i < sizeFN+2; i++) coutp[i] = ptr[i];
  coutp += sizeFN+2;
  coutp[0] = nint;
  coutp[1] = sizeco; coutp += 2;
  for (E_Int i = 0; i < sizeco; i++) coutp[i] = ptro[i];

  // Sortie  
  K_CONNECT::cleanUnreferencedFacesNGon(*cout);
  PyObject* tpl = K_ARRAY::buildArray3(*fcout, varStringC, *cout, "NGON", api);
  RELEASESHAREDU(arrayC, fc, cn);
  delete fcout; delete cout;
  return tpl;
}  

//=============================================================================
/* Rm ghost cells to a NGON array */
//=============================================================================
PyObject* K_CONVERTER::rmGhostCellsNGonBoth(PyObject* self, PyObject* args)
{
  PyObject *arrayN, *arrayC;
  E_Int depth;
  if (!PYPARSETUPLE_(args, OO_ I_, &arrayN, &arrayC, &depth)) return NULL;

  if (depth < 1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "rmGhostCells: depth must be > 0.");
    return NULL;
  }
  E_Int ni, nj, nk;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(arrayN, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(arrayN,f);
    PyErr_SetString(PyExc_TypeError, 
                    "rmGhostCells: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType,"NGON") != 0 && strcmp(eltType,"NGON*") != 0 ) 
  {
    RELEASESHAREDU(arrayN,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "rmGhostCells: array must be a NGON.");
    return NULL;
  }
    
  E_Int nic, njc, nkc;
  char* varStringC; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int resc = K_ARRAY::getFromArray3(arrayC, varStringC, fc, nic, njc, nkc, 
                                      cnc, eltTypec);
  if (resc != 2)
  {
    if (resc == 1) RELEASESHAREDS(arrayC, fc);
    PyErr_SetString(PyExc_TypeError, 
                    "rmGhostCells: array is invalid.");
    RELEASESHAREDU(arrayN, f, cn); return NULL;
  }
  if (strcmp(eltTypec, "NGON*") != 0 ) 
  {
    RELEASESHAREDU(arrayN,f,cn);
    RELEASESHAREDU(arrayC,fc,cnc);
    PyErr_SetString(PyExc_TypeError,
                    "rmGhostCells: array must be NGON*.");
    return NULL;
  }

  FldArrayI cFE; K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int* ptr = cn->begin();
  E_Int sizeFN = ptr[1]; ptr += sizeFN+2;
  E_Int ne = ptr[0];
  E_Int sizeEF = ptr[1]; ptr += 2;
  
  E_Int nfe, e1, e2, face;
  FldArrayI co(sizeEF); // elements interieurs
  E_Int nint=0; // nbre d'elts interieurs
  E_Int sizeco = 0; // taille de la connectivite d'elements interieurs
  E_Int* ptro = co.begin();
  E_Int sizengbrs, indv;
  //-----
  // Compute tag=1 : elts to be removed
  // 1st layer: external face  - Upper layers: neighbours of an elt tagged 1
  FldArrayI tag(ne); tag.setAllValuesAtNull();
  for (E_Int i = 0; i < ne; i++)
  {
    nfe = ptr[0];
    for (E_Int n = 1; n <= nfe; n++)
    {
      face = ptr[n]-1; 
      e1 = cFE1[face]; e2 = cFE2[face];
      if (e1 == 0 && e2 > 0)  tag[i] = 1;
      else if (e1 > 0 && e2 == 0) tag[i] = 1;
    }
    ptr += nfe+1;
  }


  if (depth > 1)
  {
    vector < vector<E_Int> > cEEN(ne);
    K_CONNECT::connectFE2EENbrs(cFE, cEEN);  
    for (E_Int d = 2; d <= depth; d++)
    {
      ptr = cn->begin()+sizeFN+4;
      for (E_Int i = 0; i < ne; i++)
      {
        nfe = ptr[0];
        vector<E_Int>& ngbrs = cEEN[i];
        sizengbrs = ngbrs.size();
        if (tag[i] == d-1)
        {
          for (E_Int ng = 0; ng < sizengbrs; ng++)
          {
            indv = ngbrs[ng];
            if (tag[indv] == 0) tag[indv] = d;
          }
        }
        ptr += nfe+1;
      }
    }
  }
 
  //-----
  // Select interior elements
  ptr = cn->begin()+sizeFN+4;
  for (E_Int i = 0; i < ne; i++)
  {
    nfe = ptr[0];
    if (tag[i] == 0) // recopie l'element 
    {
      ptro[0] = nfe; 
      
      for (E_Int n = 1; n <= nfe; n++) 
        ptro[n] = ptr[n];  
      
      ptro += nfe+1; sizeco += nfe+1; nint++;
    }
    ptr += nfe+1;
  }
  co.reAlloc(sizeco); ptro = co.begin();

  // Les points ajoutes sont des points crees a la fin par addGhostCells
  // si on a conserve l ordre de addGhostCells, pas besoin de remodifier la connectivite Faces/Noeuds
  E_Int api = fc->getApi();
  E_Int nfld = fc->getNfld();
  FldArrayF* fcout = new FldArrayF(nint,nfld);
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fpc = fcout->begin(eq);
    E_Float* fp = fc->begin(eq);
    for (E_Int i = 0; i < nint; i++)
    {
      fpc[i] = fp[i];
    }
  }

  // Cree la nouvelle connectivite complete
  FldArrayI* cout = new FldArrayI(sizeFN+4+sizeco);
  E_Int* coutp = cout->begin();
  ptr = cn->begin();
  for (E_Int i = 0; i < sizeFN+2; i++) coutp[i] = ptr[i];
  coutp += sizeFN+2;
  coutp[0] = nint;
  coutp[1] = sizeco; coutp += 2;
  for (E_Int i = 0; i < sizeco; i++) coutp[i] = ptro[i];

  // Sortie  
  K_CONNECT::cleanUnreferencedFacesNGon(*cout); 

  PyObject* l = PyList_New(0);
  PyObject* tpl = K_ARRAY::buildArray3(*f, varString, *cout, eltType, api);
  PyList_Append(l, tpl); Py_DECREF(tpl);
  //array en centres
  tpl = K_ARRAY::buildArray3(*fcout, varStringC, *cout, "NGON", api);
  PyList_Append(l, tpl); Py_DECREF(tpl);  delete fcout; delete cout;

  RELEASESHAREDU(arrayN, f, cn); RELEASESHAREDU(arrayC, fc, cnc);
  return l;
}  


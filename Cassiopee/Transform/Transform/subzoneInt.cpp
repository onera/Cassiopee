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

// subzoneStruct with a list of interface indices

# include "transform.h"
using namespace std;
using namespace K_FLD;

// ============================================================================
/* Subzone a structured mesh using a list of interface indices
  returns a QUAD mesh */
// ============================================================================
PyObject* K_TRANSFORM::subzoneStructInt(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* listOfInterfaces;
  if (!PyArg_ParseTuple(args,"OO", &array, &listOfInterfaces)) return NULL;
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType); 
  if ( res == 1 ) ;
  else if (res == 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: cannot be used on an unstructured array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: unknown type of array.");
    return NULL;
  }
  FldArrayI intIndices;
  E_Int ok = K_ARRAY::getFromList(listOfInterfaces, intIndices);
  if (ok == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: 2nd argument must be an integer list or a numpy.");
    RELEASESHAREDS(array,f); return NULL;
  }
  E_Int n = intIndices.getSize();
  E_Int* intIndicesp = intIndices.begin();

  char newEltType[256];
  E_Int ni1 = K_FUNC::E_max(1,ni-1); 
  E_Int nj1 = K_FUNC::E_max(1,nj-1); 
  E_Int nk1 = K_FUNC::E_max(1,nk-1);
  E_Int ni1nj1 = ni1*nj1;
  E_Int ninti  = ni*nj1*nk1;
  E_Int nintj  = ni1*nj*nk1; 
  E_Int nintij = ninti+nintj;
  E_Int nfld = f->getNfld();  
  PyObject* tpl;
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;  
  E_Int i,j,k,ind, indint, incnode;
  E_Int incdir1, incdir2;
  E_Int ni1nj = ni1*nj;
  E_Int ninj = ni*nj;
  E_Int ninj1 = ni*nj1;
  E_Int indcell=0; E_Int nov = 0;

  if (nj == 1 && nk == 1)
  {
    strcpy(newEltType, "NODE");
    FldArrayF* fnodes = new FldArrayF(n,nfld);// dimension max
    FldArrayI* connect = new FldArrayI(0,1);// nb interfaces
    for (E_Int noint = 0; noint < n; noint++)
    {
      indint = intIndicesp[noint];
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        (*fnodes)(nov,eq) = 0.5*((*f)(indint,eq)+(*f)(indint+1,eq));
      }
      nov++;
    }
    tpl = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    delete fnodes; delete connect;
  }
  else if (nk == 1)
  {
    strcpy(newEltType, "BAR");
    FldArrayF* fnodes = new FldArrayF(n*2,nfld);// dimension max
    FldArrayI* connect = new FldArrayI(n,2);// nb interfaces
    E_Int* cn1 = connect->begin(1);
    E_Int* cn2 = connect->begin(2);
    for (E_Int noint = 0; noint < n; noint++)
    {
      indint = intIndicesp[noint];
      if (indint < ninti) // i-interface 
      {
        incnode = ni;
        j = indint/ni;
        i = indint-j*ni;
      }
      else // j-interface
      {
        incnode = 1;
        j = (indint-ninti)/ni1;
        i = (indint-ninti)-j*ni1;
      } 
      
      ind = i + j*ni;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        (*fnodes)(nov,eq) = (*f)(ind,eq);
        (*fnodes)(nov+1,eq) = (*f)(ind+incnode,eq);
      }
      cn1[indcell] = nov+1;
      cn2[indcell] = nov+2;
      nov+=2;
      indcell++;
    }
    
    if (posx != 0 && posy != 0 && posz != 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                   1.e-12, newEltType, 
                                   *fnodes, *connect);
    tpl = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    delete fnodes; delete connect;
  }
  else if (ni > 1 && nj > 1 && nk > 1)
  {  
    strcpy(newEltType, "QUAD");
    FldArrayF* fnodes = new FldArrayF(n*4,nfld);// dimension max
    FldArrayI* connect = new FldArrayI(n,4);// nb interfaces
    E_Int* cn1 = connect->begin(1);
    E_Int* cn2 = connect->begin(2);
    E_Int* cn3 = connect->begin(3);
    E_Int* cn4 = connect->begin(4);

    indcell=0; nov = 0;
    for (E_Int noint = 0; noint < n; noint++)
    {
      indint = intIndicesp[noint];
      if (indint < ninti) // i-interface 
      {
        //indint = i+j*ni+k*ni*nj1;
        k = indint/(ninj1);
        j = (indint-k*ninj1)/ni;
        i = indint-j*ni-k*ninj1;
        incdir1 = ni; incdir2 = ninj;
      }
      else if ( indint < nintij)// j-interface
      {
        //indint = i+j*ni1+k*ni1*nj+ninti;
        k = (indint-ninti)/(ni1nj);
        j = (indint-ninti-k*ni1nj)/ni1;
        i = indint-ninti-k*ni1nj-j*ni1;
        incdir1 = 1; incdir2 = ninj;
      } 
      else // k-interface
      {
        //indint = i+j*ni1+k*ni1*nj1+nintij;
        k = (indint-nintij)/(ni1nj1);
        j = (indint-nintij-k*ni1nj1)/ni1;
        i = indint-nintij-k*ni1nj1-j*ni1;
        incdir1 = 1; incdir2 = ni;
      }
      ind = i + j*ni + k*ninj;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        (*fnodes)(nov,eq) = (*f)(ind,eq);
        (*fnodes)(nov+1,eq) = (*f)(ind+incdir1,eq);
        (*fnodes)(nov+2,eq) = (*f)(ind+incdir1+incdir2,eq);
        (*fnodes)(nov+3,eq) = (*f)(ind+incdir2,eq);
      }
      cn1[indcell] = nov+1;
      cn2[indcell] = nov+2;
      cn3[indcell] = nov+3;
      cn4[indcell] = nov+4;
      nov+=4; indcell++;
    }
    
    if (posx != 0 && posy != 0 && posz != 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                   1.e-12, newEltType, 
                                   *fnodes, *connect);
    tpl = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    delete fnodes; delete connect;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: 1D zones must be nk=1 and 2D zones (nj=1,nk=1).");
    RELEASESHAREDS(array,f);  return NULL;
  }
  RELEASESHAREDS(array, f);
  return tpl;
}     
// ============================================================================
/* Subzone a structured mesh using a list of interface indices
  returns a QUAD mesh */
// ============================================================================
PyObject* K_TRANSFORM::subzoneStructIntBoth(PyObject* self, PyObject* args)
{
  PyObject *arrayN, *arrayC, *listOfInterfaces;
  if (!PyArg_ParseTuple(args,"OOO", &arrayN, &arrayC, &listOfInterfaces)) return NULL;

  FldArrayI intIndices;
  E_Int ok = K_ARRAY::getFromList(listOfInterfaces, intIndices);
  if (ok == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructIntBoth: 2nd argument must be an integer list or a numpy.");
    return NULL;
  }
  E_Int n = intIndices.getSize();
  E_Int* intIndicesp = intIndices.begin();

  // Check array of nodes
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(arrayN, varString, f, ni, nj, nk, cn, eltType); 
  if (res == 1);
  else if (res == 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: cannot be used on an unstructured array.");
    RELEASESHAREDU(arrayN, f, cn); return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: unknown type of array.");
    return NULL;
  }
  // Check array of centers 
  E_Int nic, njc, nkc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  res = K_ARRAY::getFromArray3(arrayC, varStringc, fc, nic, njc, nkc, cnc, eltTypec); 
  if (res == 1);
  else if (res == 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: cannot be used on an unstructured array.");
    RELEASESHAREDS(arrayN,f); RELEASESHAREDU(arrayC, fc, cnc); return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: unknown type of array.");
    RELEASESHAREDS(arrayN,f); return NULL;
  }

  char newEltType[256];
  E_Int ni1 = K_FUNC::E_max(1,ni-1); 
  E_Int nj1 = K_FUNC::E_max(1,nj-1); 
  E_Int nk1 = K_FUNC::E_max(1,nk-1);
  E_Int ni1nj1 = ni1*nj1;
  E_Int ninti  = ni*nj1*nk1;
  E_Int nintj  = ni1*nj*nk1; 
  E_Int nintij = ninti+nintj;
  E_Int i,j,k,ind, indcell, indint, incnode;
  E_Int incdir1, incdir2;
  E_Int ni1nj = ni1*nj;
  E_Int ninj = ni*nj;
  E_Int ninj1 = ni*nj1;

  E_Int nfld = f->getNfld();  
  E_Int nfldc = fc->getNfld();  
  PyObject *tplN, *tplC;
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;  
  PyObject* l = PyList_New(0);
  E_Int nov = 0;
  E_Int noet = 0;
  if (nj == 1 && nk == 1)
  {
    strcpy(newEltType, "NODE");
    FldArrayF* fnodes = new FldArrayF(n,nfld);// dimension max
    FldArrayF* fcenters = new FldArrayF(n,nfldc);
    FldArrayI* connect = new FldArrayI(0,1);// nb interfaces
    for (E_Int noint = 0; noint < n; noint++)
    {
      indint = intIndicesp[noint];
      i = indint;
      for (E_Int eq = 1; eq <= nfld; eq++)      
        (*fnodes)(nov,eq) = (*f)(i,eq);
      if (i == 0)
      {
        for (E_Int eq = 1; eq <= nfldc; eq++)
          (*fcenters)(nov,eq) = (*fc)(i,eq); 
      }
      else if (i == ni-1)
      {
        for (E_Int eq = 1; eq <= nfldc; eq++)
          (*fcenters)(nov,eq) = (*fc)(i-1,eq); 
      }
      else
      {
        for (E_Int eq = 1; eq <= nfldc; eq++)
          (*fcenters)(nov,eq) = 0.5*((*fc)(i,eq)+(*fc)(i-1,eq));
      }
      nov++;
    }
    tplN = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    PyList_Append(l,tplN); Py_DECREF(tplN);
    tplC = K_ARRAY::buildArray(*fcenters, varStringc, *connect, -1, newEltType, true);
    PyList_Append(l,tplC); Py_DECREF(tplC);
    delete fnodes; delete fcenters; delete connect; 
  }
  else if (nk == 1)
  {
    strcpy(newEltType, "BAR");
    FldArrayF* fnodes = new FldArrayF(n*2,nfld);// dimension max
    FldArrayF* fcenters = new FldArrayF(n,nfldc);
    FldArrayI* connect = new FldArrayI(n,2);// nb interfaces
    E_Int* cn1 = connect->begin(1);
    E_Int* cn2 = connect->begin(2);
    for (E_Int noint = 0; noint < n; noint++)
    {
      indint = intIndicesp[noint];
      if (indint < ninti) // i-interface 
      {
        incnode = ni;
        j = indint/ni;
        i = indint-j*ni;
      }
      else // j-interface
      {
        incnode = 1;
        j = (indint-ninti)/ni1;
        i = (indint-ninti)-j*ni1;
      } 
      
      ind = i + j*ni;
      indcell = i+ j*ni1;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        (*fnodes)(nov,eq) = (*f)(ind,eq);
        (*fnodes)(nov+1,eq) = (*f)(ind+incnode,eq);
      }
      for (E_Int eq = 1; eq <= nfldc; eq++)      
        (*fcenters)(noet,eq) = (*fc)(indcell,eq);
      
      cn1[noet] = nov+1;
      cn2[noet] = nov+2;
      nov+=2; noet++;
    }
    
    if (posx != 0 && posy != 0 && posz != 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                   1.e-12, newEltType, 
                                   *fnodes, *connect);
    tplN = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    PyList_Append(l,tplN); Py_DECREF(tplN);
    tplC = K_ARRAY::buildArray(*fcenters, varStringc, *connect, -1, newEltType, true);
    PyList_Append(l,tplC); Py_DECREF(tplC);
    delete fnodes; delete fcenters; delete connect; 
  }
  else if (ni > 1 && nj > 1 && nk > 1)
  {  
    strcpy(newEltType, "QUAD");
    FldArrayF* fnodes = new FldArrayF(n*4,nfld);// dimension max
    FldArrayF* fcenters = new FldArrayF(n,nfldc);
    FldArrayI* connect = new FldArrayI(n,4);// nb interfaces
    E_Int* cn1 = connect->begin(1);
    E_Int* cn2 = connect->begin(2);
    E_Int* cn3 = connect->begin(3);
    E_Int* cn4 = connect->begin(4);

    indcell=0; nov = 0;
    for (E_Int noint = 0; noint < n; noint++)
    {
      indint = intIndicesp[noint];
      if (indint < ninti) // i-interface 
      {
        //indint = i+j*ni+k*ni*nj1;
        k = indint/(ninj1);
        j = (indint-k*ninj1)/ni;
        i = indint-j*ni-k*ninj1;
        incdir1 = ni; incdir2 = ninj;
      }
      else if ( indint < nintij)// j-interface
      {
        //indint = i+j*ni1+k*ni1*nj+ninti;
        k = (indint-ninti)/(ni1nj);
        j = (indint-ninti-k*ni1nj)/ni1;
        i = indint-ninti-k*ni1nj-j*ni1;
        incdir1 = 1; incdir2 = ninj;
      } 
      else // k-interface
      {
        //indint = i+j*ni1+k*ni1*nj1+nintij;
        k = (indint-nintij)/(ni1nj1);
        j = (indint-nintij-k*ni1nj1)/ni1;
        i = indint-nintij-k*ni1nj1-j*ni1;
        incdir1 = 1; incdir2 = ni;
      }
      ind = i + j*ni + k*ninj;
      indcell = i + j*ni1 + k*ni1nj1;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        (*fnodes)(nov,eq) = (*f)(ind,eq);
        (*fnodes)(nov+1,eq) = (*f)(ind+incdir1,eq);
        (*fnodes)(nov+2,eq) = (*f)(ind+incdir1+incdir2,eq);
        (*fnodes)(nov+3,eq) = (*f)(ind+incdir2,eq);
      }
      for (E_Int eq = 1; eq <= nfldc; eq++)
        (*fcenters)(noet,eq) = (*fc)(indcell,eq);
      cn1[noet] = nov+1;
      cn2[noet] = nov+2;
      cn3[noet] = nov+3;
      cn4[noet] = nov+4;
      nov+=4; noet++;
    }
    
    if (posx != 0 && posy != 0 && posz != 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                   1.e-12, newEltType, 
                                   *fnodes, *connect);
    tplN = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    PyList_Append(l,tplN); Py_DECREF(tplN);
    tplC = K_ARRAY::buildArray(*fcenters, varStringc, *connect, -1, newEltType, true);
    PyList_Append(l,tplC); Py_DECREF(tplC);
    delete fnodes; delete fcenters; delete connect; 
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneStructInt: 1D zones must be nk=1 and 2D zones (nj=1,nk=1).");
    RELEASESHAREDS(arrayN,f);  return NULL;
  }
  RELEASESHAREDS(arrayN, f); RELEASESHAREDS(arrayC, fc);
  return l;
}     


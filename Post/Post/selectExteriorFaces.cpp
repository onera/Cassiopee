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

// selectExteriorFaces

#include <stdio.h>
#include <string.h>
#include "post.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Selectionne les facettes exterieures d'un array */
// ============================================================================
PyObject* K_POST::selectExteriorFaces(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* indices;
  if (!PyArg_ParseTuple(args, "OO", &array, &indices)) return NULL;
  
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  
  PyObject* tpl = NULL;
  if (res == 1)
  {
    tpl = exteriorFacesStructured(varString, *f, ni, nj, nk, indices);
    RELEASESHAREDS(array, f);
  }
  else if (res == 2)
  {
    if (strcmp(eltType, "NGON") == 0)
    {
      if ((*cn)[2] == 2) tpl = selectExteriorFacesNGon2D(varString, *f, *cn, indices);
      else  tpl = selectExteriorFacesNGon3D(varString, *f, *cn, indices);
    }
    else tpl = exteriorFacesBasic(varString, *f, *cn, eltType, indices);
    RELEASESHAREDU(array, f, cn); 
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorFaces: invalid array.");
    RELEASESHAREDB(res, array, f, cn); 
  } 
  return tpl;
}
//=============================================================================
PyObject* K_POST::exteriorFacesStructured(char* varString, FldArrayF& f, 
                                          E_Int ni, E_Int nj, E_Int nk, 
                                          PyObject* indices)
{
  E_Int nfld = f.getNfld();
  PyObject* tpl = NULL;
  bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;
  PyObject* indir = NULL;
  E_Int* indirp = NULL;
  char newEltType[256];

  // 1D arrays
  if ((ni == 1 && nj == 1) || (nj == 1 && nk == 1) || (ni == 1 && nk == 1))
  {
    FldArrayF* fnodes = new FldArrayF(2, nfld);
    FldArrayI* connect = new FldArrayI(0, 1);
    E_Int s = f.getSize();
    for (E_Int i = 1; i <= nfld; i++)
    {
      (*fnodes)(0,i) = f(0,i);
      (*fnodes)(1,i) = f(s-1,i);
    }
    tpl = K_ARRAY::buildArray(*fnodes, varString, 
                              *connect, -1, "NODE");
    if (boolIndir == true)
    {
      indir = K_NUMPY::buildNumpyArray(2, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);
      indirp[0] = 1; indirp[1] = s;
      PyList_Append(indices, indir);  Py_DECREF(indir);
    }

    delete fnodes; delete connect;
    return tpl;
  }
  else if (ni == 1 || nj == 1 || nk == 1)//arrays 2D
  {
    strcpy(newEltType, "BAR");

    E_Int p=0, q=0;
    if (nk == 1) {p = ni; q = nj;}
    if (ni == 1) {p = nj; q = nk;}
    if (nj == 1) {p = ni; q = nk;}
    E_Int ninti = 2*(q-1); E_Int nintj = 2*(p-1);
    E_Int nfacesExt = ninti + nintj;
    if (nfacesExt==0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "exteriorFaces: set of exterior faces is empty.");
      return NULL;
    }

    FldArrayF* fnodes = new FldArrayF(nfacesExt,nfld);
    FldArrayI* connect = new FldArrayI(nfacesExt,2);      
    E_Int* cn1 = connect->begin(1);
    E_Int* cn2 = connect->begin(2);
    E_Int n = 0;

    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fnv = fnodes->begin(nv);
      n = 0;
      // Border j=1, i=1..p-1
      for (E_Int i = 0; i < p-1; i++)
      { 
        fnv[n] = fv[i]; 
        n++;
      }
      // Border i=p,j=0..q-1
      for (E_Int j = 0; j < q-1; j++)
      { 
        fnv[n] = fv[p-1+j*p];
        n++;
      }
      // Border j=jmax,i=p,0
      for (E_Int i = p-1; i > 0; i--)
      { 
        fnv[n] = fv[i+(q-1)*p]; 
        n++;
      }
      // Border i=1,j=q,1
      for (E_Int j = q-1; j > 0; j--)
      { 
        fnv[n] = fv[j*ni]; 
        n++;
      }
    }
    n = 1;
    for (E_Int noe=0; noe < nfacesExt-1; noe++)
    {
      cn1[noe]=n; cn2[noe]=n+1; n++;
    }
    if (boolIndir == true)
    {
      E_Int nbinti = p*(q-1);
      //E_Int nbintj = q*(p-1);
      indir = K_NUMPY::buildNumpyArray(nfacesExt, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);        
      E_Int noind = 0;
      // j=1 interfaces
      for (E_Int i = 0; i < p-1; i++)
      {
        indirp[noind] = nbinti+i;
        noind++;
      }
      // i=imax interfaces
      for (E_Int j = 0; j < q-1; j++)
      {
        indirp[noind] = (p-1)+ j*p;
        noind++;
      }

      // j=jmax interfaces
      for (E_Int i = 0; i < p-1; i++)
      { 
        indirp[noind] = nbinti+(p-2-i)+(q-1)*(p-1);
        noind++;
      }

      // i=0 interface
      for (E_Int j = 0; j <q-1; j++)
      { 
        indirp[noind] = (q-2-j)*p;
        noind++;
      }
      PyList_Append(indices, indir);  Py_DECREF(indir);
    }
    cn1[nfacesExt-1]=nfacesExt; cn2[nfacesExt-1]=1;
    PyObject* tpl = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    delete fnodes; delete connect;
    return tpl;
  }
  else
  {
    strcpy(newEltType, "QUAD");
    E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
    //E_Int ni1nj1 = ni1*nj1;
    E_Int nPtsExt = 2*(ni*nj+ni*nk+nj*nk);
    E_Int nfacesExt = 2*(ni1*nj1+nj1*nk1+ni1*nk1);
    FldArrayF* fnodes = new FldArrayF(nPtsExt, nfld);
    FldArrayI* connect = new FldArrayI(nfacesExt, 4);
    E_Int ninti = 0, nintj = 0, nintij = 0;
    if (boolIndir == true)
    {
      ninti  = ni*nj1*nk1;
      nintj  = ni1*nj*nk1; 
      nintij = ninti+nintj;
      indir  = K_NUMPY::buildNumpyArray(nfacesExt, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);       
    }
    E_Int ns = 0; E_Int n = 0; 
    E_Int noint = 0;
    E_Int e = 0; E_Int ind, indint;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      E_Float* fv = f.begin(nv);
      E_Float* fnv = fnodes->begin(nv);
      for (n = 0; n < ni*nj; n++) fnv[n] = fv[n];
    }
    E_Int* cn1 = connect->begin(1);
    E_Int* cn2 = connect->begin(2);
    E_Int* cn3 = connect->begin(3);
    E_Int* cn4 = connect->begin(4);
    
    for (E_Int j = 0; j < nj1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        ind = ns + i+j*ni + 1;
        cn1[e] = ind;
        cn2[e] = ind+1;
        cn3[e] = ind+1+ni;
        cn4[e] = ind+ni;
        e++;
      }

    ns = n;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      n = ns;
      E_Float* fv = f.begin(nv);
      E_Float* fnv = fnodes->begin(nv);
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        { 
          fnv[n] = fv[i+j*ni+(nk-1)*ni*nj]; n++;
        }
    }
    for (E_Int j = 0; j < nj1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        ind = ns + i+j*ni + 1;
        cn1[e] = ind;
        cn2[e] = ind+1;
        cn3[e] = ind+1+ni;
        cn4[e] = ind+ni; 
        e++;
      }

    if (boolIndir == true)
    {
      // k = 0 followed by k=nk
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          indint = i+j*ni1+nintij;
          indirp[noint] = indint;
          noint++;
        }
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          indint = i+j*ni1+nk1*ni1*nj1+nintij;
          indirp[noint] = indint;
          noint++;
        }
    }
    ns = n;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      n = ns;
      E_Float* fv = f.begin(nv);
      E_Float* fnv = fnodes->begin(nv);
      for (E_Int j = 0; j < nk; j++)
        for (E_Int i = 0; i < ni; i++)
        { 
          fnv[n] = fv[i+j*ni*nj]; n++;
        }
    }
    for (E_Int j = 0; j < nk1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        ind = ns + i+j*ni + 1;
        cn1[e] = ind;
        cn2[e] = ind+1;
        cn3[e] = ind+1+ni;
        cn4[e] = ind+ni;
        e++;
      }
    
    ns = n;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      n = ns;
      E_Float* fv = f.begin(nv);
      E_Float* fnv = fnodes->begin(nv);
      for (E_Int j = 0; j < nk; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          fnv[n] = fv[i+(nj-1)*ni+j*ni*nj]; n++;
        }
    }
    for (E_Int j = 0; j < nk1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        ind = ns + i+j*ni + 1;
        cn1[e] = ind;
        cn2[e] = ind+1;
        cn3[e] = ind+1+ni;
        cn4[e] = ind+ni; 
        e++;
      } 
    if (boolIndir == true)
    {
      // j = 0 followed by j=nj
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int i = 0; i < ni1; i++)
        {
          indint = i+k*ni1*nj+ninti;
          indirp[noint] = indint;
          noint++;
        }
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int i = 0; i < ni1; i++)
        {
          indint = i+nj1*ni1+k*ni1*nj+ninti;
          indirp[noint] = indint;
          noint++;
        }
    }
    ns = n;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      n = ns;
      E_Float* fv = f.begin(nv);
      E_Float* fnv = fnodes->begin(nv);
      for (E_Int j = 0; j < nk; j++)
        for (E_Int i = 0; i < nj; i++)
        { 
          fnv[n] = fv[i*ni+j*ni*nj]; n++;
        }
    }
    for (E_Int j = 0; j < nk1; j++)
      for (E_Int i = 0; i < nj1; i++)
      {
        ind = ns + i+j*nj + 1; 
        cn1[e] = ind;
        cn2[e] = ind+1;
        cn3[e] = ind+1+nj;
        cn4[e] = ind+nj; 
        e++;
      }

    ns = n;
    for (E_Int nv = 1; nv <= nfld; nv++)
    {
      n = ns;
      E_Float* fv = f.begin(nv);
      E_Float* fnv = fnodes->begin(nv);
      for (E_Int j = 0; j < nk; j++)
        for (E_Int i = 0; i < nj; i++)
        { 
          fnv[n] = fv[ni-1+i*ni+j*ni*nj]; n++;
        }
    }
    for (E_Int j = 0; j < nk1; j++)
      for (E_Int i = 0; i < nj1; i++)
      {
        ind = ns + i+j*nj + 1; 
        cn1[e] = ind;
        cn2[e] = ind+1;
        cn3[e] = ind+1+nj;
        cn4[e] = ind+nj;
        e++;
      }

    if (boolIndir == true)
    {
      // i = 0 followed by i=ni
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        {
          indint = j*ni+k*ni*nj1;
          indirp[noint] = indint;
          noint++;
        }
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        {
          indint = ni1+j*ni+k*ni*nj1;
          indirp[noint] = indint;
          noint++;
        }
    }      

    E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;
    if (posx != 0 && posy != 0 && posz != 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                   1.e-12, "QUAD", 
                                   *fnodes, *connect);

    if ( boolIndir==true)
    {
      PyList_Append(indices, indir);  Py_DECREF(indir);
    }
    PyObject* tpl = K_ARRAY::buildArray(*fnodes, varString, *connect, -1, newEltType);
    delete fnodes; delete connect;
    return tpl;
  }
  return NULL;
}

//=============================================================================
// exteriorFaces pour non-structure TRI, QUAD, TETRA, HEXA (Basic elts)
//=============================================================================
PyObject* K_POST::exteriorFacesBasic(char* varString, FldArrayF& f, 
                                     FldArrayI& cn, char* eltType,
                                     PyObject* indices)
{
  bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;
  char elttypeout[256];
  if (strcmp(eltType, "TRI") == 0 || strcmp(eltType, "QUAD") == 0)
    strcpy(elttypeout, "BAR");
  else if (strcmp(eltType, "TETRA") == 0)
    strcpy(elttypeout, "TRI");
  else if (strcmp(eltType, "HEXA") == 0)
    strcpy(elttypeout, "QUAD");
  else if (strcmp(eltType,"BAR") == 0 )
    strcpy(elttypeout,"NODE");
	
  E_Int nfld = f.getNfld();

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;

  //-- nouvelle version --
  if (strcmp(eltType, "BAR") == 0)
  {
    FldArrayF* fnodes = new FldArrayF(2, nfld);
    FldArrayI* connect = new FldArrayI(0, 1);
    E_Int net = cn.getSize();
    E_Int ind0 = cn(0,1)-1; E_Int ind1 = cn(net-1,2)-1;
    for (E_Int i = 1; i <= nfld; i++)
    {
      (*fnodes)(0,i) = f(ind0,i);
      (*fnodes)(1,i) = f(ind1,i);
    }
    if (posx != 0 && posy != 0 && posz != 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                   1.e-12, elttypeout, 
                                   *fnodes, *connect);
    E_Int nExtNodes = fnodes->getSize();
    if (nExtNodes == 1 ) // c est ferme : un doublon
    {
      PyErr_SetString(PyExc_ValueError,
                      "exteriorFaces: exterior face set is empty.");
      delete fnodes; delete connect; 
      return NULL;
    }
    
    PyObject* tpl = K_ARRAY::buildArray(*fnodes, varString, 
                                        *connect, -1, "NODE");
    delete fnodes; delete connect;
    if (indices != Py_None)
    { 
      PyObject* indir = K_NUMPY::buildNumpyArray(nExtNodes, 1, 1, 0);
      E_Int* indirp = K_NUMPY::getNumpyPtrI(indir);
      if ( nExtNodes > 0) {indirp[0] = ind0+1; indirp[1] = ind1+1;}
      PyList_Append(indices, indir);  Py_DECREF(indir);
    }
    return tpl;
  }
  
  FldArrayI* connect = new FldArrayI();
  short ok = exteriorFacesBasic3(f, cn, eltType, *connect, boolIndir, indices);
  
  if (ok == 0) // erreur interne
  {
    delete connect;
    return NULL;
  }
	
  if (posx != 0 && posy != 0 && posz != 0)
    K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                 1.e-12, elttypeout, 
                                 f, *connect);
  PyObject* tpl;
  if (connect->getSize() == 0)
  {
    delete connect;
    PyErr_SetString(PyExc_ValueError,
                    "exteriorFaces: exterior face set is empty.");
    tpl = NULL;
  }
  else
  {
    tpl = K_ARRAY::buildArray(f, varString, 
                              *connect, -1, elttypeout);
    delete connect;
  }

  return tpl;
}
//=============================================================================
// exteriorFaces pour les array non structures a elt basique
// Cet algorithme est fonde sur la topologie: une facette est externe
// si elle n'a qu'un seul element voisin (autre version)
//=============================================================================
short K_POST::exteriorFacesBasic3(FldArrayF& f, FldArrayI& cn,
                                  char* eltType,
                                  FldArrayI& cnn, 
                                  bool boolIndir, PyObject* indices)
{
  // Determination des elements exterieurs
  E_Int nelts = cn.getSize();
  vector< vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType, f.getSize(), cn, cEEN);

  // Nombre de voisins pour un elt interieur
  unsigned int nvoisins = 0; // nbre de voisin si interieur  
  unsigned int nfaces = 0; // nbre de faces par elements
  unsigned int nof = 0; // nbre de noeuds par face

  E_Int face[6][4];
  if (strcmp(eltType, "BAR") == 0) 
  { 
    nvoisins = 2;
    nfaces = 2; nof = 1;
    face[0][0] = 1; face[1][0] = 2;
  }
  else if (strcmp(eltType, "QUAD") == 0) 
  {
    nvoisins = 4;
    nfaces = 4; nof = 2;
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 4;
    face[3][0] = 4; face[3][1] = 1;
  }
  else if (strcmp(eltType, "TRI") == 0) 
  {
    nvoisins = 3;
    nfaces = 3; nof = 2;
    face[0][0] = 1; face[0][1] = 2;
    face[1][0] = 2; face[1][1] = 3;
    face[2][0] = 3; face[2][1] = 1;
  }
  else if (strcmp(eltType, "HEXA") == 0) 
  {
    nvoisins = 6;
    nfaces = 6; nof = 4;
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; face[0][3] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 6; face[1][3] = 5;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 7; face[2][3] = 6;
    face[3][0] = 3; face[3][1] = 4; face[3][2] = 8; face[3][3] = 7;
    face[4][0] = 1; face[4][1] = 5; face[4][2] = 8; face[4][3] = 4;
    face[5][0] = 5; face[5][1] = 6; face[5][2] = 7; face[5][3] = 8;
  }
  else if (strcmp(eltType, "TETRA") == 0) 
  {
    nvoisins = 4;
    nfaces = 4; nof = 3;
    face[0][0] = 1; face[0][1] = 3; face[0][2] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 4;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 4;
    face[3][0] = 3; face[3][1] = 1; face[3][2] = 4;
  }
  else if (strcmp(eltType, "PYRA") == 0) 
  {
    nvoisins = 5;
    nfaces = 5; nof = 3; // 2 TRIs pour la base
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3;
    face[1][0] = 3; face[1][1] = 2; face[1][2] = 1;
    face[2][0] = 1; face[2][1] = 2; face[2][2] = 5; 
    face[3][0] = 2; face[3][1] = 3; face[3][2] = 5;
    face[4][0] = 3; face[4][1] = 4; face[4][2] = 5;
    face[5][0] = 4; face[5][1] = 1; face[5][2] = 5;
  }
  else if (strcmp(eltType, "PENTA") == 0) 
  {
    nvoisins = 5;
    nfaces = 5; nof = 4; // TRI degen
    face[0][0] = 1; face[0][1] = 2; face[0][2] = 5; face[0][3] = 4;
    face[1][0] = 2; face[1][1] = 3; face[1][2] = 6; face[1][3] = 5;
    face[2][0] = 3; face[2][1] = 1; face[2][2] = 4; face[2][3] = 6;
    face[3][0] = 1; face[3][1] = 3; face[3][2] = 2; face[3][3] = 2;
    face[4][0] = 4; face[4][1] = 5; face[4][2] = 6; face[4][3] = 6;
  }
  
  E_Int nthreads = __NUMTHREADS__;
  E_Int net = nelts/nthreads+1;
  FldArrayI** cnnls = new FldArrayI* [nthreads];
  E_Int* nes = new E_Int [nthreads];
  E_Int* prev = new E_Int [nthreads];
  vector< vector<E_Int> > vectOfFaces(nthreads);

#pragma omp parallel default(shared)
  {
  E_Int ithread = __CURRENT_THREAD__;
  E_Int* indir = new E_Int [net];
  vector<E_Int>& vectOfFacesL = vectOfFaces[ithread];

  E_Int ne = 0;
  // nbre d'elements exterieurs + indirection
#pragma omp for
  for (E_Int e = 0; e < nelts; e++)
  {
    if (cEEN[e].size() != nvoisins)
    {
      indir[ne] = e; ne++;
    }
  }

  // determine les faces exterieures
  unsigned int tag;
  E_Int s; E_Int e, ev;
  E_Int nn = 0; // nbre de faces ext
  FldArrayI* cnnl = new FldArrayI(ne*nfaces,nof);
  cnnls[ithread] = cnnl;

  for (E_Int i = 0; i < ne; i++)
  {
    e = indir[i];
    s = cEEN[e].size();
    
    for (unsigned int f = 0; f < nfaces; f++)
    {
      for (E_Int v = 0; v < s; v++) // pour tous les voisins
      {
        ev = cEEN[e][v];
        for (unsigned int f2 = 0; f2 < nfaces; f2++)
        {
          // full check
          tag = 0;
          for (unsigned int n1 = 0; n1 < nof; n1++)
            for (unsigned int n2 = 0; n2 < nof; n2++)
              if (cn(e,face[f][n1]) == cn(ev,face[f2][n2])) tag += 1;
          if (tag == nof) { goto next; } 
        }
      } 

      // add face
      for (unsigned int n = 0; n < nof; n++)
      {
        (*cnnl)(nn,n+1) = cn(e,face[f][n]);
      }
      nn++;
      vectOfFacesL.push_back(e*nfaces+f+1);
      next: ;
    }
  }
  nes[ithread] = nn;
  delete [] indir;
  }
 
  // Nombre de faces total
  E_Int nntot = 0;
  for (E_Int i = 0; i < nthreads; i++) { prev[i] = nntot; nntot += nes[i]; }

  cnn.malloc(nntot, nof);

  // numero des faces
  PyObject* indicesFaces = NULL;//numpy of indices of exterior faces is created if boolIndir=true
  E_Int* indicesf = NULL;
  vector<E_Int> posFacesV(nthreads);
  E_Int sizeL = 0; 
  if (boolIndir == true)
  {
    indicesFaces = K_NUMPY::buildNumpyArray(nntot, 1, 1, 0);
    indicesf = K_NUMPY::getNumpyPtrI(indicesFaces);
    for (E_Int ithread = 0; ithread < nthreads; ithread++)
    {
      posFacesV[ithread] = sizeL;
      sizeL += vectOfFaces[ithread].size();
    }
  }
#pragma omp parallel default(shared)
  {
    E_Int  ithread = __CURRENT_THREAD__;
    E_Int ne = nes[ithread];
    E_Int p = prev[ithread];
    FldArrayI* cnnl = cnnls[ithread];
    for (size_t n = 1; n <= nof; n++)
    {
      E_Int* cnnp = cnn.begin(n);
      E_Int* cnp = cnnl->begin(n);
      for (E_Int e = 0; e < ne; e++) cnnp[e+p] = cnp[e];
    }
    if (boolIndir == true)
    {
      vector<E_Int>& vectOfFacesL = vectOfFaces[ithread];
      E_Int shiftIndices = posFacesV[ithread];
      for (size_t i = 0; i < vectOfFacesL.size(); i++)
      {
        indicesf[i+shiftIndices] = vectOfFacesL[i];
      }
    }
  }
  for (E_Int i = 0; i < nthreads; i++) delete cnnls[i];
  delete [] prev; delete [] nes; delete [] cnnls;
  if (boolIndir == true) 
  {
    PyList_Append(indices, indicesFaces);  
    Py_DECREF(indicesFaces);
  }
  return 1;
}
//=============================================================================
// IN: varString: varString de f
// IN: f: le champ
// IN: cn: connectivite NGon
// IN: outIndir: si true, sort un tableau d'indirection des faces en plus
// OUT: array des faces exterieures
// CAS NGON 3D : exterior faces are 2D zones
//==============================================================================
PyObject* K_POST::selectExteriorFacesNGon3D(char* varString, FldArrayF& f, 
                                            FldArrayI& cn, PyObject* indices)
{
  // CAS 3D
  bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;

  // cFE: connectivite face/elements.
  // Si une face n'a pas d'element gauche ou droit, retourne 0 pour 
  // cet element. 
  FldArrayI cFE; K_CONNECT::connectNG2FE(cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int* ptr = cn.begin();
  E_Int sizeFN = ptr[1]; ptr += 2; // sizeFN: taille du tableau de connectivite Face/Noeuds
  E_Int next1=0; // nbre de faces exterieures pour la nouvelle connectivite 
  E_Int next2=0; // nbre d'elements exterieurs pour la nouvelle connectivite
  E_Int sizeco1=0; // taille de la nouvelle connectivite de faces exterieures
  E_Int sizeco2=0; // taille de la nouvelle connectivite d'elements exterieurs
  E_Int e1, e2, nbnodes;

  FldArrayI posFaces;
  K_CONNECT::getPosFaces(cn, posFaces);
  E_Int* posFacep = posFaces.begin();

  // calcul du nombre d'elements et de faces pour la nouvelle connectivite
  E_Int fa = 0; E_Int numFace = 0; E_Int nfaceExt = 0;
  vector<E_Int> exteriorFaces;
  PyObject* indir = NULL;
  E_Int* indirp = NULL;

  while (fa < sizeFN) // parcours de la connectivite face/noeuds
  {
    e1 = cFE1[numFace];  // element voisin 1
    e2 = cFE2[numFace];  // element voisin 2
    nbnodes = ptr[0];

    if ((e1 == 0 && e2 != 0) || (e2 == 0 && e1 != 0))
    {
      sizeco1 += 3*nbnodes; next1 += nbnodes;
      sizeco2 += nbnodes+1; next2++;
      nfaceExt++;
      exteriorFaces.push_back(numFace);
    }
    ptr += nbnodes+1; fa += nbnodes+1; numFace++;
  }
  if (boolIndir == true)
  {
    indir = K_NUMPY::buildNumpyArray(nfaceExt, 1, 1, 0);
    indirp = K_NUMPY::getNumpyPtrI(indir);
  }
  cFE.malloc(0);
  // Calcul des nouvelles connectivites Elmt/Faces et Face/Noeuds
  FldArrayI c2n(sizeco1+sizeco2+4);
  E_Int* ptro1 = c2n.begin()+2;
  E_Int* ptro2 = c2n.begin()+sizeco1+4;

  E_Int npts = f.getSize(); E_Int nfld = f.getNfld();
  FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
  E_Int* indirNp = indirNodes.begin();
  E_Int indnewface = 1;
  E_Int indvertp1, indvertp2;
  E_Int indvertn = 0;
  vector< vector<E_Int> > indicesFaces(npts);// faces creees associees aux noeuds 
  E_Int pos, ind;
  E_Int n = exteriorFaces.size();

  for (E_Int i = 0; i < n; i++)
  {
    ind = exteriorFaces[i];
    pos = posFacep[ind];
    ptr = cn.begin()+pos;
    nbnodes = ptr[0];
    if (boolIndir == true) indirp[i] = ind+1;

    ptro2[0] = nbnodes;
    for (E_Int p = 0; p < nbnodes; p++)
    { 
      ptro1[0] = 2;
      indvertp1 = ptr[p+1]-1;
      indvertp2 = ptr[(p+1)%(nbnodes)+1]-1;
      vector<E_Int>& facesp1 = indicesFaces[indvertp1];
      vector<E_Int>& facesp2 = indicesFaces[indvertp2];

      //creation de la face P1P2 
      E_Int foundFace = -1; 
      if (facesp1.size() > 0 && facesp2.size() > 0 ) 
      {
        for (size_t fi = 0; fi < facesp1.size(); fi++)
          for (size_t fj = 0; fj < facesp2.size(); fj++)
          {
            if (facesp1[fi] == facesp2[fj] ) 
            {
              foundFace = facesp1[fi]; 
              break;
            }
            if (foundFace != -1) break;
          }        
      }
      if (foundFace == -1)
      { 
        if (indirNp[indvertp1] == -1) 
        {
          indvertn += 1;
          indirNp[indvertp1] = indvertn;
          ptro1[1] = indvertn;
        }
        else 
        {
          ptro1[1] = indirNp[indvertp1];
        }
        
        if (indirNp[indvertp2] == -1) 
        {
          indvertn += 1;
          indirNp[indvertp2] = indvertn;
          ptro1[2] = indvertn;
        }
        else 
        {
          ptro1[2] = indirNp[indvertp2];
        }
        ptro1 += 3;

        facesp1.push_back(indnewface); facesp2.push_back(indnewface);
        ptro2[p+1] = indnewface; indnewface++;
      }
      else 
      {
        ptro2[p+1] = foundFace;
      }
    }
    ptro2 += nbnodes+1;
  }
  FldArrayF* f2 = new FldArrayF(indvertn,nfld);
  E_Int indf;
  for(E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fp = f.begin(eq);
    E_Float* fnp = f2->begin(eq);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      indf = indirNp[ind]-1;
      if (indf > -1) fnp[indf] = fp[ind];
    }
  }
  // Cree la nouvelle connectivite complete
  E_Int nbFaces = indnewface-1;
  sizeFN = nbFaces*3;
  E_Int nbElts = next2;
  E_Int sizeEF = sizeco2;
  FldArrayI* cnout = new FldArrayI(sizeFN+sizeEF+4);
  E_Int* cnoutp = cnout->begin();
  ptro1 = c2n.begin()+2;
  ptro2 = c2n.begin()+4+sizeco1;
  cnoutp[0] = nbFaces; // nombre de faces
  cnoutp[1] = sizeFN; // taille du tableau de faces
  cnoutp += 2;
  for(E_Int i = 0; i < sizeFN; i++) cnoutp[i] = ptro1[i];
  cnoutp += sizeFN;
    
  cnoutp[0] = nbElts; 
  cnoutp[1] = sizeEF;
  cnoutp += 2;
  for (E_Int i = 0; i < sizeEF; i++) cnoutp[i] = ptro2[i];
  // Build array 
  char eltTypeFaces[10]; strcpy(eltTypeFaces, "NGON");
  PyObject* tpl = K_ARRAY::buildArray(*f2, varString, 
                                      *cnout, -1, eltTypeFaces);
  delete f2; delete cnout;
  if (boolIndir == true) 
  {
    PyList_Append(indices, indir); Py_DECREF(indir);
  }
  return tpl;
}

//=============================================================================
// IN: varString: varString de f
// IN: f: le champ
// IN: cn: connectivite NGon
// IN: outIndir: si true, sort un tableau d'indirection des faces en plus
// OUT: array des faces exterieures
// CAS NGON 2D : exterior faces are 1D zones
//==============================================================================
PyObject* K_POST::selectExteriorFacesNGon2D(char* varString, FldArrayF& f, 
                                            FldArrayI& cn, PyObject* indices)
{
  // CAS 2D
  bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;

  // cFE: connectivite face/elements.
  // Si une face n'a pas d'element gauche ou droit, retourne 0 pour 
  // cet element. 
  FldArrayI cFE; K_CONNECT::connectNG2FE(cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int* ptr = cn.begin();
  E_Int sizeFN = ptr[1]; ptr += 2; // sizeFN: taille du tableau de connectivite Face/Noeuds
  E_Int next1=0; // nbre de faces exterieures pour la nouvelle connectivite 
  E_Int next2=0; // nbre d'elements exterieurs pour la nouvelle connectivite
  E_Int sizeco1=0; // taille de la nouvelle connectivite de faces exterieures
  E_Int sizeco2=0; // taille de la nouvelle connectivite d'elements exterieurs
  E_Int e1, e2, nbnodes;

  FldArrayI posFaces;
  K_CONNECT::getPosFaces(cn, posFaces);
  E_Int* posFacep = posFaces.begin();

  // calcul du nombre d'elements et de faces pour la nouvelle connectivite
  E_Int fa = 0; E_Int numFace = 0; E_Int nfaceExt = 0;
  vector<E_Int> exteriorFaces;
  PyObject* indir = NULL;
  E_Int* indirp = NULL;

  while (fa < sizeFN) // parcours de la connectivite face/noeuds
  {
    e1 = cFE1[numFace];  // element voisin 1
    e2 = cFE2[numFace];  // element voisin 2
    nbnodes = ptr[0];

    if ( (e1 == 0 && e2 != 0) || (e2 == 0 && e1 != 0) )
    {
      sizeco1 += 2*nbnodes; next1 += nbnodes;
      sizeco2 += nbnodes+1; next2++;
      nfaceExt++;
      exteriorFaces.push_back(numFace);
    }
    ptr += nbnodes+1; fa += nbnodes+1; numFace++;
  }
  if (boolIndir == true)
  {
    indir = K_NUMPY::buildNumpyArray(nfaceExt, 1, 1, 0);
    indirp = K_NUMPY::getNumpyPtrI(indir);
  }
  cFE.malloc(0);
  // Calcul des nouvelles connectivites Elmt/Faces et Face/Noeuds
  FldArrayI c2n(sizeco1+sizeco2+4);
  E_Int* ptro1 = c2n.begin()+2;
  E_Int* ptro2 = c2n.begin()+sizeco1+4;

  E_Int npts = f.getSize(); E_Int nfld = f.getNfld();
  FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
  E_Int* indirNp = indirNodes.begin();
  E_Int indnewface = 1;
  E_Int indvertp1;
  E_Int indvertn = 0;
  FldArrayI indicesFaces(npts);// faces creees associees aux noeuds 
  indicesFaces.setAllValuesAt(-1);
  E_Int* indicesFacesp = indicesFaces.begin();
  E_Int pos, ind, foundFace;
  E_Int n = exteriorFaces.size();

  for (E_Int i = 0; i < n; i++)
  {
    ind = exteriorFaces[i];
    pos = posFacep[ind];
    ptr = cn.begin()+pos;
    nbnodes = ptr[0];
    if ( boolIndir == true) indirp[i] = ind+1;

    ptro2[0] = nbnodes;
    for (E_Int p = 0; p < nbnodes; p++)
    { 
      ptro1[0] = 1;
      indvertp1 = ptr[p+1]-1;
      foundFace = indicesFacesp[indvertp1];
      if ( foundFace == -1)
      { 
        if ( indirNp[indvertp1] == -1 ) 
        {
          indvertn += 1;
          indirNp[indvertp1] = indvertn;
          ptro1[1] = indvertn;
        }
        else 
        {
          ptro1[1] = indirNp[indvertp1];
        }             
        ptro1 += 2;
        indicesFacesp[indvertp1] = indnewface;
        ptro2[p+1] = indnewface; indnewface++;
      }
      else 
      {
        ptro2[p+1] = foundFace;
      }
    }
    ptro2 += nbnodes+1;
  }
  FldArrayF* f2 = new FldArrayF(indvertn,nfld);    
  for(E_Int eq = 1; eq <= nfld; eq++)
  {
    E_Float* fp = f.begin(eq);
    E_Float* fnp = f2->begin(eq);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      E_Int indf = indirNp[ind]-1;
      if ( indf> -1) fnp[indf] = fp[ind];
    }
  }
  // Cree la nouvelle connectivite complete
  E_Int nbFaces = indnewface-1;
  sizeFN = nbFaces*2;
  E_Int nbElts = next2;
  E_Int sizeEF = sizeco2;
  FldArrayI* cnout = new FldArrayI(sizeFN+sizeEF+4);
  E_Int* cnoutp = cnout->begin();
  ptro1 = c2n.begin()+2;
  ptro2 = c2n.begin()+4+sizeco1;
  cnoutp[0] = nbFaces; // nombre de faces
  cnoutp[1] = sizeFN; // taille du tableau de faces
  cnoutp+=2;
  for(E_Int i = 0; i < sizeFN; i++)
    cnoutp[i] = ptro1[i];

  cnoutp+=sizeFN;
    
  cnoutp[0] = nbElts; 
  cnoutp[1] = sizeEF;
  cnoutp+=2;
  for(E_Int i = 0; i < sizeEF; i++)
    cnoutp[i] = ptro2[i];
  // Build array 
  char eltTypeFaces[10]; strcpy(eltTypeFaces, "NGON");
  PyObject* tpl = K_ARRAY::buildArray(*f2, varString, 
                                      *cnout, -1, eltTypeFaces);
  delete f2; delete cnout;
  if ( boolIndir == true) 
  {
    PyList_Append(indices, indir);  Py_DECREF(indir);
  }
  return tpl;
 
}
//===========================================================================
/* Creation du tableau des sommets de chaque face pour un element donne 
   ex: TRI: 1 face=1 arete=2 sommets
       TETRA: 1 face=1 triangle=3 sommets
   IN: et: no de l'element
   IN: cn: connectivite elements->noeuds
   OUT: face: indices des noeuds dans une face (face 1, face 2, face 3,...)
*/
//===========================================================================
short K_POST::buildFaceInfo(E_Int et, FldArrayI& cn, FldArrayI& face)
{
  E_Int nfaces = face.getSize(); // nb de faces par elt
  E_Int nvertex = face.getNfld(); // nb de sommet par face
  E_Int nvtot = cn.getNfld(); // nb de vertex par element
  vector<E_Int> indv(nvtot);
  for (E_Int v = 0; v < nvtot; v++) indv[v] = cn(et, v+1);

  E_Int* facep = face.begin();
  if (nfaces == 3 && nvertex == 2) // TRI
  {
    facep[0] = indv[0]; facep[3] = indv[1];
    facep[1] = indv[1]; facep[4] = indv[2];
    facep[2] = indv[2]; facep[5] = indv[0];
  }
  else if (nfaces == 4 && nvertex == 2) // QUAD
  {
    facep[0] = indv[0]; facep[4] = indv[1]; 
    facep[1] = indv[1]; facep[5] = indv[2];
    facep[2] = indv[2]; facep[6] = indv[3];
    facep[3] = indv[3]; facep[7] = indv[0];
  }
  else if (nfaces == 4 && nvertex == 3) // TETRA
  {
    facep[0] = indv[0]; facep[4] = indv[1]; facep[8] = indv[3];//124 
    facep[1] = indv[1]; facep[5] = indv[2]; facep[9] = indv[3]; //234
    facep[2] = indv[2]; facep[6] = indv[0]; facep[10] = indv[3]; //314
    facep[3] = indv[0]; facep[7] = indv[1]; facep[11] = indv[2]; //123
  }
  else if (nfaces == 6 && nvertex == 4) // HEXA
  {
    //1234
    E_Int no = 0;
    facep[no] = indv[0];
    facep[no+nfaces] = indv[1];
    facep[no+2*nfaces] = indv[2]; 
    facep[no+3*nfaces] = indv[3]; 
    //5678
    no = 1;
    facep[no] = indv[4]; 
    facep[no+nfaces] = indv[5];
    facep[no+2*nfaces] = indv[6]; 
    facep[no+3*nfaces] = indv[7];  
    //1485
    no = 2;
    facep[no] = indv[0]; 
    facep[no+nfaces] = indv[3];
    facep[no+2*nfaces] = indv[7]; 
    facep[no+3*nfaces] = indv[4]; 
    //2376
    no = 3;
    facep[no] = indv[1]; 
    facep[no+nfaces] = indv[2];
    facep[no+2*nfaces] = indv[6]; 
    facep[no+3*nfaces] = indv[5];
    //1265
    no = 4;
    facep[no] = indv[0]; 
    facep[no+nfaces] = indv[1];
    facep[no+2*nfaces] = indv[5]; 
    facep[no+3*nfaces] = indv[4];
    //4378
    no = 5;
    facep[no] = indv[3]; 
    facep[no+nfaces] = indv[2];
    facep[no+2*nfaces] = indv[6]; 
    facep[no+3*nfaces] = indv[7];
  }
  else 
  {
    printf("Warning: buildFaceInfo: element type not implemented.");
    return 0;
  }
  return 1;
}

//===========================================================================
/* Creation du tableau des sommets de chaque face
   ex: TRI: 1 face=1 arete=2 sommets
       TETRA: 1 face=1 triangle=3 sommets
   IN: cn: connectivite elements->noeuds
   OUT: face: indices dans le tableau de connectivite
*/
//===========================================================================
short K_POST::buildFaceInfo2(FldArrayI& cn, FldArrayI& face)
{
  E_Int nfaces = face.getSize(); // nb de faces par elt
  E_Int nvertex = face.getNfld(); // nb de sommet par face
  E_Int* facep = face.begin();

  if (nfaces == 3 && nvertex == 2) // TRI
  {
    facep[0] = 1; facep[3] = 2;
    facep[1] = 2; facep[4] = 3;
    facep[2] = 3; facep[5] = 1;
  }
  else if (nfaces == 4 && nvertex == 2) // QUAD
  {
    facep[0] = 1; facep[4] = 2; 
    facep[1] = 2; facep[5] = 3;
    facep[2] = 3; facep[6] = 4;
    facep[3] = 4; facep[7] = 1;
  }
  else if (nfaces == 4 && nvertex == 3) // TETRA
  {
    facep[0] = 1; facep[4] = 2; facep[8] = 4;//124 
    facep[1] = 2; facep[5] = 3; facep[9] = 4; //234
    facep[2] = 3; facep[6] = 1; facep[10] = 4; //314
    facep[3] = 1; facep[7] = 2; facep[11] = 3; //123
  }
  else if (nfaces == 6 && nvertex == 4) // HEXA
  {
    //1234
    E_Int no = 0;
    facep[no] = 1;
    facep[no+nfaces] = 2;
    facep[no+2*nfaces] = 3; 
    facep[no+3*nfaces] = 4; 
    //5678
    no = 1;
    facep[no] = 5; 
    facep[no+nfaces] = 6;
    facep[no+2*nfaces] = 7; 
    facep[no+3*nfaces] = 8;
    //1485
    no = 2;
    facep[no] = 1; 
    facep[no+nfaces] = 4;
    facep[no+2*nfaces] = 8; 
    facep[no+3*nfaces] = 5;
    //2376
    no = 3;
    facep[no] = 2; 
    facep[no+nfaces] = 3;
    facep[no+2*nfaces] = 7; 
    facep[no+3*nfaces] = 6;
    //1265
    no = 4;
    facep[no] = 1; 
    facep[no+nfaces] = 2;
    facep[no+2*nfaces] = 6; 
    facep[no+3*nfaces] = 5;
    //4378
    no = 5;
    facep[no] = 4; 
    facep[no+nfaces] = 3;
    facep[no+2*nfaces] = 7; 
    facep[no+3*nfaces] = 8;
  }
  else 
  {
    printf("Warning: buildFaceInfo2: element type not implemented.");
    return 0;
  }
  return 1;
}

//===========================================================================
/* Retourne 1 si toutes les facettes de face1 sont deja traitees */
//===========================================================================
short K_POST::testCommonFaces(FldArrayI& face1, FldArrayI& face2,
                              FldArrayIS& tag1)
{
  E_Int nfaces = face1.getSize();
  E_Int nvertex = face1.getNfld();
  E_Int* face1p = face1.begin();
  E_Int* face2p = face2.begin();

  short* tag1p = tag1.begin();
  // parcours de toutes les facettes d'un elt1
  short cnt;
  short found = 1;
  E_Int t;
  for (E_Int f1 = 0; f1 < nfaces; f1++)
  {
    if (tag1p[f1] == 0)
    {
      // parcours de toutes les facettes des autres elts
      for (E_Int f2 = 0; f2 < nfaces; f2++)
      {
        // test des sommets un par un
        cnt = 0;
        for (E_Int v1 = 0; v1 < nvertex; v1++)
        {
          t = face1p[f1+v1*nfaces];
          for (E_Int v2 = 0; v2 < nvertex; v2++)
          {
            if (t == face2p[f2+v2*nfaces])
            {cnt++; break;}
          }
        }
        
        if (cnt == nvertex)
        {
          tag1p[f1] = 1; goto next;
        }
      }
      found = 0;
    }
    next: continue;
  }
  
  return found;
}

//=============================================================================
// exteriorFaces pour les array non structures conformes
// Cet algorithme est fonde sur un critere geometrique: on essaie de
// retrouver les centres des facettes enregistres dans un KdTree
// Enregistre 2 fois: facette interne
// Enregistre 1 fois: facette externe
//=============================================================================
short K_POST::exteriorFacesBasic2(E_Int nfaces, E_Int nvertex,
                                  FldArrayF& f, FldArrayI& cn,
                                  E_Int posx, E_Int posy, E_Int posz,
                                  FldArrayF& fnodes, FldArrayI& connect,
                                  PyObject* indices)
{
  E_Int ne = cn.getSize(); // nbre d'elements
  E_Int nfld = f.getNfld(); // nb de champs
  E_Float* fx = f.begin(posx);
  E_Float* fy = f.begin(posy);
  E_Float* fz = f.begin(posz);

  E_Int nfacesmax = ne*nfaces;
  FldArrayI face(nfaces, nvertex);
  K_POST::buildFaceInfo2(cn, face);

  FldArrayF interfaceCenters(nfacesmax, 3);
  E_Float* interfacex = interfaceCenters.begin(1);
  E_Float* interfacey = interfaceCenters.begin(2);
  E_Float* interfacez = interfaceCenters.begin(3);
  E_Float fvertex = 1./nvertex;

#pragma omp parallel default(shared)
  {
    E_Float xm, ym, zm;
    E_Int ind, ff;
#pragma omp for
    for (E_Int et = 0; et < ne; et++)
    {
      for (E_Int j = 0; j < nfaces; j++)
      {
        xm = 0.; ym = 0.; zm = 0.;
        for (E_Int k = 0; k < nvertex; k++) 
        {
          ind = cn(et, face[j+k*nfaces])-1;
          xm += fx[ind]; ym += fy[ind]; zm += fz[ind];
        }
        xm = xm * fvertex; ym = ym * fvertex; zm = zm * fvertex;
        ff = j+nfaces*et;
        interfacex[ff] = xm;
        interfacey[ff] = ym;
        interfacez[ff] = zm;
      }
    }
  }

  ArrayAccessor<FldArrayF> coordAcc(interfaceCenters, 1, 2, 3);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);
  
  E_Int ff = 0;
  E_Float pt[3];
  E_Float dx, dy, dz;
  E_Int ind, indp;
  short* tag = new short [nfacesmax];
  for (E_Int i = 0; i < nfacesmax; i++) tag[i] = 0;
  E_Int nfaceExt = 0;
  for (E_Int et = 0; et < ne; et++)
  {
    for (E_Int j = 0; j < nfaces; j++)
    {
      ff = et*nfaces+j;
      pt[0] = interfacex[ff];
      pt[1] = interfacey[ff];
      pt[2] = interfacez[ff]; 
      ind = globalKdt.getClosest(pt);
      indp = globalKdt.getClosest(ind);
      dx = interfacex[indp]-interfacex[ind];
      dy = interfacey[indp]-interfacey[ind];
      dz = interfacez[indp]-interfacez[ind];
      if (dx*dx + dy*dy + dz*dz > 1.e-24) // exterieur
      {
        tag[ff] = 1; nfaceExt++;
      }
    }
  }
  interfaceCenters.malloc(0);
  
  connect.malloc(nfaceExt, nvertex);
  fnodes.malloc(nfaceExt*nvertex, nfld);
  E_Int ee = 0; E_Int nn = 0;

  for (E_Int et = 0; et < ne; et++)
  {
    for (E_Int j = 0; j < nfaces; j++)
    {
      if (tag[et*nfaces+j] > 0.5) // exterieur
      {
        for (E_Int k = 0; k < nvertex; k++)
        {
          ind = cn(et, face[j+k*nfaces])-1;
          for (E_Int t = 1; t <= nfld; t++) fnodes(nn, t) = f(ind, t);
          connect(ee, k+1) = nn+1;
          nn++;
        }
        ee++;
      }
    }
  }
  delete [] tag;

  return 1;
}

//=============================================================================
// exteriorFaces pour les array non structures a elt basique
// Cet algorithme est fonde sur la topologie: une facette est externe
// si elle n'a qu'un seul voisin
//=============================================================================
short K_POST::exteriorFacesBasic(E_Int nfaces, E_Int nvertex,
                                 FldArrayF& f, FldArrayI& cn, 
                                 FldArrayF& fnodes, FldArrayI& connect)
{
  E_Int ne = cn.getSize(); // nbre d'elements
  E_Int nfld = f.getNfld(); // nb de champs
  
  E_Int nfacesmax = ne*nfaces;
  fnodes.malloc(nvertex*nfacesmax, nfld);
  connect.malloc(nfacesmax, nvertex);
  
  FldArrayI face1(nfaces, nvertex); // indices des sommets des faces
  FldArrayI face2(nfaces, nvertex);
  FldArrayIS tag1(nfaces);
  E_Int cf = 0; //compteur facettes
  E_Int cc = 0; //compteur connect
  short fin = 0;

  // Calcul des bounding box de chaque cellule
  FldArrayF bbox;
  K_COMPGEOM::boundingBoxOfUnstrCells(cn,
                                      f.begin(1), f.begin(2), f.begin(3),
                                      bbox);

  E_Float* bbox1 = bbox.begin(1);
  E_Float* bbox2 = bbox.begin(2);
  E_Float* bbox3 = bbox.begin(3);
  E_Float* bbox4 = bbox.begin(4);
  E_Float* bbox5 = bbox.begin(5);
  E_Float* bbox6 = bbox.begin(6);

  for (E_Int e1 = 0; e1 < ne; e1++)
  {  
    fin = 0;
    // construction du tableau face1
    short ok = buildFaceInfo(e1, cn, face1);
    if (ok == 0) return 0;
    tag1.setAllValuesAtNull();

    for (E_Int e2 = 0; e2 < ne; e2++)
    {
      if (e1 != e2 && 
          K_COMPGEOM::compBoundingBoxIntersection(
            bbox1[e1],bbox4[e1],bbox2[e1],bbox5[e1],bbox3[e1],bbox6[e1],
            bbox1[e2],bbox4[e2],bbox2[e2],bbox5[e2],bbox3[e2],bbox6[e2]) 
          == 1)
      {
        ok = buildFaceInfo(e2, cn, face2);
        if (ok == 0) return 0;
        // tag1 vaut 1 pour les faces communes
        fin = testCommonFaces(face1, face2, tag1);
        if (fin == 1) goto next;
      }
    }
    
    // recuperation des faces non interieures
    for (E_Int fi = 0; fi < nfaces; fi++)
    {
      if (tag1[fi] == 0) // pas commune a un autre elt
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fnodesp = fnodes.begin(n);
          E_Float* fp = f.begin(n);
        
          for (E_Int v1 = 0; v1 < nvertex; v1++)
          {
            E_Int s1 = face1(fi, v1+1)-1;
            fnodesp[cf+v1] = fp[s1];
          }
        }
        for (E_Int v1 = 1; v1 <= nvertex; v1++)
          connect(cc, v1) = cf + v1;

        cc += 1;
        cf += nvertex;
      }
    }
    next:;
  }
  fnodes.reAllocMat(cf, nfld);
  connect.reAllocMat(cc, nvertex);
  return 1;
}

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

// selectExteriorFaces

#include <stdio.h>
#include <string.h>
#include <map>
#include <unordered_map>
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
  if (!PYPARSETUPLE_(args, OO_, &array, &indices)) return NULL;
  
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
      E_Int nv;
      E_Int* ngon = cn->getNGon();
      E_Int* indPG = cn->getIndPG();
      cn->getFace(0, nv, ngon, indPG);
      if (nv == 2) tpl = selectExteriorFacesNGon2D(varString, *f, *cn, indices);
      else tpl = selectExteriorFacesNGon3D(varString, *f, *cn, indices);
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
  E_Int api = 1;
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
    tpl = K_ARRAY::buildArray3(nfld, varString, 2, 0, "NODE", 0, api);
    FldArrayF* fnodes; K_ARRAY::getFromArray3(tpl, fnodes);
    E_Int s = f.getSize();
    for (E_Int i = 1; i <= nfld; i++)
    {
      (*fnodes)(0,i) = f(0,i);
      (*fnodes)(1,i) = f(s-1,i);
    }
    
    if (boolIndir)
    {
      indir = K_NUMPY::buildNumpyArray(2, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);
      indirp[0] = 1; indirp[1] = s;
      PyList_Append(indices, indir); Py_DECREF(indir);
    }

    RELEASESHAREDS(tpl, fnodes);
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

    tpl = K_ARRAY::buildArray3(nfld, varString, nfacesExt, nfacesExt, newEltType, 0, api);
    FldArrayF* fnodes; FldArrayI* connect;
    K_ARRAY::getFromArray3(tpl, fnodes, connect);
    FldArrayI& cm = *(connect->getConnect(0));

    if (boolIndir)
    {
      indir = K_NUMPY::buildNumpyArray(nfacesExt, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);
    }

    #pragma omp parallel
    {
      E_Int e;
      #pragma omp for
      for (E_Int noe = 0; noe < nfacesExt-1; noe++)
      {
        cm(noe,1) = noe+1; cm(noe,2) = noe+2;
      }

      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        E_Float* fv = f.begin(nv);
        E_Float* fnv = fnodes->begin(nv);
        // Border j=1, i=1..p-1
        #pragma omp for
        for (E_Int i = 0; i < p-1; i++)
        { 
          fnv[i] = fv[i];
        }
        // Border i=p,j=0..q-1
        #pragma omp for
        for (E_Int j = 0; j < q-1; j++)
        { 
          e = p-1+j;
          fnv[e] = fv[p-1+j*p];
        }
        // Border j=jmax,i=p,0
        #pragma omp for
        for (E_Int i = p-1; i > 0; i--)
        { 
          e = nintj+q-1-i;
          fnv[e] = fv[i+(q-1)*p];
        }
        // Border i=1,j=q,1
        #pragma omp for
        for (E_Int j = q-1; j > 0; j--)
        { 
          e = nfacesExt-j;
          fnv[e] = fv[j*ni];
        }
      }

      if (boolIndir)
      {
        E_Int nbinti = p*(q-1);
        // j=1 interfaces
        #pragma omp for
        for (E_Int i = 0; i < p-1; i++)
        {
          indirp[i] = nbinti+i;
        }
        // i=imax interfaces
        #pragma omp for
        for (E_Int j = 0; j < q-1; j++)
        {
          e = p-1+j;
          indirp[e] = (p-1)+ j*p;
        }

        // j=jmax interfaces
        #pragma omp for
        for (E_Int i = 0; i < p-1; i++)
        { 
          e = p+q-2+i;
          indirp[e] = nbinti+(p-2-i)+(q-1)*(p-1);
        }

        // i=0 interface
        #pragma omp for
        for (E_Int j = 0; j < q-1; j++)
        { 
          e = nintj+q-1+j;
          indirp[e] = (q-2-j)*p;
        }
      }
    }

    cm(nfacesExt-1, 1) = nfacesExt; cm(nfacesExt-1, 2) = 1;
    if (boolIndir)
    {
      PyList_Append(indices, indir); Py_DECREF(indir);
    }
    RELEASESHAREDU(tpl, fnodes, connect);
    return tpl;
  }
  else
  {
    strcpy(newEltType, "QUAD");
    E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
    E_Int ni1nj1 = ni1*nj1;
    E_Int nPtsExt = 2*(ni*nj+ni*nk+nj*nk);
    E_Int nfacesExt = 2*(ni1nj1+nj1*nk1+ni1*nk1);

    tpl = K_ARRAY::buildArray3(nfld, varString, nPtsExt, nfacesExt, newEltType, 0, api);
    FldArrayF* fnodes; FldArrayI* connect;
    K_ARRAY::getFromArray3(tpl, fnodes, connect);
    FldArrayI& cm = *(connect->getConnect(0));

    if (boolIndir)
    {
      indir  = K_NUMPY::buildNumpyArray(nfacesExt, 1, 1, 0);
      indirp = K_NUMPY::getNumpyPtrI(indir);       
    }

    #pragma omp parallel
    {
      E_Int ninti  = ni*nj1*nk1;
      E_Int nintj  = ni1*nj*nk1; 
      E_Int nintij = ninti+nintj;
      E_Int e, ne, ind, indint, ns, n, noint;
      
      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        E_Float* fv = f.begin(nv);
        E_Float* fnv = fnodes->begin(nv);
        #pragma omp for
        for (n = 0; n < ni*nj; n++) fnv[n] = fv[n];
      }
    
      #pragma omp for
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          e = j*ni1 + i;
          ind = i+j*ni + 1;
          cm(e,1) = ind;
          cm(e,2) = ind+1;
          cm(e,3) = ind+1+ni;
          cm(e,4) = ind+ni;
        }

      ns = ni*nj;
      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        E_Float* fv = f.begin(nv);
        E_Float* fnv = fnodes->begin(nv);
        #pragma omp for
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni; i++)
          { 
            n = ns + j*ni + i;
            fnv[n] = fv[i+j*ni+nk1*ni*nj];
          }
      }
      
      #pragma omp for
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          e = ni1nj1 + j*ni1 + i;
          ind = ns + i+j*ni + 1;
          cm(e,1) = ind;
          cm(e,2) = ind+1;
          cm(e,3) = ind+1+ni;
          cm(e,4) = ind+ni;
        }

      if (boolIndir)
      {
        // k = 0 followed by k=nk
        #pragma omp for
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            noint = j*ni1 + i;
            indint = i+j*ni1+nintij;
            indirp[noint] = indint;
          }
        #pragma omp for
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            noint = ni1nj1 + j*ni1 + i;
            indint = i+j*ni1+nk1*ni1nj1+nintij;
            indirp[noint] = indint;
          }
      }

      ns = 2*ni*nj;
      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        E_Float* fv = f.begin(nv);
        E_Float* fnv = fnodes->begin(nv);
        #pragma omp for
        for (E_Int j = 0; j < nk; j++)
          for (E_Int i = 0; i < ni; i++)
          { 
            n = ns + j*ni + i;
            fnv[n] = fv[i+j*ni*nj];
          }
      }
      
      ne = 2*ni1nj1;
      for (E_Int j = 0; j < nk1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          e = ne + j*ni1 + i;
          ind = ns + i+j*ni + 1;
          cm(e,1) = ind;
          cm(e,2) = ind+1;
          cm(e,3) = ind+1+ni;
          cm(e,4) = ind+ni;
        }
      
      ns = 2*ni*nj + ni*nk;
      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        E_Float* fv = f.begin(nv);
        E_Float* fnv = fnodes->begin(nv);
        #pragma omp for
        for (E_Int j = 0; j < nk; j++)
          for (E_Int i = 0; i < ni; i++)
          {
            n = ns + j*ni + i;
            fnv[n] = fv[i+(nj-1)*ni+j*ni*nj];
          }
      }

      ne = 2*ni1nj1 + ni1*nk1;
      #pragma omp for
      for (E_Int j = 0; j < nk1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          e = ne + j*ni1 + i;
          ind = ns + i+j*ni + 1;
          cm(e,1) = ind;
          cm(e,2) = ind+1;
          cm(e,3) = ind+1+ni;
          cm(e,4) = ind+ni;
        } 

      if (boolIndir)
      {
        // j = 0 followed by j=nj
        ne = 2*ni1nj1;
        #pragma omp for
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            noint = ne + k*ni1 + i;
            indint = i+k*ni1*nj+ninti;
            indirp[noint] = indint;
          }
        ne = 2*ni1nj1 + ni1*nk1;
        #pragma omp for
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            noint = ne + k*ni1 + i;
            indint = i+nj1*ni1+k*ni1*nj+ninti;
            indirp[noint] = indint;
          }
      }

      ns = 2*(ni*nj + ni*nk);
      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        E_Float* fv = f.begin(nv);
        E_Float* fnv = fnodes->begin(nv);
        #pragma omp for
        for (E_Int j = 0; j < nk; j++)
          for (E_Int i = 0; i < nj; i++)
          { 
            n = ns + j*nj + i;
            fnv[n] = fv[i*ni+j*ni*nj];
          }
      }

      ne = 2*(ni1nj1 + ni1*nk1);
      #pragma omp for
      for (E_Int j = 0; j < nk1; j++)
        for (E_Int i = 0; i < nj1; i++)
        {
          e = ne + j*nj1 + i;
          ind = ns + i+j*nj + 1; 
          cm(e,1) = ind;
          cm(e,2) = ind+1;
          cm(e,3) = ind+1+nj;
          cm(e,4) = ind+nj;
        }

      ns = 2*(ni*nj + ni*nk) + nj*nk;
      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        E_Float* fv = f.begin(nv);
        E_Float* fnv = fnodes->begin(nv);
        #pragma omp for
        for (E_Int j = 0; j < nk; j++)
          for (E_Int i = 0; i < nj; i++)
          { 
            n = ns + j*nj + i;
            fnv[n] = fv[ni-1+i*ni+j*ni*nj];
          }
      }
      
      ne = 2*(ni1nj1 + ni1*nk1) + nj1*nk1;
      #pragma omp for
      for (E_Int j = 0; j < nk1; j++)
        for (E_Int i = 0; i < nj1; i++)
        {
          e = ne + j*nj1 + i;
          ind = ns + i+j*nj + 1; 
          cm(e,1) = ind;
          cm(e,2) = ind+1;
          cm(e,3) = ind+1+nj;
          cm(e,4) = ind+nj;
        }

      if (boolIndir)
      {
        // i = 0 followed by i=ni
        ne = 2*(ni1nj1 + ni1*nk1);
        #pragma omp for
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj1; j++)
          {
            noint = ne + k*nj1 + j;
            indint = j*ni+k*ni*nj1;
            indirp[noint] = indint;
          }
        ne = 2*(ni1nj1 + ni1*nk1) + nj1*nk1;
        #pragma omp for
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int j = 0; j < nj1; j++)
          {
            noint = ne + k*nj1 + j;
            indint = ni1+j*ni+k*ni*nj1;
            indirp[noint] = indint;
          }
      }
    }

    // Clean connectivity
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;
    if (posx != 0 && posy != 0 && posz != 0)
      K_CONNECT::cleanConnectivity(posx, posy, posz, 
                                   1.e-12, newEltType,
                                   *fnodes, *connect);

    PyObject* tpl2 = K_ARRAY::buildArray3(*fnodes, varString, *connect, newEltType);
    RELEASESHAREDU(tpl, fnodes, connect);
    if (boolIndir)
    {
      PyList_Append(indices, indir); Py_DECREF(indir);
    }
    return tpl2;
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
  if (strcmp(eltType, "TRI") == 0 || strcmp(eltType, "QUAD") == 0) strcpy(elttypeout, "BAR");
  else if (strcmp(eltType, "TETRA") == 0) strcpy(elttypeout, "TRI");
  else if (strcmp(eltType, "HEXA") == 0) strcpy(elttypeout, "QUAD");
  else if (strcmp(eltType, "PYRA") == 0) strcpy(elttypeout, "QUAD");
  else if (strcmp(eltType, "PENTA") == 0) strcpy(elttypeout, "QUAD");
  else if (strcmp(eltType, "BAR") == 0) strcpy(elttypeout, "NODE");
	
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
      if (nExtNodes > 0) {indirp[0] = ind0+1; indirp[1] = ind1+1;}
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
  size_t nvoisins = 0; // nbre de voisin si interieur  
  E_Int nfaces = 0; // nbre de faces par elements
  E_Int nof = 0; // nbre de noeuds par face

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
    nfaces = 5; nof = 4; // QUAD degen
    face[0][0] = 1; face[0][1] = 4; face[0][2] = 3; face[0][3] = 2;
    face[1][0] = 1; face[1][1] = 2; face[1][2] = 5; face[1][3] = 1;
    face[2][0] = 2; face[2][1] = 3; face[2][2] = 5; face[2][3] = 2;
    face[3][0] = 3; face[3][1] = 4; face[3][2] = 5; face[3][3] = 3;
    face[4][0] = 4; face[4][1] = 1; face[4][2] = 5; face[4][3] = 4;
  }
  else if (strcmp(eltType, "PENTA") == 0) 
  {
    nvoisins = 5;
    nfaces = 5; nof = 4; // QUAD degen
    face[0][0] = 1; face[0][1] = 2; face[0][2] = 5; face[0][3] = 4;
    face[1][0] = 2; face[1][1] = 3; face[1][2] = 6; face[1][3] = 5;
    face[2][0] = 3; face[2][1] = 1; face[2][2] = 4; face[2][3] = 6;
    face[3][0] = 1; face[3][1] = 3; face[3][2] = 2; face[3][3] = 1;
    face[4][0] = 4; face[4][1] = 5; face[4][2] = 6; face[4][3] = 4;
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
  E_Int tag;
  E_Int s; E_Int e, ev;
  E_Int nn = 0; // nbre de faces ext
  FldArrayI* cnnl = new FldArrayI(ne*nfaces,nof);
  cnnls[ithread] = cnnl;

  for (E_Int i = 0; i < ne; i++)
  {
    e = indir[i];
    s = cEEN[e].size();
    
    for (E_Int f = 0; f < nfaces; f++)
    {
      for (E_Int v = 0; v < s; v++) // pour tous les voisins
      {
        ev = cEEN[e][v];
        for (E_Int f2 = 0; f2 < nfaces; f2++)
        {
          // full check
          tag = 0;
          for (E_Int n1 = 0; n1 < nof; n1++)
            for (E_Int n2 = 0; n2 < nof; n2++)
              if (cn(e,face[f][n1]) == cn(ev,face[f2][n2])) tag += 1;
          if (tag == nof) { goto next; } 
        }
      } 

      // add face
      for (E_Int n = 0; n < nof; n++)
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
  PyObject* indicesFaces = NULL; //numpy of indices of exterior faces is created if boolIndir=true
  E_Int* indicesf = NULL;
  vector<E_Int> posFacesV(nthreads);
  E_Int sizeL = 0; 
  if (boolIndir)
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
    for (E_Int n = 1; n <= nof; n++)
    {
      E_Int* cnnp = cnn.begin(n);
      E_Int* cnp = cnnl->begin(n);
      for (E_Int e = 0; e < ne; e++) cnnp[e+p] = cnp[e];
    }
    if (boolIndir)
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
  if (boolIndir) 
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
// IN: boolIndir: si true, sort un tableau d'indirection des faces en plus
// OUT: array des faces exterieures
// CAS NGON 3D : exterior faces are 2D zones
//==============================================================================
PyObject* K_POST::selectExteriorFacesNGon3D(char* varString, FldArrayF& f, 
                                            FldArrayI& cn, PyObject* indices)
{
  // CAS 3D
  bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;

  // Acces non universel sur le ptrs
  E_Int* ngon = cn.getNGon();
  E_Int* indPG = cn.getIndPG();
  E_Int nfaces = cn.getNFaces();
  E_Int npts = f.getSize(), nfld = f.getNfld();
  E_Int shift = 1, api = f.getApi();
  if (api == 3) shift = 0;

  // cFE: connectivite face/elements.
  // Si une face n'a pas d'element gauche ou droit, retourne 0 pour 
  // cet element.
  FldArrayI cFE; K_CONNECT::connectNG2FE(cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int sizeFN2, sizeEF2 = 0;
  E_Int fidx, e1, e2, nbnodes, nptsExt;
  E_Int indvertn = 1, indedgen = 1;

  // Calcul du nombre de points uniques, aretes uniques, et faces dans la
  // nouvelle connectivite 2D
  E_Int nedgesExt = 0, nfacesExt = 0;
  vector<E_Int> indirVertices(npts, -1);
  vector<E_Int> exteriorFaces;

  // Les vertices sont mappes pour calculer leur table d'indirection sachant 
  // que leur nombre n'est pas connu a priori et qu'ils ne sont pas parcourus
  // dans l'ordre.
  std::unordered_map<E_Int, E_Int> vertexMap;
  // Les aretes sont hashees pour determiner le nombre unique d'aretes et
  // construire la nouvelle connectivite 2D
  vector<E_Int> edge(2);
  std::pair<E_Int, E_Bool> initEdge(-1, false);
  //TopologyOpt E; std::map<TopologyOpt, std::pair<E_Int, E_Bool> > edgeMap;
  //Topology E; std::unordered_map<Topology, std::pair<E_Int, E_Bool>, JenkinsHash<Topology> > edgeMap;
  TopologyOpt E; std::unordered_map<TopologyOpt, std::pair<E_Int, E_Bool>, JenkinsHash<TopologyOpt> > edgeMap;

  for (E_Int i = 0; i < nfaces; i++)
  {
    e1 = cFE1[i]; // element voisin 1
    e2 = cFE2[i]; // element voisin 2
    if ((e1 == 0 && e2 != 0) || (e2 == 0 && e1 != 0))
    {
      E_Int* face = cn.getFace(i, nbnodes, ngon, indPG);
      sizeEF2 += nbnodes+shift; nfacesExt++;
      exteriorFaces.push_back(i+1);
      
      for (E_Int p = 0; p < nbnodes; p++)
      {
        edge[0] = face[p]-1;
        edge[1] = face[(p+1)%nbnodes]-1;
        //E.set(edge); // version with Topology
        E.set(edge.data(), 2); // version with TopologyOpt
        // Ensure V and E are initially mapped to an initial value if either
        // doesn't exist
        auto resV = vertexMap.insert(std::make_pair(edge[0], 0));
        auto resE = edgeMap.insert(std::make_pair(E, initEdge));
        // Increment the value associated with V. If it is 1, then first
        // time this vertex is met, set indirVertices
        if (++resV.first->second == 1)
        {
          indirVertices[edge[0]] = indvertn;
          indvertn++;
        }
        // If the value associated with E is -1, then first time this edge
        // is encountered, set to current unique edge count
        if (resE.first->second.first == -1)
        {
          resE.first->second.first = indedgen;
          indedgen++;
        }
      }
    }
  }

  nptsExt = vertexMap.size();
  nedgesExt = edgeMap.size();
  sizeFN2 = (2+shift)*nedgesExt;
    
  PyObject* indir = NULL;
  E_Int* indirp = NULL;
  if (boolIndir)
  {
    indir = K_NUMPY::buildNumpyArray(nfacesExt, 1, 1, 0);
    indirp = K_NUMPY::getNumpyPtrI(indir);
  }
  cFE.malloc(0);
  
  // Calcul des nouvelles connectivites Elmt/Faces et Face/Noeuds
  E_Int ngonType = 1; // CGNSv3 compact array1
  if (api == 2) ngonType = 2; // CGNSv3, array2
  else if (api == 3) ngonType = 3; // force CGNSv4, array3
  E_Boolean center = false;
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, nptsExt, nfacesExt,
                                       nedgesExt, "NGON", sizeFN2, sizeEF2,
                                       ngonType, center, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);
  E_Int* ngon2 = cn2->getNGon();
  E_Int* nface2 = cn2->getNFace();
  E_Int *indPG2 = NULL, *indPH2 = NULL;
  
  if (api == 2 || api == 3)
  {
    indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
  }

  E_Int c1 = 0, c2 = 0; // positions in ngon2 and nface2
  for (E_Int i = 0; i < nfacesExt; i++) 
  {
    fidx = exteriorFaces[i];
    E_Int* face = cn.getFace(fidx-1, nbnodes, ngon, indPG);
    if (boolIndir) indirp[i] = fidx;
    
    nface2[c2] = nbnodes;
    if (api == 2 || api == 3) indPH2[i] = nbnodes;
    
    for (E_Int p = 0; p < nbnodes; p++)
    {
      edge[0] = face[p]-1; 
      edge[1] = face[(p+1)%nbnodes]-1;
    
      //E.set(edge); // version with Topology
      E.set(edge.data(), 2); // version with TopologyOpt

      // Find the edge in the edge map
      auto resE = edgeMap.find(E);
      //if (resE == edgeMap.end()) { printf("not found\n"); fflush(stdout); }
      
      if (not resE->second.second)
      {
        resE->second.second = true;
        ngon2[c1] = 2;
        ngon2[c1+shift] = indirVertices[edge[0]];
        ngon2[c1+1+shift] = indirVertices[edge[1]];
        c1 += 2+shift;
      }
      nface2[c2+p+shift] = resE->second.first;
    }
    c2 += nbnodes+shift;
  }
  
  #pragma omp parallel
  {
    E_Int indf;
    if (api == 2 || api == 3)
    {
      #pragma omp for
      for(E_Int i = 0; i < nedgesExt; i++) indPG2[i] = 2;
    }
  
    for(E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fp = f.begin(eq);
      E_Float* f2p = f2->begin(eq);
      #pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        indf = indirVertices[ind]-1;
        if (indf > -1) f2p[indf] = fp[ind];
      }
    }
  }
  
  if (boolIndir)
  {
    PyList_Append(indices, indir); Py_DECREF(indir);
  }
  RELEASESHAREDU(tpl, f2, cn2);
  return tpl;
}

//=============================================================================
// IN: varString: varString de f
// IN: f: le champ
// IN: cn: connectivite NGon
// IN: boolIndir: si true, sort un tableau d'indirection des faces en plus
// OUT: array des faces exterieures
// CAS NGON 2D : exterior faces are 1D zones
//==============================================================================
PyObject* K_POST::selectExteriorFacesNGon2D(char* varString, FldArrayF& f, 
                                            FldArrayI& cn, PyObject* indices)
{
  // CAS 2D
  bool boolIndir = false;
  if (indices != Py_None) boolIndir = true;

  // Acces non universel sur le ptrs
  E_Int* ngon = cn.getNGon();
  E_Int* indPG = cn.getIndPG();
  E_Int nfaces = cn.getNFaces();
  E_Int npts = f.getSize(), nfld = f.getNfld();
  E_Int shift = 1, api = f.getApi();
  if (api == 3) shift = 0;

  // cFE: connectivite face/elements.
  // Si une face n'a pas d'element gauche ou droit, retourne 0 pour 
  // cet element. 
  FldArrayI cFE; K_CONNECT::connectNG2FE(cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int sizeFN2, sizeEF2 = 0;
  E_Int e1, e2, dummy, nptsExt;
  E_Int indvertn = 1;

  // Calcul du nombre de points uniques et aretes uniques dans la
  // nouvelle connectivite 1D
  vector<E_Int> indirVertices(npts, -1);

  // Les vertices sont mappes pour calculer leur table d'indirection sachant 
  // que leur nombre n'est pas connu a priori et qu'ils ne sont pas parcourus
  // dans l'ordre.
  std::unordered_map<E_Int, E_Int > vertexMap;
  // Les aretes sont hashees pour determiner le nombre unique d'aretes et
  // ainsi construire les "elements" de la nouvelle connectivite 1D
  E_Int nedgesExt = 0;
  vector<E_Int> edge(2);
  vector<E_Int> exteriorEdges;
  TopologyOpt E;
  std::unordered_map<TopologyOpt, E_Int, JenkinsHash<TopologyOpt> > edgeMap;
  //Topology E;
  //std::unordered_map<Topology, E_Int, JenkinsHash<Topology> > edgeMap;


  for (E_Int i = 0; i < nfaces; i++)
  {
    e1 = cFE1[i]; // element voisin 1
    e2 = cFE2[i]; // element voisin 2
    E_Int* face = cn.getFace(i, dummy, ngon, indPG);
    if ((e1 == 0 && e2 != 0) || (e2 == 0 && e1 != 0))
    {
      // Increment the value associated with V/E. If it is 1, then first
      // time this vertex/edge is encountered
      edge[0] = face[0]-1;
      auto resV = vertexMap.insert(std::make_pair(edge[0], 0));
      if (++resV.first->second == 1)
      {
        indirVertices[edge[0]] = indvertn;
        indvertn++;
      }

      edge[1] = face[1]-1;
      resV = vertexMap.insert(std::make_pair(edge[1], 0));
      if (++resV.first->second == 1)
      {
        indirVertices[edge[1]] = indvertn;
        indvertn++;
      }

      E.set(edge.data(), 2);
      //E.set(edge);
      
      auto resE = edgeMap.insert(std::make_pair(E, 0));
      if (++resE.first->second == 1)
      {
        exteriorEdges.push_back(i+1);
      }
    }
  }

  nptsExt = vertexMap.size();
  sizeFN2 = (1+shift)*nptsExt;
  nedgesExt = edgeMap.size();
  sizeEF2 = (2+shift)*nedgesExt;

  PyObject* indir = NULL;
  E_Int* indirp = NULL;
  if (boolIndir)
  {
    indir = K_NUMPY::buildNumpyArray(nedgesExt, 1, 1, 0);
    indirp = K_NUMPY::getNumpyPtrI(indir);
  }
  cFE.malloc(0);
  // Calcul des nouvelles connectivites Elmt/Faces et Face/Noeuds
  E_Int ngonType = 1; // CGNSv3 compact array1
  if (api == 2) ngonType = 2; // CGNSv3, array2
  else if (api == 3) ngonType = 3; // force CGNSv4, array3
  E_Boolean center = false;
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, nptsExt, nedgesExt,
                                       nptsExt, "NGON", sizeFN2, sizeEF2,
                                       ngonType, center, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);
  E_Int* ngon2 = cn2->getNGon();
  E_Int* nface2 = cn2->getNFace();
  E_Int *indPG2 = NULL, *indPH2 = NULL;
  if (api == 2 || api == 3)
  {
    indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
  }
  
  #pragma omp parallel
  {
    E_Int ind, fidx, v1, v2;

    #pragma omp for
    for(E_Int i = 0; i < nptsExt; i++)
    {
      ind = i*(1+shift);
      ngon2[ind] = 1;
      ngon2[ind+shift] = i+1;
    }
    #pragma omp for
    for(E_Int i = 0; i < nedgesExt; i++)
    {
      fidx = exteriorEdges[i];
      E_Int* face = cn.getFace(fidx-1, dummy, ngon, indPG);
      if (boolIndir) indirp[i] = fidx;
      v1 = face[0]-1; v2 = face[1]-1;

      ind = i*(2+shift);
      nface2[ind] = 2;
      nface2[ind+shift] = indirVertices[v1];
      nface2[ind+1+shift] = indirVertices[v2];
    }

    if (api == 2 || api == 3)
    {
      #pragma omp for
      for(E_Int i = 0; i < nptsExt; i++) indPG2[i] = 1;
      #pragma omp for
      for(E_Int i = 0; i < nedgesExt; i++) indPH2[i] = 2;
    }
  
    for(E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fp = f.begin(eq);
      E_Float* f2p = f2->begin(eq);
      #pragma omp for
      for (E_Int i = 0; i < npts; i++)
      {
        ind = indirVertices[i]-1;
        if (ind > -1) f2p[ind] = fp[i];
      }
    }
  }

  RELEASESHAREDU(tpl, f2, cn2);
  if (boolIndir) 
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
  short found = 1;
 
  #pragma omp parallel
  {
    E_Int t;
    short cnt;
    // parcours de toutes les facettes d'un elt1
    #pragma omp for
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
  E_Int nfacesExt = 0;
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
        tag[ff] = 1; nfacesExt++;
      }
    }
  }
  interfaceCenters.malloc(0);
  
  connect.malloc(nfacesExt, nvertex);
  fnodes.malloc(nfacesExt*nvertex, nfld);
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

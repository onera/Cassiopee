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

# include <stdio.h>
# include <stdlib.h>
# include "converter.h"

using namespace std;
using namespace K_FLD;

// 1: array1, 2: array1/2, 3: array1/2/3
#define ARRAYCODE 3

//=============================================================================
/* Copy the contain of an array in another array */
//=============================================================================
PyObject* K_CONVERTER::copy(PyObject* self, PyObject* args)
{
  PyObject* array;  
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

#if ARRAYCODE == 1
/* array1 only code */
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, eltType, true);

  if (res == 1)
  { 
    PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, 
                                        ni, nj, nk);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    K_KCORE::memcpy__(fnp, f->begin(), f->getSize()*f->getNfld());
    //memcpy(fnp, f->begin(), f->getSize()*f->getNfld()*sizeof(E_Float));
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    E_Int csize = cn->getSize()*cn->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString,
                                        f->getSize(), cn->getSize(), 
                                        -1, eltType, false, csize);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    //memcpy(fnp, f->begin(), f->getSize()*f->getNfld()*sizeof(E_Float));
    K_KCORE::memcpy__(fnp, f->begin(), f->getSize()*f->getNfld());
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    return NULL;
  }
#endif

#if ARRAYCODE == 2
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray2(array, varString, f, ni, nj, nk, cn, eltType);
  E_Int api = f->getApi();
  E_Int nfld = f->getNfld(); E_Int npts = f->getSize();

  if (res == 1)
  { 
    // for all arrays
    PyObject* tpl = K_ARRAY::buildArray2(nfld, varString, ni, nj, nk, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray2(tpl, varString, f2, ni, nj, nk, cn2, eltType);

    #pragma omp parallel
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f->begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for
        for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      }
    } 
    RELEASESHAREDS(array, f); RELEASESHAREDS(tpl, f2);
    return tpl;
  }
  else if (res == 2)
  {
    E_Boolean center = false;
    if (strchr(eltType, '*') != NULL) center = true;
    E_Int nfld = f->getNfld(); E_Int npts = f->getSize();
    E_Int nelts=0, nfaces=0, sizeNGon=0, sizeNFace=0;
    if (strcmp(eltType, "NGON") == 0 || strcmp(eltType, "NGON*") == 0)
    {
      nelts = cn->getNElts();
      nfaces = cn->getNFaces();
      sizeNGon = cn->getSizeNGon();
      sizeNFace = cn->getSizeNFace();
      //printf("sizengon=%d, sizenface=%d\n", sizeNGon, sizeNFace);
    }
    else
    {
      nelts = cn->getSize();
    }
    //printf("eltType %s, centers=%d\n", eltType, center); 
    PyObject* tpl = K_ARRAY::buildArray2(
      nfld, varString, npts, nelts,-1, 
      eltType, center, sizeNGon, sizeNFace, nfaces, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray2(tpl, varString, f2, ni, nj, nk, cn2, eltType);

    // copie des champs
    #pragma omp parallel
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f->begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for
        for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      }
    }
    // copie de la connectivite
    if (sizeNGon > 0 && cn->isNGon() == 2) // NGon pour Array2
    { 
      E_Int* ngon = cn->getNGon();
      E_Int* nface = cn->getNFace();
      E_Int* indPG = cn->getIndPG();
      E_Int* indPH = cn->getIndPH();
      E_Int* ngon2 = cn2->getNGon();
      E_Int* nface2 = cn2->getNFace();
      E_Int* indPG2 = cn2->getIndPG(); // exists??
      E_Int* indPH2 = cn2->getIndPH();
      
      #pragma omp parallel
      {
        #pragma omp for
        for (E_Int i = 0; i < sizeNGon; i++) ngon2[i] = ngon[i];
        #pragma omp for
        for (E_Int i = 0; i < sizeNFace; i++) nface2[i] = nface[i];
        if (indPG != NULL)
        {
          #pragma omp for
          for (E_Int i = 0; i < nfaces; i++) indPG2[i] = indPG[i];
        }
        if (indPH != NULL)
        {
          #pragma omp for
          for (E_Int i = 0; i < nelts; i++) indPH2[i] = indPH[i];
        }
      }
      if (indPG == NULL)
      {
          // Compute from NGON - not threaded !!
          E_Int pos = 0;
          for (E_Int i = 0; i < nfaces; i++) { indPG2[i] = pos; pos += ngon[pos]+1; }
      }  
      if (indPH == NULL)
      {
          // Compute from NFACE - not threaded !!
          E_Int pos = 0;
          for (E_Int i = 0; i < nelts; i++) { indPH2[i] = pos; pos += nface[pos]+1; }
      }
    }
    else if (sizeNGon > 0) // NGon pour Array1
    {
      E_Int* cnp = cn->begin();
      E_Int* cn2p = cn2->begin();
      #pragma omp parallel
      {
        #pragma omp for
        for (E_Int i = 0; i < sizeNGon+sizeNFace+4; i++) cn2p[i] = cnp[i]; 
      }
    }
    else // elements basiques tout Array
    {
      E_Int* cnp = cn->begin();
      E_Int* cn2p = cn2->begin();
      E_Int nppe = cn->getNfld();
      #pragma omp parallel
      {
        #pragma omp for
        for (E_Int i = 0; i < nelts*nppe; i++) cn2p[i] = cnp[i];
      }
    } 
    RELEASESHAREDU(array, f, cn); RELEASESHAREDU(tpl, f2, cn2);
    return tpl;
  }
  else
  {
    return NULL;
  }
#endif

#if ARRAYCODE == 3
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  E_Int api = f->getApi();
  if (api == 2) api = 3;
  E_Int nfld = f->getNfld(); E_Int npts = f->getSize();

  if (res == 1)
  { 
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, ni, nj, nk, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);

    #pragma omp parallel
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f->begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for
        for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
      }
    } 
    RELEASESHAREDS(array, f); RELEASESHAREDS(tpl, f2);
    return tpl;
  }
  else if (res == 2)
  {
    PyObject* tpl = NULL;
    if (strcmp(eltType, "NGON") == 0 || strcmp(eltType, "NGON*") == 0)
    {
      E_Int dim; 
      E_Int ngonType = cn->isNGon();
      E_Int nelts = cn->getNElts();
      E_Int nfaces = cn->getNFaces();
      E_Int sizeNGon = cn->getSizeNGon();
      E_Int sizeNFace = cn->getSizeNFace();
      E_Boolean center = false;
      if (strchr(eltType, '*') != NULL) center = true;
      PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts, nfaces, 
        eltType, sizeNGon, sizeNFace, ngonType, center, api);
      FldArrayF* f2; FldArrayI* cn2;
      K_ARRAY::getFromArray3(tpl, f2, cn2);

      if (center == true) dim = nelts;
      else dim = npts;
      #pragma omp parallel
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fp = f->begin(n);
          E_Float* f2p = f2->begin(n);
          #pragma omp for
          for (E_Int i = 0; i < dim; i++) f2p[i] = fp[i];
        }
        E_Int* ngonp = cn->getNGon();
        E_Int* ngon2p = cn2->getNGon();
        #pragma omp for
        for (E_Int i = 0; i < sizeNGon; i++) ngon2p[i] = ngonp[i];
        E_Int* nfacep = cn->getNFace();
        E_Int* nface2p = cn2->getNFace();
        #pragma omp for
        for (E_Int i = 0; i < sizeNFace; i++) nface2p[i] = nfacep[i];
        if (api > 1)
        {
          E_Int* indPGp = cn->getIndPG();
          E_Int* indPG2p = cn2->getIndPG();
          //if (indPGp == NULL) { printf("indPG is null\n"); fflush(stdout); }
          if (ngonType == 2) dim = nfaces;
          else dim = nfaces+1;
          #pragma omp for
          for (E_Int i = 0; i < dim; i++) indPG2p[i] = indPGp[i];
          E_Int* indPHp = cn->getIndPH();
          E_Int* indPH2p = cn2->getIndPH();
          if (ngonType == 2) dim = nelts;
          else dim = nelts+1;
          #pragma omp for
          for (E_Int i = 0; i < dim; i++) indPH2p[i] = indPHp[i];
        } 
      }
      RELEASESHAREDU(array, f, cn); RELEASESHAREDU(tpl, f2, cn2);
      return tpl;
    }
    else // BE
    {
      E_Int ncon = cn->getNConnect();
      E_Boolean center = false;
      if (strchr(eltType, '*') != NULL) center = true;
      vector< E_Int > neltsPerConnect(ncon);
      E_Int nelts = 0;
      for (E_Int i = 0; i < ncon; i++)
      { FldArrayI& cm = *(cn->getConnect(i));
        neltsPerConnect[i] = cm.getSize(); 
        nelts += cm.getSize(); }
      tpl = K_ARRAY::buildArray3(nfld, varString, npts, neltsPerConnect, eltType, center, api);
      FldArrayF* f2; FldArrayI* cn2;  
      K_ARRAY::getFromArray3(tpl, f2, cn2);
      
      // copie des champs
      E_Int dim;
      if (center == true) dim = nelts;
      else dim = npts;
      #pragma omp parallel
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fp = f->begin(n);
          E_Float* f2p = f2->begin(n);
          #pragma omp for
          for (E_Int i = 0; i < dim; i++) f2p[i] = fp[i];
        }
      
        for (E_Int i = 0; i < ncon; i++)
        { 
          FldArrayI& cm = *(cn->getConnect(i));
          FldArrayI& cm2 = *(cn2->getConnect(i));
          E_Int* cmp = cm.begin();
          E_Int* cm2p = cm2.begin();
          #pragma omp for
          for (E_Int i = 0; i < cm.getSize()*cm.getNfld(); i++) cm2p[i] = cmp[i];
        }
      }
      RELEASESHAREDU(array, f, cn); RELEASESHAREDU(tpl, f2, cn2);
      return tpl;
    }
  }
  else
  {
    return NULL;
  }
#endif
}

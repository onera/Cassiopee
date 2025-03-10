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
# include <string.h>
# include "post.h"

using namespace K_FLD;
using namespace std;

extern "C"
{
    void k6compstructmetric_(
    const E_Int& im, const E_Int& jm, const E_Int& km,
    const E_Int& nbcells, const E_Int& nintt,
    const E_Int& ninti, const E_Int& nintj,
    const E_Int& nintk,
    E_Float* x, E_Float* y, E_Float* z,
    E_Float* vol, E_Float* surfx, E_Float* surfy, E_Float* surfz,
    E_Float* snorm, E_Float* cix, E_Float* ciy, E_Float* ciz);
}

//=============================================================================
/* Calcul du gradient d'un ensemble de champs en centres
   Le gradient est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeGrad2NGon(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayc;
  PyObject* volc; PyObject* cellNc;
  PyObject* indices; PyObject* field;
  if (!PyArg_ParseTuple(args, "OOOOOO", &array, &arrayc, &volc, &cellNc,
                        &indices, &field)) return NULL;

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn,
                                     eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeGrad2: invalid array.");
    return NULL;
  }
  if (res == 1 || strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDB(res,array,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "computeGrad2: only for NGons.");
    return NULL;
  }

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeGrad2: coordinates not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posx++; posy++; posz++;

  // Check arrayc
  char* varStringc; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int nic, njc, nkc; // number of points of array
  res = K_ARRAY::getFromArray3(arrayc, varStringc, fc, nic, njc, nkc, cnc,
                               eltTypec);

  // Extract cellN if any
  E_Float* cellNp = NULL;
  E_Int ncells = fc->getSize();
  if (cellNc != Py_None) K_NUMPY::getFromNumpyArray(cellNc, cellNp, ncells, true);

  // Nombre de variable dont il faut calculer le gradient
  E_Int nfld = fc->getNfld();
  vector<char*> vars;
  K_ARRAY::extractVars(varStringc, vars);

  vector<char*> varStrings;
  for (E_Int i = 0; i < nfld; i++)
  {
    char* local;
    computeGradVarsString(vars[i], local);
    varStrings.push_back(local);
  }
  
  E_Int size=0;
  for (E_Int i = 0; i < nfld; i++)
  {
    size += strlen(varStrings[i])+1;
  }
  char* varStringOut = new char [size];
  char* pt = varStringOut;
  for (E_Int i = 0; i < nfld; i++)
  {
    char* v = varStrings[i];
    for (size_t j = 0; j < strlen(v); j++)
    {
      *pt = v[j]; pt++;
    }
    *pt = ','; pt++;
    delete [] varStrings[i];
  }
  pt--; *pt = '\0';
  for (size_t i = 0; i < vars.size(); i++) delete [] vars[i];

  // Calcul FE
  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);

  // Calcul du champ aux faces
  E_Int nfaces = cn->getNFaces();
  E_Int nelts = cn->getNElts();
  
  FldArrayF faceField(nfaces, nfld);
  if (cellNp == NULL)
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = faceField.begin(n);
      E_Float* s = fc->begin(n);
      #pragma omp parallel
      {
        E_Int i1, i2;
        #pragma omp for
        for (E_Int i = 0; i < nfaces; i++)
        {
          i1 = cFE1[i] - 1; i2 = cFE2[i] - 1;
          if (i2 != -1) fp[i] = 0.5*(s[i1] + s[i2]);
          else fp[i] = s[i1];
        }
      }
    }
  }
  else // cellN
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = faceField.begin(n);
      E_Float* s = fc->begin(n);
      #pragma omp parallel
      {
        E_Int i1, i2;
        #pragma omp for
        for (E_Int i = 0; i < nfaces; i++)
        {
          i1 = cFE1[i] - 1; i2 = cFE2[i] - 1;
          if (i2 != -1)
          {
            if (cellNp[i1] == 0) fp[i] = s[i2];
            else fp[i] = 0.5*(s[i1] + s[i2]);
          }
          else fp[i] = s[i1];
        }
      }
    }
  }

  // Replace DataSet
  FldArrayI* inds=NULL; FldArrayF* bfield=NULL;
  if (indices != Py_None && field != Py_None)
  {
    K_NUMPY::getFromNumpyArray(indices, inds, true);
    K_NUMPY::getFromNumpyArray(field, bfield, true);

    E_Int n = inds->getSize()*inds->getNfld();
    E_Int* pind = inds->begin();
    E_Float* pf = bfield->begin();
    E_Float* fp = faceField.begin(1);
    
    #pragma omp parallel
    {
      E_Int ind;
      #pragma omp for
      for (E_Int i = 0; i < n; i++)
      {
        ind = pind[i]-1;
        fp[ind] = pf[i];
      }
    }
    
    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(field, bfield);
  }

  // Build unstructured NGON array from existing connectivity & empty fields
  FldArrayF* gp = new FldArrayF(nelts, nfld*3, true); gp->setAllValuesAtNull();
  PyObject* tpl = K_ARRAY::buildArray3(*gp, varStringOut, *cn, "NGON");
  delete gp; K_ARRAY::getFromArray3(tpl, gp);

  FldArrayF surf(nfaces, 4);
  E_Float* sxp = surf.begin(1);
  E_Float* syp = surf.begin(2);
  E_Float* szp = surf.begin(3);
  E_Float* snp = surf.begin(4);
  K_METRIC::compNGonFacesSurf(f->begin(posx), f->begin(posy),
                              f->begin(posz), *cn,
                              sxp, syp, szp, snp, &cFE);
  // gradient
  E_Int i1, i2;
  E_Float ff;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpx = gp->begin(3*n+1);
    E_Float* gpy = gp->begin(3*n+2);
    E_Float* gpz = gp->begin(3*n+3);
    E_Float* fp = faceField.begin(n+1);
    for (E_Int i = 0; i < nfaces; i++)
    {
      i1 = cFE1[i] - 1; i2 = cFE2[i] - 1;
      ff = fp[i];
      if (i1 != -1)
      {
        gpx[i1] += ff*sxp[i];
        gpy[i1] += ff*syp[i];
        gpz[i1] += ff*szp[i];
      }
      if (i2 != -1)
      {
        gpx[i2] -= ff*sxp[i];
        gpy[i2] -= ff*syp[i];
        gpz[i2] -= ff*szp[i];
      }
    }
  }
  
  // free mem
  surf.malloc(0); faceField.malloc(0);
  
  FldArrayF vol(nelts);
  E_Float* volp = vol.begin(1);
  if (volc == Py_None)
  {  
    K_METRIC::CompNGonVol(f->begin(posx), f->begin(posy),
  			  f->begin(posz), *cn, volp);
  }
  else
  {
    FldArrayF* vols=NULL; 
    K_NUMPY::getFromNumpyArray(volc, vols, true);
    volp = vols->begin();
    RELEASESHAREDN(volc, vols);
  }

  #pragma omp parallel
  {
    E_Float voli;
    for (E_Int n = 0; n < nfld; n++)
    {
      E_Float* gpx = gp->begin(3*n+1);
      E_Float* gpy = gp->begin(3*n+2);
      E_Float* gpz = gp->begin(3*n+3);
      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        voli = 1./K_FUNC::E_max(volp[i], K_CONST::E_MIN_VOL);
        gpx[i] = gpx[i] * voli;
        gpy[i] = gpy[i] * voli;
        gpz[i] = gpz[i] * voli;
      }
    }
  }

  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDU(arrayc, fc, cnc);
  RELEASESHAREDS(tpl, gp);
  if (cellNc != Py_None) Py_DECREF(cellNc);

  delete [] varStringOut;
  return tpl;
}

//=============================================================================
/* Calcul du gradient d'un ensemble de champs en centres
   Le gradient est fourni aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeGrad2Struct(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayc;
  PyObject* indices; PyObject* field; PyObject* cellNc;
  if (!PyArg_ParseTuple(args, "OOOOO", &array, &arrayc, &cellNc,
                        &indices, &field)) return NULL;

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn,
                                     eltType);
  if (res != 1)
  {
    if (res == 2) RELEASESHAREDU(array,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "computeGrad2: invalid array.");
    return NULL;
  }
  E_Int dimPb = 3;
  if (nk == 1)
  {
    if (nj > 1) dimPb = 2;
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeGrad2: not valid for 1D structured arrays.");
      RELEASESHAREDS(array,f); return NULL;
    }
  }
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeGrad2: coordinates not found in array.");
    RELEASESHAREDS(array,f); return NULL;
  }
  posx++; posy++; posz++;

  // Check arrayc
  char* varStringc; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int nic, njc, nkc; // number of points of array
  res = K_ARRAY::getFromArray3(arrayc, varStringc, fc, nic, njc, nkc, cnc,
                               eltTypec);


  // Extract cellN if any
  E_Float* cellNp = NULL;
  E_Int ncells = fc->getSize();
  if (cellNc != Py_None) K_NUMPY::getFromNumpyArray(cellNc, cellNp, ncells, true);

  // Nombre de variables dont il faut calculer le gradient
  E_Int nfld = fc->getNfld();
  vector<char*> vars;
  K_ARRAY::extractVars(varStringc, vars);

  vector<char*> varStrings;
  for (E_Int i = 0; i < nfld; i++)
  {
    char* local;
    computeGradVarsString(vars[i], local);
    varStrings.push_back(local);
  }
  E_Int size=0;
  for (E_Int i = 0; i < nfld; i++)
  {
    size += strlen(varStrings[i])+1;
  }
  char* varStringOut = new char [size];
  char* pt = varStringOut;
  for (E_Int i = 0; i < nfld; i++)
  {
    char* v = varStrings[i];
    for (size_t j = 0; j < strlen(v); j++)
    {
      *pt = v[j]; pt++;
    }
    *pt = ','; pt++;
    delete [] varStrings[i];
  }
  pt--; *pt = '\0';
  for (size_t i = 0; i < vars.size(); i++) delete [] vars[i];

  E_Int nicnjc = nic*njc;
  E_Int ninjc = ni*njc;
  E_Int nicnj = nic*nj;
  E_Int nbIntI = ninjc*nkc;
  E_Int nbIntJ = nicnj*nkc;
  E_Int nbIntK = nicnjc*nk;
  E_Int nbIntIJ = nbIntI+nbIntJ;
  E_Int nbIntTot = nbIntIJ+nbIntK;
  
  // Build
  FldArrayF faceField(nbIntTot,nfld); faceField.setAllValuesAtNull();
  FldArrayI voisins(nbIntTot,2); voisins.setAllValuesAt(-1);
  E_Int* cellG = voisins.begin(1);
  E_Int* cellD = voisins.begin(2);
  PyObject* tpl;
  if (dimPb == 2)
    tpl = computeGrad2Struct2D(ni, nj, nic, njc, varStringOut, cellNp, 
                               f->begin(posx), f->begin(posy), f->begin(posz),
                               *fc, faceField, cellG, cellD, indices, field);
  else if (dimPb == 3)
    tpl = computeGrad2Struct3D(ni, nj, nk, nic, njc, nkc, varStringOut, cellNp, 
                               f->begin(posx), f->begin(posy), f->begin(posz),
                               *fc, faceField, cellG, cellD, indices, field);
  else
  { 
    PyErr_SetString(PyExc_TypeError,
                    "computeGrad2: invalid pb dimension.");
    tpl = NULL;
  }
  delete [] varStringOut;
  RELEASESHAREDS(array,f);
  RELEASESHAREDS(arrayc, fc);
  if (cellNc != Py_None) Py_DECREF(cellNc);
  return tpl;
}

//=============================================================================
PyObject* K_POST::computeGrad2Struct3D(E_Int ni, E_Int nj, E_Int nk,
                                       E_Int nic, E_Int njc, E_Int nkc,
                                       const char* varStringOut, E_Float* cellNp,
                                       E_Float* xt, E_Float* yt, E_Float* zt,
                                       FldArrayF& fc, FldArrayF& faceField,
                                       E_Int* cellG, E_Int* cellD,
                                       PyObject* indices, PyObject* field)
{
  E_Int nicnjc = nic*njc;
  E_Int ninjc = ni*njc;
  E_Int nicnj = nic*nj;
  E_Int nic1 = nic - 1; E_Int njc1 = njc - 1; E_Int nkc1 = nkc - 1;
  
  E_Int nbIntI = ninjc*nkc;
  E_Int nbIntJ = nicnj*nkc;
  E_Int nbIntK = nicnjc*nk;
  E_Int nbIntIJ = nbIntI+nbIntJ;
  E_Int nbIntTot = nbIntIJ+nbIntK;
  
  E_Int nfld = fc.getNfld(); E_Int nfldg = nfld*3;
  E_Int ncells = nicnjc*nkc;

  if (cellNp == NULL)
  {
    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fcn = fc.begin(eq);
      E_Float* fintp = faceField.begin(eq);
      
      #pragma omp parallel
      {
        E_Int i, j, k;
        E_Int indint, indcellg, indcelld;
        
        // faces en i
        #pragma omp for
        for (E_Int idx = 0; idx < nkc*njc*nic1; idx++) 
        {
          i = idx%nic1 + 1;
          j = (idx/nic1)%njc;
          k = idx/(nic1*njc);

          indint = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          indcelld = indcellg + 1;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en i
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc*njc; idx++) 
        {
          j = idx%njc;
          k = idx/njc;

          i = 0;
          indint = i + j*ni + k*ninjc;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          i = nic;
          indint = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en j
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc*njc1*nic; idx++) 
        {
          i = idx%nic;
          j = (idx/nic)%njc1 + 1;
          k = idx/(nic*njc1);

          indint = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          indcelld = indcellg + nic;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en j
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc*nic; idx++) 
        {
          i = idx%nic;
          k = idx/nic;

          j = 0;
          indint = i + j*nic + k*nicnj + nbIntI;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indint] = -1;
          cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          j = njc;
          indint = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en k
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc1*njc*nic; idx++) 
        {
          i = idx%nic;
          j = (idx/nic)%njc;
          k = idx/(nic*njc) + 1;

          indint = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          indcelld = indcellg + nicnjc;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }

        // bords des faces en k
        #pragma omp for nowait
        for (E_Int idx = 0; idx < njc*nic; idx++) 
        {
          i = idx%nic;
          j = idx/nic;

          k = 0;
          indint = i + j*nic + k*nicnjc + nbIntIJ;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld];

          k = nkc;
          indint = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg];
        }
      }
    }
  }
  else // cellN
  {
    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fcn = fc.begin(eq);
      E_Float* fintp = faceField.begin(eq);
      
      #pragma omp parallel
      {
        E_Int i, j, k;
        E_Int indint, indcellg, indcelld;
        
        // faces en i
        #pragma omp for
        for (E_Int idx = 0; idx < nkc*njc*nic1; idx++) 
        {
          i = idx%nic1 + 1;
          j = (idx/nic1)%njc;
          k = idx/(nic1*njc);

          indint = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          indcelld = indcellg + 1;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indint] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indint] = fcn[indcellg];
          else fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en i
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc*njc; idx++) 
        {
          j = idx%njc;
          k = idx/njc;

          i = 0;
          indint = i + j*ni + k*ninjc;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          i = nic;
          indint = i + j*ni + k*ninjc;
          indcellg = (i - 1) + j*nic + k*nicnjc;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en j
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc*njc1*nic; idx++) 
        {
          i = idx%nic;
          j = (idx/nic)%njc1 + 1;
          k = idx/(nic*njc1);

          indint = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          indcelld = indcellg + nic;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indint] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indint] = fcn[indcellg];
          else fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en j
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc*nic; idx++) 
        {
          i = idx%nic;
          k = idx/nic;

          j = 0;
          indint = i + j*nic + k*nicnj + nbIntI;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indint] = -1;
          cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          j = njc;
          indint = i + j*nic + k*nicnj + nbIntI;
          indcellg = i + (j - 1)*nic + k*nicnjc;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }
        
        // faces en k
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nkc1*njc*nic; idx++) 
        {
          i = idx%nic;
          j = (idx/nic)%njc;
          k = idx/(nic*njc) + 1;

          indint = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          indcelld = indcellg + nicnjc;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indint] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indint] = fcn[indcellg];
          else fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en k
        #pragma omp for nowait
        for (E_Int idx = 0; idx < njc*nic; idx++) 
        {
          i = idx%nic;
          j = idx/nic;

          k = 0;
          indint = i + j*nic + k*nicnjc + nbIntIJ;
          indcelld = i + j*nic + k*nicnjc;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld];

          k = nkc;
          indint = i + j*nic + k*nicnjc + nbIntIJ;
          indcellg = i + j*nic + (k - 1)*nicnjc;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg];
        }
      }
    }
  }

  // Replace DataSet
  if (indices != Py_None && field != Py_None)
  {
    FldArrayI* inds=NULL; FldArrayF* bfield=NULL;
    K_NUMPY::getFromNumpyArray(indices, inds, true);
    K_NUMPY::getFromNumpyArray(field, bfield, true);

    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pindint = inds->begin();
    E_Float* pf = bfield->begin();
    
    #pragma omp parallel
    {
      E_Int indint;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fintp = faceField.begin(eq);
        #pragma omp for
        for (E_Int noint = 0; noint < ninterfaces; noint++)
        {
          indint = pindint[noint];
          fintp[indint] = pf[noint];
        }
      }
    }
    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(field, bfield);
  }

  PyObject* tpl = K_ARRAY::buildArray3(nfldg,varStringOut, nic, njc, nkc);
  E_Float* gnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF gp(ncells, nfldg, gnp, true); gp.setAllValuesAtNull();
  FldArrayF surf(nbIntTot,3);
  FldArrayF centerInt(nbIntTot,3);
  E_Float* sxp = surf.begin(1);
  E_Float* syp = surf.begin(2);
  E_Float* szp = surf.begin(3);
  FldArrayF surfnorm(nbIntTot);
  E_Float* snp = surfnorm.begin();
  FldArrayF vol(ncells); E_Float* volp = vol.begin();

  k6compstructmetric_(ni, nj, nk, ncells, nbIntTot, nbIntI, nbIntJ, nbIntK,
                      xt, yt, zt,
                      volp, sxp, syp, szp, snp,
                      centerInt.begin(1), centerInt.begin(2),centerInt.begin(3));
  
  // free mem
  centerInt.malloc(0); surfnorm.malloc(0);

  // gradient
  E_Int indcellg, indcelld;
  E_Float ff;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpx = gp.begin(3*n+1);
    E_Float* gpy = gp.begin(3*n+2);
    E_Float* gpz = gp.begin(3*n+3);
    E_Float* fp = faceField.begin(n+1);

    for (E_Int indint = 0; indint < nbIntTot; indint++)
    {
      indcellg = cellG[indint]; indcelld = cellD[indint];
      ff = fp[indint];
      if (indcellg != -1)
      {
        gpx[indcellg] += ff*sxp[indint];
        gpy[indcellg] += ff*syp[indint];
        gpz[indcellg] += ff*szp[indint];
      }
      if (indcelld != -1)
      {
        gpx[indcelld] -= ff*sxp[indint];
        gpy[indcelld] -= ff*syp[indint];
        gpz[indcelld] -= ff*szp[indint];
      }
    }
  }
  
  // free mem
  surf.malloc(0);
    
  #pragma omp parallel
  {
    E_Float voli;
    
    for (E_Int n = 0; n < nfld; n++)
    {
      E_Float* gpx = gp.begin(3*n+1);
      E_Float* gpy = gp.begin(3*n+2);
      E_Float* gpz = gp.begin(3*n+3);
      
      #pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        voli = 1./K_FUNC::E_max(volp[i], K_CONST::E_MIN_VOL);
        gpx[i] *= voli; gpy[i] *= voli; gpz[i] *= voli;
      }
    }
  }
  
  return tpl;
}
//=============================================================================
PyObject* K_POST::computeGrad2Struct2D(E_Int ni, E_Int nj, E_Int nic, E_Int njc,
                                       const char* varStringOut, E_Float* cellNp,
                                       E_Float* xt, E_Float* yt, E_Float* zt,
                                       FldArrayF& fc, FldArrayF& faceField,
                                       E_Int* cellG, E_Int* cellD,
                                       PyObject* indices, PyObject* field)
{
  E_Int nkc = 1;
  E_Int nicnjc = nic*njc;
  E_Int ninjc = ni*njc;
  E_Int nicnj = nic*nj;
  E_Int nic1 = nic - 1; E_Int njc1 = njc - 1;
  
  E_Int nbIntI = ninjc;
  E_Int nbIntJ = nicnj;
  E_Int nbIntIJ = nbIntI+nbIntJ;
  E_Int nbIntTot = nbIntIJ;
  
  E_Int nfld = fc.getNfld(); E_Int nfldg = nfld*3;
  E_Int ncells = nicnjc;
  FldArrayF sint(nbIntTot,3); sint.setAllValuesAtNull();
  E_Float* sxint = sint.begin(1);
  E_Float* syint = sint.begin(2);
  E_Float* szint = sint.begin(3);
  
  // Compute length of faces
  #pragma omp parallel
  {
    E_Int i, j, indint;
    E_Int indm, indp;
    E_Float d13x, d13y, d13z, d24x, d24y, d24z;
    
    #pragma omp for
    for (E_Int idx = 0; idx < ninjc; idx++)
    {
      i = idx%ni;
      j = idx/ni;
      
      indm = i + j*ni; indp = indm + ni;
      d13x = xt[indp] - xt[indm];
      d13y = yt[indp] - yt[indm];
      d13z = 1;
      d24x = xt[indm] - xt[indp];
      d24y = yt[indm] - yt[indp];
      d24z = 1;

      sxint[idx] = 0.5*(d13y*d24z - d13z*d24y);
      syint[idx] = 0.5*(d13z*d24x - d13x*d24z);
      szint[idx] = 0.5*(d13x*d24y - d13y*d24x);
    }
    
    #pragma omp for nowait
    for (E_Int idx = 0; idx < nicnj; idx++)
    {
      i = idx%nic;
      j = idx/nic;
      indint = ninjc + idx;
      
      indm = i + j*ni; indp = indm + 1;
      d13x = xt[indp] - xt[indm];
      d13y = yt[indp] - yt[indm];
      d13z = 1;
      d24x = xt[indp] - xt[indm];
      d24y = yt[indp] - yt[indm];
      d24z = -1;

      sxint[indint] = 0.5*(d13y*d24z - d13z*d24y);
      syint[indint] = 0.5*(d13z*d24x - d13x*d24z);
      szint[indint] = 0.5*(d13x*d24y - d13y*d24x);
    }
  }

  if (cellNp == NULL)
  {
    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fcn = fc.begin(eq);
      E_Float* fintp = faceField.begin(eq);
        
      #pragma omp parallel
      {
        E_Int i, j;
        E_Int indint, indcellg, indcelld;
        
        // faces en i internes
        #pragma omp for nowait
        for (E_Int idx = 0; idx < njc*nic1; idx++) 
        {
          i = idx%nic1 + 1;
          j = idx/nic1;
          
          indint = i + j*ni;
          indcellg = (i - 1) + j*nic; indcelld = indcellg + 1;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en i
        #pragma omp for nowait
        for (E_Int j = 0; j < njc; j++)
        {
          // faces i = 0
          indint = j*ni;
          indcelld = j*nic;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          // faces i = ni
          indint = (ni - 1) + j*ni;
          indcellg = (nic - 1) + j*nic;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }

        // faces en j internes
        #pragma omp for
        for (E_Int idx = 0; idx < nic*njc1; idx++) 
        {
          i = idx%nic;
          j = idx/nic + 1;
    
          indint = i + j*nic + nbIntI;
          indcellg = i + (j - 1)*nic; indcelld = indcellg + nic;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
        }
        
        // bords des faces en j
        #pragma omp for
        for (E_Int i = 0; i < nic; i++)
        {
          // faces j = 0
          indint = i + nbIntI;
          indcelld = i;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          // faces j = jmax
          indint = i + njc*nic + nbIntI;
          indcellg = i + (njc - 1)*nic;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }
      }
    }
  }
  else // cellN
  {
    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fcn = fc.begin(eq);
      E_Float* fintp = faceField.begin(eq);
        
      #pragma omp parallel
      {
        E_Int i, j;
        E_Int indint, indcellg, indcelld;
        
        // faces en i internes
        #pragma omp for
        for (E_Int idx = 0; idx < njc*nic1; idx++) 
        {
          i = idx%nic1 + 1;
          j = idx/nic1;
          
          indint = i + j*ni;
          indcellg = (i - 1) + j*nic; indcelld = indcellg + 1;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indint] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indint] = fcn[indcellg];
          else fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en i
        #pragma omp for nowait
        for (E_Int j = 0; j < njc; j++)
        {
          // faces i = 0
          indint = j*ni;
          indcelld = j*nic;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          // faces i = ni
          indint = (ni - 1) + j*ni;
          indcellg = (nic - 1) + j*nic;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }

        // faces en j internes
        #pragma omp for nowait
        for (E_Int idx = 0; idx < nic*njc1; idx++) 
        {
          i = idx%nic;
          j = idx/nic + 1;
    
          indint = i + j*nic + nbIntI;
          indcellg = i + (j - 1)*nic; indcelld = indcellg + nic;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          if (cellNp[indcellg] == 0) fintp[indint] = fcn[indcelld];
          else if (cellNp[indcelld] == 0) fintp[indint] = fcn[indcellg];
          else fintp[indint] = 0.5*(fcn[indcellg] + fcn[indcelld]);
        }
        
        // bords des faces en j
        #pragma omp for nowait
        for (E_Int i = 0; i < nic; i++)
        {
          // faces j = 0
          indint = i + nbIntI;
          indcelld = i;
          cellG[indint] = -1; cellD[indint] = indcelld;
          fintp[indint] = fcn[indcelld]; // Extrapolation de l interieur

          // faces j = jmax
          indint = i + njc*nic + nbIntI;
          indcellg = i + (njc - 1)*nic;
          cellG[indint] = indcellg; cellD[indint] = -1;
          fintp[indint] = fcn[indcellg]; // Extrapolation de l interieur
        }
      }
    }
  }
  
  // Replace DataSet
  if (indices != Py_None && field != Py_None)
  {
    FldArrayI* inds=NULL; FldArrayF* bfield=NULL;
    K_NUMPY::getFromNumpyArray(indices, inds, true);
    K_NUMPY::getFromNumpyArray(field, bfield, true);

    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pindint = inds->begin();
    E_Float* pf = bfield->begin();
    
    #pragma omp parallel
    {
      E_Int indint;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fintp = faceField.begin(eq);
        #pragma omp for
        for (E_Int noint = 0; noint < ninterfaces; noint++)
        {
          indint = pindint[noint];
          fintp[indint] = pf[noint];
        }
      }
    }
    
    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(field, bfield);
  }
  PyObject* tpl = K_ARRAY::buildArray3(nfldg,varStringOut, nic, njc, nkc);
  E_Float* gnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF gp(ncells, nfldg, gnp, true); gp.setAllValuesAtNull();

  // gradient
  E_Float ff;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpx = gp.begin(3*n+1);
    E_Float* gpy = gp.begin(3*n+2);
    E_Float* gpz = gp.begin(3*n+3);
    E_Float* fp = faceField.begin(n+1);
    for (E_Int indint = 0; indint < nbIntTot; indint++)
    {
      E_Int indcellg = cellG[indint]; E_Int indcelld = cellD[indint];
      ff = fp[indint];
      if (indcellg != -1)
      {
        gpx[indcellg] += ff*sxint[indint];
        gpy[indcellg] += ff*syint[indint];
        gpz[indcellg] += ff*szint[indint];
      }
      if (indcelld != -1)
      {
        gpx[indcelld] -= ff*sxint[indint];
        gpy[indcelld] -= ff*syint[indint];
        gpz[indcelld] -= ff*szint[indint];
      }
    }
  }
  
  // free mem
  sint.malloc(0);
    
  #pragma omp parallel
  {
    E_Float voli;
    for (E_Int n = 0; n < nfld; n++)
    {
      E_Float* gpx = gp.begin(3*n+1);
      E_Float* gpy = gp.begin(3*n+2);
      E_Float* gpz = gp.begin(3*n+3);
      #pragma omp for
      for (E_Int indcell = 0; indcell < ncells; indcell++)
      {
        voli = K_METRIC::compVolOfStructCell2D(ni, nj, xt, yt, zt, indcell, -1);
        voli = 1./K_FUNC::E_max(voli, K_CONST::E_MIN_VOL);
        gpx[indcell] *= voli; gpy[indcell] *= voli; gpz[indcell] *= voli;
      }
    }
  }
  
  return tpl;
}

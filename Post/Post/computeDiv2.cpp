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
/* Compute the divergence of a set of vector fields given in cell centers
   The divergence is given on cell centers. */
//=============================================================================
PyObject* K_POST::computeDiv2NGon(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayc;
  PyObject* indices; PyObject* fieldX; PyObject* fieldY; PyObject* fieldZ;
  if (!PyArg_ParseTuple(args, "OOOOOO", &array, &arrayc,
                        &indices, &fieldX, &fieldY, &fieldZ)) return NULL;

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn,
                                    eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: invalid array.");
    return NULL;
  }
  if (res == 1 || strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDB(res,array,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: only for NGons.");
    return NULL;
  }
  E_Int dimPb = 3; // Consider 3-dimensional vector fields even for 2D NGon

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: coordinates not found in array.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  }
  posx++; posy++; posz++;

  // Check arrayc
  char* varStringc; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int nic, njc, nkc; // number of points of array
  res = K_ARRAY::getFromArray(arrayc, varStringc, fc, nic, njc, nkc, cnc,
                              eltTypec, true);

  E_Int npts = f->getSize();

  // Number of vector fields whose divergence to compute (three components for each)
  E_Int nfld = fc->getNfld(); // total number of scalar fields
  vector<char*> vars;
  K_ARRAY::extractVars(varStringc, vars);
  if (nfld % dimPb != 0)
  {
    RELEASESHAREDB(res,array,f,cn);
    RELEASESHAREDB(res,arrayc,fc,cnc);
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: not all components were found for each vector field.");
    return NULL;
  }
  else nfld /= dimPb;

  vector<char*> varStrings;
  for (E_Int i = 0; i < nfld; i++)
  {
    for (E_Int m = 0; m < dimPb-1; m++)
    {
      if (strncmp(vars[dimPb*i+m], vars[dimPb*i+m+1], strlen(vars[dimPb*i+m])-1) != 0)
      {
        RELEASESHAREDB(res,array,f,cn);
        RELEASESHAREDB(res,arrayc,fc,cnc);
        PyErr_SetString(PyExc_TypeError,
                        "computeDiv2: invalid names for vector component fields.");
        return NULL;
      }
    }
    char* sv0 = vars[dimPb*i]; char* sv1 = vars[dimPb*i+1]; char* sv2 = vars[dimPb*i+2];
    char s0 = sv0[strlen(sv0)-1];
    char s1 = sv1[strlen(sv1)-1];
    char s2 = sv2[strlen(sv2)-1];
    if (s0 != 'X' || s1 != 'Y' || s2 != 'Z')
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeDiv2: error with the order of given scalar fields.");
      return NULL;
    }
    char* local;
    computeDivVarsString(vars[i*dimPb], local);
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

  // Compute FE connectivity
  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);

  // Compute field on element faces
  E_Int* cnp = cn->begin();
  E_Int nfaces = cnp[0];
  E_Int nelts = cnp[cnp[1]+2];
  // printf("nfaces=%d\n", nfaces);
  //printf("nelts=%d\n", nelts);

  FldArrayF faceField(nfaces, dimPb*nfld);
  E_Int i1, i2;
  for (E_Int n = 1; n <= dimPb*nfld; n++)
  {
    E_Float* fp = faceField.begin(n);
    E_Float* s = fc->begin(n);
    for (E_Int i = 0; i < nfaces; i++)
    {
      i1 = cFE1[i]-1; i2 = cFE2[i]-1;
      if (i1 != -1 && i2 != -1) fp[i] = 0.5*(s[i1]+s[i2]);
      else if (i1 == -1) fp[i] = s[i2];
      else if (i2 == -1) fp[i] = s[i1];
      else fp[i] = 0.;
    }
  }

  // Replace DataSet
  FldArrayI* inds=NULL; FldArrayF* bfieldX=NULL;
  FldArrayF* bfieldY=NULL; FldArrayF* bfieldZ=NULL;
  if (indices != Py_None && fieldX != Py_None && fieldY != Py_None
                                              && fieldZ != Py_None )
  {
    K_NUMPY::getFromNumpyArray(indices, inds, true);
    K_NUMPY::getFromNumpyArray(fieldX, bfieldX, true);
    K_NUMPY::getFromNumpyArray(fieldY, bfieldY, true);
    K_NUMPY::getFromNumpyArray(fieldZ, bfieldZ, true);

    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pind = inds->begin();

    E_Float* pf[3];
    pf[0] = bfieldX->begin();
    pf[1] = bfieldY->begin();
    pf[2] = bfieldZ->begin();

    for (E_Int n = 1; n <= dimPb*nfld; n++)
    {
      // printf("n : %d \n",n);
      E_Float* fp = faceField.begin(n);
      E_Int ind;

      for (E_Int i = 0; i < ninterfaces; i++)
      {
        ind = pind[i]-1;
        fp[ind] = pf[n-1][i];
      }
    }
  }
  // Build empty unstructured array
  PyObject* tpl = K_ARRAY::buildArray(nfld, varStringOut, npts,
                                      nelts, -1, eltType, true,
                                      cn->getSize()*cn->getNfld());

  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
  E_Float* gnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF gp(nelts, nfld, gnp, true); gp.setAllValuesAtNull();

  FldArrayF surf(nfaces, 4);
  E_Float* sxp = surf.begin(1);
  E_Float* syp = surf.begin(2);
  E_Float* szp = surf.begin(3);
  E_Float* snp = surf.begin(4);
  K_METRIC::compNGonFacesSurf(f->begin(posx), f->begin(posy),
                              f->begin(posz), *cn,
                              sxp, syp, szp, snp, &cFE);

  // divergence
  E_Float ffx, ffy, ffz;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpdv = gp.begin(n+1);
    E_Float* fpx = faceField.begin(3*n+1);
    E_Float* fpy = faceField.begin(3*n+2);
    E_Float* fpz = faceField.begin(3*n+3);
    for (E_Int i = 0; i < nfaces; i++)
    {
      i1 = cFE1[i]-1; i2 = cFE2[i]-1;
      ffx = fpx[i]; ffy = fpy[i]; ffz = fpz[i];
      //printf("sx=%f %f %f\n", sxp[i],syp[i],szp[i]);
      if (i1 != -1)
      {
        gpdv[i1] += ffx*sxp[i] + ffy*syp[i] + ffz*szp[i];
      }
      if (i2 != -1)
      {
        gpdv[i2] -= ffx*sxp[i] + ffy*syp[i] + ffz*szp[i];
      }
    }
  }
  // free mem
  surf.malloc(0); faceField.malloc(0);

  FldArrayF vol(nelts);
  E_Float* volp = vol.begin(1);
  K_METRIC::CompNGonVol(f->begin(posx), f->begin(posy),
                        f->begin(posz), *cn, volp);

  E_Float voli;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpdv = gp.begin(n+1);
    for (E_Int i = 0; i < nelts; i++)
    {
      //printf("vol=%f\n", volp[i]);
      voli = 1./K_FUNC::E_max(volp[i], 1.e-12);
      gpdv[i] *= voli;
    }
  }

  if (indices != Py_None && fieldX != Py_None )
  {
    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(fieldX, bfieldX);
    RELEASESHAREDN(fieldY, bfieldY);
    RELEASESHAREDN(fieldZ, bfieldZ);
  }

  RELEASESHAREDB(res,array,f,cn);
  RELEASESHAREDB(res,arrayc,fc,cnc);
  delete [] varStringOut;
  return tpl;
}

//=============================================================================
/* Compute the divergence of a set of vector fields given in cell centers
   The divergence is given on cell centers. */
//=============================================================================
PyObject* K_POST::computeDiv2Struct(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayc;
  PyObject* indices; PyObject* fieldX; PyObject* fieldY; PyObject* fieldZ;
  if (!PyArg_ParseTuple(args, "OOOOOO", &array, &arrayc,
                        &indices, &fieldX, &fieldY, &fieldZ)) return NULL;

  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn,
                                    eltType, true);
  if (res != 1)
  {
    if ( res == 2) RELEASESHAREDB(res,array,f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: invalid array.");
    return NULL;
  }
  E_Int dimPb = 3;
  if (nk == 1)
  {
    if (nj > 1) dimPb = 2;
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeDiv2: not valid for 1D structured arrays.");
      RELEASESHAREDS(array,f); return NULL;
    }
  }
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: coordinates not found in array.");
    RELEASESHAREDS(array,f); return NULL;
  }
  posx++; posy++; posz++;

  // Check arrayc
  char* varStringc; char* eltTypec;
  FldArrayF* fc; FldArrayI* cnc;
  E_Int nic, njc, nkc; // number of points of array
  res = K_ARRAY::getFromArray(arrayc, varStringc, fc, nic, njc, nkc, cnc,
                              eltTypec, true);

  // Number of vector fields whose divergence to compute (dimPb components for each)
  E_Int nfld = fc->getNfld(); // total number of scalar fields
  vector<char*> vars;
  K_ARRAY::extractVars(varStringc, vars);
  if (nfld % dimPb != 0)
  {
    RELEASESHAREDB(res,array,f,cn);
    RELEASESHAREDB(res,arrayc,fc,cnc);
    PyErr_SetString(PyExc_TypeError,
                    "computeDiv2: not all components were found for each vector field.");
    return NULL;
  }
  else nfld /= dimPb;

  // Check xyz-plane based only on the first of the given vector fields
  E_Int ixyz = 0; // =0 XY-plane, =1 XZ-plane, =2 YZ-plane
  if (dimPb == 2)
  {
    char* sv0 = vars[0]; char* sv1 = vars[1];
    char s0 = sv0[strlen(sv0)-1];
    char s1 = sv1[strlen(sv1)-1];
    if (s0 == 'X' && s1 == 'Y') ixyz = 0;
    else if (s0 == 'X' && s1 == 'Z') ixyz = 1;
    else if (s0 == 'Y' && s1 == 'Z') ixyz = 2;
    else 
    {
      PyErr_SetString(PyExc_TypeError,
                      "computeDiv2: error with the order of given scalar fields.");
      return NULL;
    }
  }
  else if (dimPb == 3)
  {
    for (E_Int i = 0; i < nfld; i++)
    {
      char* sv0 = vars[dimPb*i]; char* sv1 = vars[dimPb*i+1]; char* sv2 = vars[dimPb*i+2];
      char s0 = sv0[strlen(sv0)-1];
      char s1 = sv1[strlen(sv1)-1];
      char s2 = sv2[strlen(sv2)-1];
      if (s0 != 'X' || s1 != 'Y' || s2 != 'Z') {
        PyErr_SetString(PyExc_TypeError,
                        "computeDiv2: error with the order of given scalar fields.");
        return NULL;
      }
    }
  }

  vector<char*> varStrings;
  for (E_Int i = 0; i < nfld; i++)
  {
    for (E_Int m = 0; m < dimPb-1; m++)
    {
      if (strncmp(vars[dimPb*i+m], vars[dimPb*i+m+1], strlen(vars[dimPb*i+m])-1) != 0)
      {
        RELEASESHAREDB(res,array,f,cn);
        RELEASESHAREDB(res,arrayc,fc,cnc);
        PyErr_SetString(PyExc_TypeError,
                        "computeDiv2: invalid names for vector component fields.");
        return NULL;
      }
    }
    char* local;
    computeDivVarsString(vars[i*dimPb], local);
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
  FldArrayF faceField(nbIntTot,dimPb*nfld); faceField.setAllValuesAtNull();
  FldArrayI voisins(nbIntTot,2); voisins.setAllValuesAt(-1);
  E_Int* cellG = voisins.begin(1);
  E_Int* cellD = voisins.begin(2);
  PyObject* tpl;
  if (dimPb == 2)
    tpl = computeDiv2Struct2D(ni, nj, nic, njc, ixyz, varStringOut, f->begin(posx),  f->begin(posy), f->begin(posz),
                              *fc, faceField, cellG, cellD, indices, fieldX, fieldY, fieldZ); // ATTENTION !!!!!
  else if (dimPb == 3)
    tpl = computeDiv2Struct3D(ni, nj, nk, nic, njc, nkc, varStringOut, f->begin(posx), f->begin(posy), f->begin(posz),
                              *fc, faceField, cellG, cellD, indices, fieldX, fieldY, fieldZ);
  delete [] varStringOut;
  RELEASESHAREDS(array,f);
  RELEASESHAREDS(arrayc,fc);
  return tpl;
}

//=============================================================================
PyObject* K_POST::computeDiv2Struct3D(E_Int ni, E_Int nj, E_Int nk,
                                      E_Int nic, E_Int njc, E_Int nkc,
                                      const char* varStringOut,
                                      E_Float* xt, E_Float* yt, E_Float* zt,
                                      FldArrayF& fc, FldArrayF& faceField,
                                      E_Int* cellG, E_Int* cellD,
                                      PyObject* indices, PyObject* fieldX,
                                      PyObject* fieldY, PyObject* fieldZ )
{
  E_Int indint, indcellg, indcelld;
  E_Int nicnjc = nic*njc;
  E_Int ninjc = ni*njc;
  E_Int nicnj = nic*nj;
  E_Int nbIntI = ninjc*nkc;
  E_Int nbIntJ = nicnj*nkc;
  E_Int nbIntK = nicnjc*nk;
  E_Int nbIntIJ = nbIntI+nbIntJ;
  E_Int nbIntTot = nbIntIJ+nbIntK;
  E_Int nfld = fc.getNfld();
  nfld /= 3; // nfld: num of vector fields
  E_Int nfldg = nfld*3; // nfldg: num of scalar components
  E_Int ncells = nicnjc*nkc;
  for (E_Int eq = 0; eq < nfldg; eq++)
  {
    E_Float* fcn = fc.begin(eq+1);
    E_Float* fintp = faceField.begin(eq+1);
    //faces en i
    for (E_Int k = 0; k < nkc; k++)
      for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 1; i < nic; i++)
        {
          indint = i+j*ni+k*ninjc;
          indcellg = (i-1)+j*nic+k*nicnjc;
          indcelld = indcellg+1;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
        }
    //bords des faces en i
    for (E_Int k = 0; k < nkc; k++)
      for (E_Int j = 0; j < njc; j++)
      {
        E_Int i = 0;
        indint = i+j*ni+k*ninjc;
        indcelld = i + j*nic+k*nicnjc;
        fintp[indint] = fcn[indcelld];//extrapolation de l interieur
        cellG[indint] = -1; cellD[indint] = indcelld;
        i = nic;
        indint = i+j*ni+k*ninjc;
        indcellg = (i-1) + j*nic+k*nicnjc;
        fintp[indint] = fcn[indcellg];//extrapolation de l interieur
        cellG[indint] = indcellg; cellD[indint] = -1;
      }
    //faces en j
    for (E_Int k = 0; k < nkc; k++)
      for (E_Int j = 1; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          indint = i+j*nic+k*nicnj + nbIntI;
          indcellg = i+(j-1)*nic+k*nicnjc;
          indcelld = indcellg+nic;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
        }
    // bords des faces en j
    for (E_Int k = 0; k < nkc; k++)
      for (E_Int i = 0; i < nic; i++)
      {
        E_Int j = 0;
        indint = i+j*nic+k*nicnj+nbIntI;
        indcelld = i+j*nic+k*nicnjc;
        cellG[indint] = -1; cellD[indint] = indcelld;
        fintp[indint] = fcn[indcelld];//extrapolation de l'interieur

        j = njc;
        indint = i+j*nic+k*nicnj+nbIntI;
        indcellg = i+(j-1)*nic+k*nicnjc;
        cellG[indint] = indcellg; cellD[indint] = -1;
        fintp[indint] = fcn[indcellg]; //extrapolation de l'interieur
      }

    //faces en k
    for (E_Int k = 1; k < nkc; k++)
      for (E_Int j = 0; j < njc; j++)
        for (E_Int i = 0; i < nic; i++)
        {
          indint = i+j*nic+k*nicnjc + nbIntIJ;
          indcellg = i+j*nic+(k-1)*nicnjc;
          indcelld = indcellg + nicnjc;
          cellG[indint] = indcellg; cellD[indint] = indcelld;
          fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
        }
    for (E_Int j = 0; j < njc; j++)
      for (E_Int i = 0; i < nic; i++)
      {
        E_Int k = 0;
        indint = i+j*nic+k*nicnjc + nbIntIJ;
        indcelld = i+j*nic+k*nicnjc;;
        cellG[indint] = -1; cellD[indint] = indcelld;
        fintp[indint] = fcn[indcelld];

        k = nkc;
        indint = i+j*nic+k*nicnjc + nbIntIJ;
        indcellg = i+j*nic+(k-1)*nicnjc;
        cellG[indint] = indcellg; cellD[indint] = -1;
        fintp[indint] = fcn[indcellg];
      }
  }
  // Replace DataSet
  if (indices != Py_None && fieldX != Py_None && fieldY != Py_None
                                              && fieldZ != Py_None )
  {
    FldArrayI* inds=NULL; FldArrayF* bfieldX=NULL;
    FldArrayF* bfieldY=NULL; FldArrayF* bfieldZ=NULL;
    K_NUMPY::getFromNumpyArray(indices, inds, true);
    K_NUMPY::getFromNumpyArray(fieldX, bfieldX, true);
    K_NUMPY::getFromNumpyArray(fieldY, bfieldY, true);
    K_NUMPY::getFromNumpyArray(fieldZ, bfieldZ, true);


    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pindint = inds->begin();

    E_Float* pf[3];
    pf[0] = bfieldX->begin();
    pf[1] = bfieldY->begin();
    pf[2] = bfieldZ->begin();

    E_Int indint;
    for (E_Int eq = 0; eq < nfldg; eq++)
    {
      E_Float* fintp = faceField.begin(eq+1);
      for (E_Int noint = 0; noint < ninterfaces; noint++)
      {
        indint = pindint[noint];
        fintp[indint] = pf[eq][noint];
      }
    }
    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(fieldX, bfieldX);
    RELEASESHAREDN(fieldY, bfieldY);
    RELEASESHAREDN(fieldZ, bfieldZ);
  }

  // Build empty array
  PyObject* tpl = K_ARRAY::buildArray(nfld,varStringOut, nic, njc, nkc);
  E_Float* gnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF gp(ncells, nfld, gnp, true); gp.setAllValuesAtNull();

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
  centerInt.malloc(0); surfnorm.malloc(0);

  // divergence
  E_Float ffx, ffy, ffz, voli;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpdv = gp.begin(n+1);
    E_Float* fpx = faceField.begin(3*n+1);
    E_Float* fpy = faceField.begin(3*n+2);
    E_Float* fpz = faceField.begin(3*n+3);
    for (E_Int indint = 0; indint < nbIntTot; indint++)
    {
      E_Int indcellg = cellG[indint]; E_Int indcelld = cellD[indint];
      ffx = fpx[indint];
      ffy = fpy[indint];
      ffz = fpz[indint];
      if (indcellg != -1)
      {
        gpdv[indcellg] += ffx*sxp[indint] + ffy*syp[indint] + ffz*szp[indint];
      }
      if (indcelld != -1)
      {
        gpdv[indcelld] -= ffx*sxp[indint] + ffy*syp[indint] + ffz*szp[indint];
      }
    }
  }
  // free mem
  surf.malloc(0); faceField.malloc(0);

  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpdv = gp.begin(n+1);
    for (E_Int i = 0; i < ncells; i++)
    {
      voli = 1./K_FUNC::E_max(volp[i], 1.e-12);
      gpdv[i] *= voli;
    }
  }
  return tpl;
}
//=============================================================================
PyObject* K_POST::computeDiv2Struct2D(E_Int ni, E_Int nj, E_Int nic, E_Int njc,
                                      E_Int ixyz, const char* varStringOut,
                                      E_Float* xt, E_Float* yt, E_Float* zt,
                                      FldArrayF& fc, FldArrayF& faceField,
                                      E_Int* cellG, E_Int* cellD,
                                      PyObject* indices, PyObject* fieldX,
                                      PyObject* fieldY, PyObject* fieldZ)
{
  E_Int nkc = 1;
  E_Int indint, indcellg, indcelld, indm, indp;
  E_Int nicnjc = nic*njc;
  E_Int ninjc = ni*njc;
  E_Int nicnj = nic*nj;
  E_Int nbIntI = ninjc;
  E_Int nbIntJ = nicnj;
  E_Int nbIntIJ = nbIntI+nbIntJ;
  E_Int nbIntTot = nbIntIJ;
  E_Int nfld = fc.getNfld();
  nfld /= 2; // nfld: num of vector fields
  E_Int nfldg = nfld*2; // nfldg: num of scalar components
  E_Int ncells = nicnjc;
  FldArrayF sint(nbIntTot,3); sint.setAllValuesAtNull();
  E_Float* sxint = sint.begin(1);
  E_Float* syint = sint.begin(2);
  E_Float* szint = sint.begin(3);

 // Compute length of faces
  E_Float d13x, d13y, d13z, d24x, d24y, d24z;

  indint = 0;
  for (E_Int j = 0; j < njc; j++)
    for (E_Int i = 0; i < ni; i++)
    {
      indm = i+j*ni; indp = indm+ni;
      d13x = xt[indp]-xt[indm];
      d13y = yt[indp]-yt[indm];
      d13z = 1;
      d24x = xt[indm]-xt[indp];
      d24y = yt[indm]-yt[indp];
      d24z = 1;

      sxint[indint] = 0.5*(d13y*d24z-d13z*d24y);
      syint[indint] = 0.5*(d13z*d24x-d13x*d24z);
      szint[indint] = 0.5*(d13x*d24y-d13y*d24x);
      indint+=1;
    }
  for (E_Int j = 0; j < nj; j++)
    for (E_Int i = 0; i < nic; i++)
    {
      indm = i+j*ni; indp = indm+1;
      d13x = xt[indp]-xt[indm];
      d13y = yt[indp]-yt[indm];
      d13z = 1;
      d24x = xt[indp]-xt[indm];
      d24y = yt[indp]-yt[indm];
      d24z = -1;

      sxint[indint] = 0.5*(d13y*d24z-d13z*d24y);
      syint[indint] = 0.5*(d13z*d24x-d13x*d24z);
      szint[indint] = 0.5*(d13x*d24y-d13y*d24x);
      indint+=1;
    }

  for (E_Int eq = 0; eq < nfldg; eq++)
  {
    E_Float* fcn = fc.begin(eq+1);
    E_Float* fintp = faceField.begin(eq+1);
    //faces en i internes
    for (E_Int j = 0; j < njc; j++)
      for (E_Int i = 1; i < nic; i++)
      {
        indint = i+j*ni;
        indcellg = (i-1)+j*nic;
        indcelld = indcellg+1;
        cellG[indint] = indcellg; cellD[indint] = indcelld;
        fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
      }
    //bords des faces en i
    for (E_Int j = 0; j < njc; j++)
    {
      //faces i = 0
      indint = j*ni;
      indcelld = j*nic;
      fintp[indint] = fcn[indcelld];//extrapolation de l interieur
      cellG[indint] = -1; cellD[indint] = indcelld;

      //faces i=ni
      indint = (ni-1)+j*ni;
      indcellg = (nic-1) + j*nic;
      fintp[indint] = fcn[indcellg];//extrapolation de l interieur
      cellG[indint] = indcellg; cellD[indint] = -1;
    }

    //faces en j interieures
    for (E_Int j = 1; j < njc; j++)
      for (E_Int i = 0; i < nic; i++)
      {
        indint = i+j*nic+nbIntI;
        indcellg = i+(j-1)*nic;
        indcelld = indcellg+nic;
        cellG[indint] = indcellg; cellD[indint] = indcelld;
        fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
      }
    //bords des faces en j
    for (E_Int i = 0; i < nic; i++)
    {
      //faces j=0
      indint = i+nbIntI;
      indcelld = i;
      cellG[indint] = -1; cellD[indint] = indcelld;
      fintp[indint] = fcn[indcelld];//extrapolation de l interieur

      //faces j=jmax
      indint = i+njc*nic+nbIntI;
      indcellg = i+(njc-1)*nic;
      cellG[indint] = indcellg; cellD[indint] = -1;
      fintp[indint] = fcn[indcellg];//extrapolation de l interieur
    }
  }

  // // Compute length of faces
  // indint = 0;
  // for (E_Int j = 0; j < njc; j++)
  //   for (E_Int i = 0; i < ni; i++)
  //   {
  //     indm = (i-1)+j*ni; indp = indm+ni;
  //     sxint[indint] = xt[indp]-xt[indm];
  //     syint[indint] = yt[indp]-yt[indm];
  //     szint[indint] = zt[indp]-zt[indm];
  //     indint+=1;
  //   }
  // for (E_Int j = 0; j < nj; j++)
  //   for (E_Int i = 0; i < nic; i++)
  //   {
  //     indm = (i-1)+j*ni; indp = indm+1;
  //     sxint[indint] = xt[indp]-xt[indm];
  //     syint[indint] = yt[indp]-yt[indm];
  //     szint[indint] = zt[indp]-zt[indm];
  //     indint+=1;
  //   }
  // for (E_Int eq = 0; eq < nfldg; eq++)
  // {
  //   E_Float* fcn = fc.begin(eq+1);
  //   E_Float* fintp = faceField.begin(eq+1);
  //   //faces en i
  //   for (E_Int j = 0; j < njc; j++)
  //     for (E_Int i = 1; i < nic; i++)
  //     {
  //       indint = i+j*ni;
  //       indcellg = (i-1)+j*nic;
  //       indcelld = indcellg+1;
  //       cellG[indint] = indcellg; cellD[indint] = indcelld;
  //       fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
  //     }
  //   //bords des faces en i
  //   for (E_Int j = 0; j < njc; j++)
  //   {
  //     E_Int i = 0;
  //     indint = i+j*ni;
  //     indcelld = i + j*nic;
  //     fintp[indint] = fcn[indcelld];//extrapolation de l interieur
  //     cellG[indint] = -1; cellD[indint] = indcelld;
  //     i = nic;
  //     indint = i+j*ni;
  //     indcellg = (i-1) + j*nic;
  //     fintp[indint] = fcn[indcellg];//extrapolation de l interieur
  //     cellG[indint] = indcellg; cellD[indint] = -1;
  //   }
  //   //faces en j
  //   for (E_Int j = 1; j < njc; j++)
  //     for (E_Int i = 0; i < nic; i++)
  //     {
  //       indint = i+j*nic+nbIntI;
  //       indcellg = i+(j-1)*nic;
  //       indcelld = indcellg+nic;
  //       cellG[indint] = indcellg; cellD[indint] = indcelld;
  //       fintp[indint] = 0.5*(fcn[indcellg]+fcn[indcelld]);
  //     }
  //   //bords des faces en j
  //   for (E_Int i = 0; i < nic; i++)
  //   {
  //     E_Int j = 0;
  //     indint = i+j*nic+nbIntI;
  //     indcellg = i+(j-1)*nic;
  //     indcelld = indcellg+nic;
  //     cellG[indint] = -1; cellD[indint] = indcelld;
  //     fintp[indint] = fcn[indcelld];//extrapolation de l interieur

  //     j = njc;
  //     indint = i+j*nic+nbIntI;
  //     indcellg = i+(j-1)*nic;
  //     indcelld = indcellg+nic;
  //     cellG[indint] = indcellg; cellD[indint] = -1;
  //     fintp[indint] = fcn[indcellg];//extrapolation de l interieur
  //   }
  // }



  // Replace DataSet
  if (indices != Py_None && (fieldX != Py_None | fieldY != Py_None
			     | fieldZ != Py_None) )
  {
    FldArrayI* inds=NULL; FldArrayF* bfieldX=NULL;
    FldArrayF* bfieldY=NULL; FldArrayF* bfieldZ=NULL;
    K_NUMPY::getFromNumpyArray(indices, inds, true);
    K_NUMPY::getFromNumpyArray(fieldX, bfieldX, true);
    K_NUMPY::getFromNumpyArray(fieldY, bfieldY, true);
    K_NUMPY::getFromNumpyArray(fieldZ, bfieldZ, true);

    E_Int ninterfaces = inds->getSize()*inds->getNfld();
    E_Int* pindint = inds->begin();

    E_Float* pf[2];

    if (ixyz == 0) { pf[0] =  bfieldX->begin(); pf[1] = bfieldY->begin(); }
    if (ixyz == 1) { pf[0] =  bfieldX->begin(); pf[1] = bfieldZ->begin(); }
    if (ixyz == 2) { pf[0] =  bfieldY->begin(); pf[1] = bfieldZ->begin(); }

    E_Int indint;
    for (E_Int eq = 0; eq < nfldg; eq++)
    {
      E_Float* fintp = faceField.begin(eq+1);
      for (E_Int noint = 0; noint < ninterfaces; noint++)
      {
        indint = pindint[noint];
        fintp[indint] = pf[eq][noint];
      }
    }
    RELEASESHAREDN(indices, inds);
    RELEASESHAREDN(fieldX, bfieldX);
    RELEASESHAREDN(fieldY, bfieldY);
    RELEASESHAREDN(fieldZ, bfieldZ);
  }

  // Build empty array
  PyObject* tpl = K_ARRAY::buildArray(nfld,varStringOut, nic, njc, nkc);
  E_Float* gnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF gp(ncells, nfld, gnp, true); gp.setAllValuesAtNull();

  E_Int inti; E_Int intj;
  if (ixyz == 0) {inti = 1; intj = 2;}
  else if (ixyz == 1) {inti = 1; intj = 3;}
  else if (ixyz == 2) {inti = 2; intj = 3;}
  E_Float* siint = sint.begin(inti);
  E_Float* sjint = sint.begin(intj);

  // divergence
  E_Float ffi, ffj;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpdv = gp.begin(n+1);
    E_Float* fpi = faceField.begin(2*n+1);
    E_Float* fpj = faceField.begin(2*n+2);
    for (E_Int indint = 0; indint < nbIntTot; indint++)
    {
      E_Int indcellg = cellG[indint]; E_Int indcelld = cellD[indint];
      ffi = fpi[indint];
      ffj = fpj[indint];
      if (indcellg != -1)
      {
        gpdv[indcellg] += ffi*siint[indint] + ffj*sjint[indint];
      }
      if (indcelld != -1)
      {
        gpdv[indcelld] -= ffi*siint[indint] + ffj*sjint[indint];
      }
    }
  }
  // free mem
  sint.malloc(0); faceField.malloc(0);

  E_Float voli;
  for (E_Int n = 0; n < nfld; n++)
  {
    E_Float* gpdv = gp.begin(n+1);
    for (E_Int indcell = 0; indcell < ncells; indcell++)
    {
      voli =  K_METRIC::compVolOfStructCell2D(ni, nj, xt, yt, zt, indcell, -1);
      voli = 1./K_FUNC::E_max(voli, 1.e-12);
      gpdv[indcell] *= voli;
    }
  }
  return tpl;
}


//=============================================================================
/* From the initial chain of variables: (x,y,z,var1,var2,...)
   Create the chain (divvar1, divvar2, ....)
   This routine allocates varStringOut */
//=============================================================================
void K_POST::computeDivVarsString(char* varString, char*& varStringOut)
{
  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);
  E_Int c = -1;
  E_Int varsSize = vars.size();
  E_Int sizeVarStringOut = 0;
  for (E_Int v = 0; v < varsSize; v++)
  {
    E_Int vsize = strlen(vars[v]);
    sizeVarStringOut+=vsize+4;//divvarString,
  }
  varStringOut = new char [sizeVarStringOut];

  for (E_Int v = 0; v < varsSize; v++)
  {
    char*& var0 = vars[v];
    if (strcmp(var0, "x") != 0 &&
        strcmp(var0, "y") != 0 &&
        strcmp(var0, "z") != 0)
    {
      if (c == -1)
      {
        strcpy(varStringOut, "div");
        c = 1;
      }
      else strcat(varStringOut, ",div");
      char nme[strlen(var0)-1];
      strncpy(nme, var0, strlen(var0)-1);
      nme[strlen(var0)-1] = '\0';
      strcat(varStringOut, nme);
    }
  }
  for (E_Int v = 0; v < varsSize; v++) delete [] vars[v];
}

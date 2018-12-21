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

// Geometry search preconditionning (hook)
 
# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;
# include "ExtArith/quad_double.hpp"
using namespace ExtendedArithmetics;

// ============================================================================
/* Enregistre les centres des faces de a dans un KdTree 
   hook type=0 
   IN: a: array avec coordonnees */
// ============================================================================
PyObject* K_CONVERTER::registerFaces(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray2(array, varString, 
                               f, nil, njl, nkl, cnl, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createHook: array is invalid.");
    return NULL;
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "createHook: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Parcours les faces, calcule les centres, les enregistre dans le KdTree
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  FldArrayF* centers = NULL;
  if (res == 1) // structure
  {
    E_Int ni = nil; E_Int nj = njl; E_Int nk = nkl;
    E_Int nij = ni*nj;
    E_Int ni1 = E_max(ni-1,1);
    E_Int nj1 = E_max(nj-1,1);
    E_Int nk1 = E_max(nk-1,1);
    E_Int nfaces = 0;
    if (nj > 1 && nk > 1) nfaces = ni*nj1*nk1 + ni1*nj*nk1 + ni1*nj1*nk;
    else if (nj > 1) nfaces = ni*nj1*nk1 + ni1*nj*nk1;
    else nfaces = ni*nj1*nk1;
    E_Int ninti = ni*nj1*nk1;
    E_Int nintj = ni1*nj*nk1;
    centers = new FldArrayF(nfaces, 3);
    E_Float* cx = centers->begin(1);
    E_Float* cy = centers->begin(2);
    E_Float* cz = centers->begin(3);
    E_Float inv0 = E_Float(0.25);
    quad_double qinv0 = quad_double(0.25);

#pragma omp parallel default(shared)
    {
    E_Int ind,ind1,ind2,ind3,ind4,ip,jp,kp;
      
    // interface en i
    for (E_Int k = 0; k < nk1; k++)
      for (E_Int j = 0; j < nj1; j++)
#pragma omp for
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i+j*ni+k*ni*nj1;
          ip = E_min(i+1,ni-1);
          jp = E_min(j+1,nj-1);
          kp = E_min(k+1,nk-1);
          ind1 = i+j*ni+k*nij;
          ind2 = i+jp*ni+k*nij;
          ind3 = i+j*ni+kp*nij;
          ind4 = i+jp*ni+kp*nij;

#ifdef QUADDOUBLE
          quad_double qxp[4], qyp[4], qzp[4], qcx, qcy, qcz;
          qxp[0] = quad_double(xp[ind1]);
          qxp[1] = quad_double(xp[ind2]);
          qxp[2] = quad_double(xp[ind3]);
          qxp[3] = quad_double(xp[ind4]);
          qcx = qinv0*(qxp[0]+qxp[1]+qxp[2]+qxp[3]);
          qyp[0] = quad_double(yp[ind1]);
          qyp[1] = quad_double(yp[ind2]);
          qyp[2] = quad_double(yp[ind3]);
          qyp[3] = quad_double(yp[ind4]);
          qcy = qinv0*(qyp[0]+qyp[1]+qyp[2]+qyp[3]);
          qzp[0] = quad_double(zp[ind1]);
          qzp[1] = quad_double(zp[ind2]);
          qzp[2] = quad_double(zp[ind3]);
          qzp[3] = quad_double(zp[ind4]);
          qcz = qinv0*(qzp[0]+qzp[1]+qzp[2]+qzp[3]);
          cx[ind] = E_Float(qcx);
          cy[ind] = E_Float(qcy);
          cz[ind] = E_Float(qcz);
#else
          {
          #ifdef __INTEL_COMPILER
          #pragma float_control(precise, on)
          #endif
          cx[ind] = inv0*(xp[ind1]+xp[ind2]+xp[ind3]+xp[ind4]);
          cy[ind] = inv0*(yp[ind1]+yp[ind2]+yp[ind3]+yp[ind4]);
          cz[ind] = inv0*(zp[ind1]+zp[ind2]+zp[ind3]+zp[ind4]);
          }
#endif
        }
    // interface en j
    if (nj > 1)
    {
      for (E_Int k = 0; k < nk1; k++)
#pragma omp for
        for (E_Int j = 0; j < nj; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind = ninti+i+j*ni1+k*ni1*nj;
            ip = E_min(i+1,ni-1);
            jp = E_min(j+1,nj-1);
            kp = E_min(k+1,nk-1);
            ind1 = i+j*ni+k*nij;
            ind2 = ip+j*ni+k*nij;
            ind3 = i+j*ni+kp*nij;
            ind4 = ip+j*ni+kp*nij;

#ifdef QUADDOUBLE
          quad_double qxp[4], qyp[4], qzp[4], qcx, qcy, qcz;
          qxp[0] = quad_double(xp[ind1]);
          qxp[1] = quad_double(xp[ind2]);
          qxp[2] = quad_double(xp[ind3]);
          qxp[3] = quad_double(xp[ind4]);
          qcx = qinv0*(qxp[0]+qxp[1]+qxp[2]+qxp[3]);
          qyp[0] = quad_double(yp[ind1]);
          qyp[1] = quad_double(yp[ind2]);
          qyp[2] = quad_double(yp[ind3]);
          qyp[3] = quad_double(yp[ind4]);
          qcy = qinv0*(qyp[0]+qyp[1]+qyp[2]+qyp[3]);
          qzp[0] = quad_double(zp[ind1]);
          qzp[1] = quad_double(zp[ind2]);
          qzp[2] = quad_double(zp[ind3]);
          qzp[3] = quad_double(zp[ind4]);
          qcz = qinv0*(qzp[0]+qzp[1]+qzp[2]+qzp[3]);
          cx[ind] = E_Float(qcx);
          cy[ind] = E_Float(qcy);
          cz[ind] = E_Float(qcz);
#else
          {
            #ifdef __INTEL_COMPILER
            #pragma float_control(precise, on)
            #endif
            cx[ind] = inv0*(xp[ind1]+xp[ind2]+xp[ind3]+xp[ind4]);
            cy[ind] = inv0*(yp[ind1]+yp[ind2]+yp[ind3]+yp[ind4]);
            cz[ind] = inv0*(zp[ind1]+zp[ind2]+zp[ind3]+zp[ind4]);
          }
#endif
          }
    }
    // interface en k
    if (nk > 1)
    {
#pragma omp for
      for (E_Int k = 0; k < nk; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind = ninti+nintj+i+j*ni1+k*ni1*nj1;
            ip = E_min(i+1,ni-1);
            jp = E_min(j+1,nj-1);
            kp = E_min(k+1,nk-1);
            ind1 = i+j*ni+k*nij;
            ind2 = ip+j*ni+k*nij;
            ind3 = i+jp*ni+k*nij;
            ind4 = ip+jp*ni+k*nij;
    
#ifdef QUADDOUBLE
          quad_double qxp[4], qyp[4], qzp[4], qcx, qcy, qcz;
          qxp[0] = quad_double(xp[ind1]);
          qxp[1] = quad_double(xp[ind2]);
          qxp[2] = quad_double(xp[ind3]);
          qxp[3] = quad_double(xp[ind4]);
          qcx = qinv0*(qxp[0]+qxp[1]+qxp[2]+qxp[3]);
          qyp[0] = quad_double(yp[ind1]);
          qyp[1] = quad_double(yp[ind2]);
          qyp[2] = quad_double(yp[ind3]);
          qyp[3] = quad_double(yp[ind4]);
          qcy = qinv0*(qyp[0]+qyp[1]+qyp[2]+qyp[3]);
          qzp[0] = quad_double(zp[ind1]);
          qzp[1] = quad_double(zp[ind2]);
          qzp[2] = quad_double(zp[ind3]);
          qzp[3] = quad_double(zp[ind4]);
          qcz = qinv0*(qzp[0]+qzp[1]+qzp[2]+qzp[3]);
          cx[ind] = E_Float(qcx);
          cy[ind] = E_Float(qcy);
          cz[ind] = E_Float(qcz);            
#else
          {
            #ifdef __INTEL_COMPILER
            #pragma float_control(precise, on)
            #endif
            cx[ind] = inv0*(xp[ind1]+xp[ind2]+xp[ind3]+xp[ind4]);
            cy[ind] = inv0*(yp[ind1]+yp[ind2]+yp[ind3]+yp[ind4]);
            cz[ind] = inv0*(zp[ind1]+zp[ind2]+zp[ind3]+zp[ind4]);
          }
#endif
          }
    }
    }

  }
  else if (res == 2 && strcmp(eltType, "NGON") == 0) // NGON
  {
    E_Int nfaces = cnl->getNFaces();
    centers = new FldArrayF(nfaces, 3);
    E_Float* cx = centers->begin(1);
    E_Float* cy = centers->begin(2);
    E_Float* cz = centers->begin(3);
    E_Int* ptr;
    if (cnl->isNGon() == 2) ptr = cnl->getNGon();
    else ptr = cnl->begin();
    FldArrayI posFace; K_CONNECT::getPosFaces(*cnl, posFace);

# pragma omp parallel for default(shared)
    for (E_Int i = 0; i < nfaces; i++)
    {
      E_Int posf = posFace[i];
      E_Int* ptrFace = &ptr[posf];
      E_Int nv = ptrFace[0];
      E_Float xf=0., yf=0., zf=0.;
      
#ifdef QUADDOUBLE
      quad_double qxf, qyf, qzf;
      quad_double qinv = quad_double(nv);
      for (E_Int n = 1; n <= nv; n++)
      { 
        E_Int ind = ptrFace[n]-1; 
        qxf = qxf+quad_double(xp[ind]); 
        qyf = qyf+quad_double(yp[ind]); 
        qzf = qzf+quad_double(zp[ind]); 
      }
      qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
#else
      {
      #ifdef __INTEL_COMPILER
      #pragma float_control(precise, on)
      #endif
      for (E_Int n = 1; n <= nv; n++)
      { 
        E_Int ind = ptrFace[n]-1; xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      }
      E_Float inv = 1./E_Float(nv); xf *= inv; yf *= inv; zf *= inv;
      }
#endif

      cx[i] = xf; cy[i] = yf; cz[i] = zf;
    }
  }
  else // Basic elements
  {
    E_Int nelts = cnl->getSize();
    E_Int nof = 0; E_Int nfaces = 0;
    E_Int face[6][4];
    if (strcmp(eltType, "BAR") == 0) 
    {
      nfaces = 2; nof = 1;
      face[0][0] = 1; face[1][0] = 2;
    }
    else if (strcmp(eltType, "QUAD") == 0) 
    {
      nfaces = 4; nof = 2;
      face[0][0] = 1; face[0][1] = 2;
      face[1][0] = 2; face[1][1] = 3;
      face[2][0] = 3; face[2][1] = 4;
      face[3][0] = 4; face[3][1] = 1;
    }
    else if (strcmp(eltType, "TRI") == 0) 
    {
      nfaces = 3; nof = 2;
      face[0][0] = 1; face[0][1] = 2;
      face[1][0] = 2; face[1][1] = 3;
      face[2][0] = 3; face[2][1] = 1;
    }
    else if (strcmp(eltType, "HEXA") == 0) 
    {
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
      nfaces = 4; nof = 3;
      face[0][0] = 1; face[0][1] = 3; face[0][2] = 2;
      face[1][0] = 1; face[1][1] = 2; face[1][2] = 4;
      face[2][0] = 2; face[2][1] = 3; face[2][2] = 4;
      face[3][0] = 3; face[3][1] = 1; face[3][2] = 4;
    }
    else if (strcmp(eltType, "PYRA") == 0) 
    {
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
      nfaces = 5; nof = 4; // TRI degen
      face[0][0] = 1; face[0][1] = 2; face[0][2] = 5; face[0][3] = 4;
      face[1][0] = 2; face[1][1] = 3; face[1][2] = 6; face[1][3] = 5;
      face[2][0] = 3; face[2][1] = 1; face[2][2] = 4; face[2][3] = 6;
      face[3][0] = 1; face[3][1] = 3; face[3][2] = 2; face[3][3] = 2;
      face[4][0] = 4; face[4][1] = 5; face[4][2] = 6; face[4][3] = 6;
    }
    
    centers = new FldArrayF(nelts*nfaces, 3);
    E_Float* cx = centers->begin(1);
    E_Float* cy = centers->begin(2);
    E_Float* cz = centers->begin(3);
    E_Float inv = E_Float(1./nof);
    quad_double qinv = quad_double(nof);

#pragma omp parallel default(shared)
    {
      E_Int ind, indl;

#pragma omp for  
    for (E_Int i = 0; i < nelts; i++)
    {
      for (E_Int f = 0; f < nfaces; f++)
      {
        ind = f + i*nfaces;

#ifdef QUADDOUBLE
        quad_double qcx, qcy, qcz;
        for (E_Int n = 0; n < nof; n++)
        {
          indl = (*cnl)(i,face[f][n])-1;
          qcx = qcx+quad_double(xp[indl]); 
          qcy = qcy+quad_double(yp[indl]); 
          qcz = qcz+quad_double(zp[indl]);
        }
        qcx = qcx/qinv; qcy = qcy/qinv; qcz = qcz/qinv;
        cx[ind] = E_Float(qcx);
        cy[ind] = E_Float(qcy);
        cz[ind] = E_Float(qcz);
#else
        {
        #ifdef __INTEL_COMPILER
        #pragma float_control(precise, on)
        #endif
        cx[ind] = 0.; cy[ind] = 0.; cz[ind] = 0.;
        for (E_Int n = 0; n < nof; n++)
        {
          indl = (*cnl)(i,face[f][n])-1;
          cx[ind] += xp[indl]; cy[ind] += yp[indl]; cz[ind] += zp[indl];
        }
        cx[ind] *= inv; cy[ind] *= inv; cz[ind] *= inv;
        }
#endif
      }// loop on faces 
    }// loop on elts
    }//omp
  }

  ArrayAccessor<FldArrayF>* coordAcc = 
    new ArrayAccessor<FldArrayF>(*centers, 1, 2, 3); // ref sur centers
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    new K_SEARCH::KdTree<FldArrayF>(*coordAcc); // ref sur coordAcc

  PyObject* hook;
  E_Int* type = new E_Int [1]; type[0] = 0;
  E_Int sizePacket = 4;
  void** packet = new void* [sizePacket];
  packet[0] = type; // hook type
  packet[1] = centers;
  packet[2] = coordAcc;
  packet[3] = globalKdt;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif

  RELEASESHAREDB(res, array, f, cnl);
  return hook;
}

// ============================================================================
/* Enregistre les BB des cellules de a dans un ADT 
   hook type=1
   IN: arrays: structure ou non-structure (TETRA) */
// ============================================================================
PyObject* K_CONVERTER::registerCells(PyObject* self, PyObject* args)
{
  PyObject* listFields;
  if (!PyArg_ParseTuple(args, "O", &listFields)) return NULL;

  if (PyList_Check(listFields) == false)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createHook: arg must be a list of arrays.");
    return NULL;
  }

  E_Int  posxi, posyi, poszi;
  // Extract infos from listFields
  vector<E_Int> resl;  vector<char*> varString;
  vector<FldArrayF*> fields;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  E_Boolean skipNoCoord = true; E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false; E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(listFields, resl, varString, fields, 
                                      a2, a3, a4, objs, skipDiffVars, skipNoCoord, 
                                      skipStructured, skipUnstructured, true);
  E_Int nzones = resl.size();
  if (isOk == -1)
  {
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);  
    PyErr_SetString(PyExc_TypeError,
                    "createHook: invalid list of arrays.");
    return NULL;
  }
  E_Int nfldTot = fields[0]->getNfld();
  
  // Verification du nb de champs donnes
  for (E_Int i = 0; i < nzones; i++)
  {
    E_Int nfld0 = fields[i]->getNfld();
    if (nfld0 != nfldTot)
    {
      for (E_Int no = 0; no < nzones; no++)
        RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]); 
      PyErr_SetString(PyExc_TypeError,
                      "createHook: input arrays must have the same number of variables.");  
      return NULL;
    }
  }

  // Verification de posxi, posyi, poszi dans listFields: vars 1,2,3 imposees
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs; 
  for (E_Int no = 0; no < nzones; no++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(varString[no]); posxi++;
    posyi = K_ARRAY::isCoordinateYPresent(varString[no]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(varString[no]); poszi++;
    if (posxi != 1 || posyi != 2 || poszi != 3)
    {
      for (E_Int noi = 0; noi < nzones; noi++)
        RELEASESHAREDA(resl[noi],objs[noi],fields[noi],a2[noi],a3[noi],a4[noi]); 
      PyErr_SetString(PyExc_TypeError,
                      "createHook: coordinates in input arrays must be variables 1,2,3.");
      return NULL;
    }
    posxs.push_back(posxi); posys.push_back(posyi); poszs.push_back(poszi); 
  }
  // Verif des zones pour l adt
  for (E_Int no = 0; no < nzones; no++)
  {
    if (resl[no] == 2) 
    {
      char* eltType0 = (char*)a3[no];
      if (K_STRING::cmp(eltType0,"TETRA")!= 0)
      {
        for (E_Int noi = 0; noi < nzones; noi++)
          RELEASESHAREDA(resl[noi],objs[noi],fields[noi],a2[noi],a3[noi],a4[noi]); 
        PyErr_SetString(PyExc_TypeError,
                        "setInterpData: unstructured zones must be TETRA.");
        return NULL;
      }
    }
    else if (resl[no] == 1) 
    {
      if (*(E_Int*)a2[no]<2 || *(E_Int*)a3[no]<2 ) 
      {
        for (E_Int noi = 0; noi < nzones; noi++)
          RELEASESHAREDA(resl[noi],objs[noi],fields[noi],a2[noi],a3[noi],a4[noi]); 
        PyErr_SetString(PyExc_TypeError,
                        "setInterpData: structured donor zones must be 3D or nk=1.");
        return NULL;
      }
    }
  }
  // Liste des interpDatas
  vector<K_INTERP::InterpAdt*> interpDatas;
  E_Int isBuilt;
  for (E_Int no = 0; no < nzones; no++)
  {
    K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(
      fields[no]->getSize(), 
      fields[no]->begin(posxs[no]),
      fields[no]->begin(posys[no]),
      fields[no]->begin(poszs[no]),
      a2[no], a3[no], a4[no], isBuilt);
    if (isBuilt == 1) 
      interpDatas.push_back(adt);
    else 
    {
      for (E_Int noi = 0; noi < nzones; noi++)
        RELEASESHAREDA(resl[noi],objs[noi],fields[noi],a2[noi],a3[noi],a4[noi]); 
      for (size_t noi = 0; noi < interpDatas.size(); noi++)
        delete interpDatas[noi];
      PyErr_SetString(PyExc_TypeError,
                      "createHook: 2D structured zones must be z=constant.");
      return NULL;
    }
  }
  // Build hook
  E_Int s1 = nzones;
  PyObject* hook;
  E_Int* type = new E_Int [2]; type[0] = 1; type[1] = s1; 
  E_Int sizePacket = s1 + 1;
  void** packet = new void* [sizePacket];
  packet[0] = type; // hook type
  for (E_Int i = 0; i < s1; i++) packet[i+1] = (void*)interpDatas[i];

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif

  for (E_Int no = 0; no < nzones; no++)
    RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]); 
  return hook;
}

// ============================================================================
/* Enregistre les noeuds a dans un KdTree 
   hook type=2
   IN: a: tout type d'array avec coordonnees */
// ============================================================================
PyObject* K_CONVERTER::registerNodes(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray2(array, varString, 
                              f, nil, njl, nkl, cnl, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createHook: array is invalid.");
    return NULL;
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "createHook: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Parcours noeuds, les enregistre dans le KdTree
  E_Int nodes = f->getSize();
  FldArrayF* coords = new FldArrayF(nodes, 3);
  E_Float* cx = coords->begin(1);
  E_Float* cy = coords->begin(2);
  E_Float* cz = coords->begin(3);
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);

#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < nodes; i++)
  {
    cx[i] = xp[i]; cy[i] = yp[i]; cz[i] = zp[i];
  }
  ArrayAccessor<FldArrayF>* coordAcc = 
    new ArrayAccessor<FldArrayF>(*coords, 1, 2, 3); // ref sur coords
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    new K_SEARCH::KdTree<FldArrayF>(*coordAcc); // ref sur coordAcc

  PyObject* hook;
  E_Int* type = new E_Int [1]; type[0] = 2;
  E_Int sizePacket = 4;
  void** packet = new void* [sizePacket];
  packet[0] = type; // hook type
  packet[1] = coords;
  packet[2] = coordAcc;
  packet[3] = globalKdt;

#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  RELEASESHAREDB(res, array, f, cnl);
  return hook;
}

// ============================================================================
/* Enregistre les centres des elements de a dans un KdTree 
   hook type=3 
   IN: a: array NGON ou Basic Elements */
// ============================================================================
PyObject* K_CONVERTER::registerElements(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray2(array, varString, 
                               f, nil, njl, nkl, cnl, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createHook: array is invalid.");
    return NULL;
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "createHook: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  FldArrayF* centers = NULL;
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  // Parcours les elements, calcule les centres, les enregistre dans le KdTree
  if (res == 1) // structure
  {
    E_Int ni1 = E_max(nil-1,1);
    E_Int nj1 = E_max(njl-1,1);
    E_Int ni1nj1 = ni1*nj1;
    E_Int nk1 = E_max(nkl-1,1);
    E_Int nelts = ni1nj1*nk1;
    E_Int nij = nil*njl;
    centers = new FldArrayF(nelts,3);
    E_Float* cx = centers->begin(1);
    E_Float* cy = centers->begin(2);
    E_Float* cz = centers->begin(3);
    E_Float inv = E_Float(0.125);
    quad_double qinv = quad_double(8.);

#pragma omp parallel for default(shared)
    for (E_Int k = 0; k < nk1; k++)
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          E_Int ip = E_min(i+1,nil-1); 
          E_Int jp = E_min(j+1,njl-1); 
          E_Int kp = E_min(k+1,nkl-1);

          E_Int indT[8];
          E_Int indcell = i+j*ni1+k*ni1nj1;
          indT[0] = i  + j*nil  + k*nij;
          indT[1] = ip + j*nil  + k*nij;
          indT[2] = i  + jp*nil + k*nij;
          indT[3] = ip + jp*nil + k*nij;
          indT[4] = i  + j*nil  + kp*nij;
          indT[5] = ip + j*nil  + kp*nij;
          indT[6] = i  + jp*nil + kp*nij;
          indT[7] = ip + jp*nil + kp*nij;

          E_Float xf=0., yf=0.,zf=0.;
          quad_double qxf, qyf, qzf;

#ifdef QUADDOUBLE
          for (E_Int nov = 0; nov < 8; nov++)
          {
            E_Int indv = indT[nov];
            qxf = qxf+quad_double(xp[indv]); 
            qyf = qyf+quad_double(yp[indv]); 
            qzf = qzf+quad_double(zp[indv]); 
          }
          qxf =qxf/qinv;
          qyf =qyf/qinv;
          qzf =qzf/qinv;
          cx[indcell]=E_Float(qxf);
          cy[indcell]=E_Float(qyf);
          cz[indcell]=E_Float(qzf);
#else
          {
          #ifdef __INTEL_COMPILER
          #pragma float_control(precise, on)
          #endif
          for (E_Int nov = 0; nov < 8; nov++)
          {
            E_Int indv = indT[nov];
            xf += xp[indv]; yf += yp[indv]; zf+= zp[indv];
          }
          xf*=inv; yf*=inv; zf*=inv;
          cx[indcell]=xf; cy[indcell]=yf; cz[indcell]=zf;
          }
#endif     
        }
  }
  else if (res == 2 && strcmp(eltType, "NGON") == 0)
  {
    E_Int nelts = cnl->getNElts();
    centers = new FldArrayF(nelts, 3);
    E_Float* cx = centers->begin(1);
    E_Float* cy = centers->begin(2);
    E_Float* cz = centers->begin(3);
    E_Int* ptrf; E_Int* ptre;
    if (cnl->isNGon() == 2) ptrf = cnl->getNGon();
    else ptrf = cnl->begin();
    if (cnl->isNGon() == 2) ptre = cnl->getNFace();
    else ptre = cnl->begin();
    FldArrayI posFace; K_CONNECT::getPosFaces(*cnl, posFace);
    FldArrayI posElt; K_CONNECT::getPosElts(*cnl, posElt);
    
# pragma omp parallel for default(shared)
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int pose = posElt[i];
      E_Int* ptrElt = &ptre[pose];
      E_Int nf = ptrElt[0];
      E_Float xf=0., yf=0.,zf=0.;
      E_Int c = 0;
      quad_double qxf, qyf, qzf;

      for (E_Int n = 1; n <= nf; n++)
      { 
        E_Int ind = ptrElt[n]-1;
        E_Int pos = posFace[ind];
        E_Int* ptrFace = &ptrf[pos];
        E_Int nv = ptrFace[0];

#ifdef QUADDOUBLE
        for (E_Int p = 1; p <= nv; p++)
        {
          ind = ptrFace[p]-1; 
          qxf = qxf+quad_double(xp[ind]); 
          qyf = qyf+quad_double(yp[ind]); 
          qzf = qzf+quad_double(zp[ind]); 
          c++;
        }
      }
      quad_double qinv = quad_double(c);
      qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);

#else
        //#pragma float_control(precise, on)
        for (E_Int p = 1; p <= nv; p++)
        {
          ind = ptrFace[p]-1; xf += xp[ind]; yf += yp[ind]; zf += zp[ind]; c++;
        }
      }
      E_Float inv = 1./E_Float(c); xf *= inv; yf *= inv; zf *= inv;
#endif
      
      cx[i] = xf; cy[i] = yf; cz[i] = zf;
  } // loop on elts
  }// NGON
  else
  {
    E_Int nelts = cnl->getSize();
    E_Int nv = cnl->getNfld();
    E_Float inv = E_Float(1./nv);
    quad_double qinv = quad_double(nv);
    centers = new FldArrayF(nelts, 3);
    E_Float* cx = centers->begin(1);
    E_Float* cy = centers->begin(2);
    E_Float* cz = centers->begin(3);

#pragma omp parallel default(shared)
    {
      E_Int ind;

#pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Float xf=0., yf=0., zf=0.;

#ifdef QUADDOUBLE
      quad_double qxf, qyf, qzf;
      for (E_Int j = 1; j <= nv; j++)
      {
        ind = (*cnl)(i,j)-1;
        qxf = qxf+quad_double(xp[ind]); 
        qyf = qyf+quad_double(yp[ind]); 
        qzf = qzf+quad_double(zp[ind]); 
      }
      qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      cx[i] = E_Float(qxf); cy[i] = E_Float(qyf); cz[i] = E_Float(qzf);
#else
      //#pragma float_control(precise, on) 
      for (E_Int j = 1; j <= nv; j++)
      {
        ind = (*cnl)(i,j)-1;
        xf += xp[ind]; yf+= yp[ind]; zf += zp[ind];
      }
      xf *= inv; yf*= inv; zf*= inv; 
      cx[i] = xf; cy[i] = yf; cz[i] = zf;
#endif
    }//loop for elts
    }// instruct omp
  }//EB
  
  ArrayAccessor<FldArrayF>* coordAcc = 
    new ArrayAccessor<FldArrayF>(*centers, 1, 2, 3); // ref sur centers
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    new K_SEARCH::KdTree<FldArrayF>(*coordAcc); // ref sur coordAcc

  PyObject* hook;
  E_Int* type = new E_Int [1]; type[0] = 3;

  E_Int sizePacket = 4;
  void** packet = new void* [sizePacket];
  packet[0] = type; // hook type
  packet[1] = centers;
  packet[2] = coordAcc;
  packet[3] = globalKdt;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  hook = PyCObject_FromVoidPtr(packet, NULL);
#else
  hook = PyCapsule_New(packet, NULL, NULL);
#endif
  RELEASESHAREDB(res, array, f, cnl);
  return hook;
}

//=============================================================================
/* 
   Fonction generale de free hook
 */
//=============================================================================
PyObject* K_CONVERTER::freeHook(PyObject* self, PyObject* args)
{
  PyObject* hook;
  if (!PyArg_ParseTuple(args, "O", &hook)) return NULL;
  
  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* typep = (E_Int*)packet[0]; // type of hook
  E_Int type = *typep;
  
  switch (type)
  {
    case 0:
    case 2:
    case 3:
    case 100:
    case 102:
    case 103:
    {
      // KDT (0;2;3)
      FldArrayF* pt = (FldArrayF*)packet[1];
      //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
      //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
      K_SEARCH::KdTree<FldArrayF>* globalKdt = 
        (K_SEARCH::KdTree<FldArrayF>*) packet[3];
      
      delete pt;
      delete globalKdt; // delete tout
    }
    break;

    case 1:
    {
      // Cells BB in ADT (1)
      K_INTERP::InterpAdt* d;
      E_Int s1 = typep[1]; 
      for (E_Int i = 0; i < s1; i++) 
      {
        d = (K_INTERP::InterpAdt*)(packet[i+1]);
        delete d; // delete kmesh also
      }
    }
    break;

    case 4:
    {
      // METRIC Structuree: sx, sy, sz, surfno (4)
      FldArrayF* sx = (FldArrayF*)packet[1];
      FldArrayF* sy = (FldArrayF*)packet[2];
      FldArrayF* sz = (FldArrayF*)packet[3];
      FldArrayF* sno = (FldArrayF*)packet[4]; 
      delete sx; delete sy; delete sz; delete sno;
    }
    break;

    case 5:
    {
      // METRIC NGON: sx, sy, sz, surfno, cFE (5)
      FldArrayF* sx = (FldArrayF*)packet[1];
      FldArrayF* sy = (FldArrayF*)packet[2];
      FldArrayF* sz = (FldArrayF*)packet[3];
      FldArrayF* sno = (FldArrayF*)packet[4];
      FldArrayI* cFE = (FldArrayI*)packet[5];
      delete sx; delete sy; delete sz; delete sno; delete cFE;
    }
    break;
  }

  delete [] typep;
  delete [] packet;
  Py_INCREF(Py_None);
  return Py_None;
}

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

// Routines d'identification geometrique

# include "converter.h"
using namespace K_FUNC;
using namespace K_FLD;
# include "ExtArith/quad_double.hpp"
using namespace ExtendedArithmetics;

// ============================================================================
/* Identifie les noeuds de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: zone a identfier
   Retourne la liste des noeuds dans la numerotation du KDT. */
// ============================================================================
PyObject* K_CONVERTER::identifyNodes(PyObject* self, PyObject* args)
{  
  PyObject* array; PyObject* hook; E_Float tol;
  if (!PYPARSETUPLEF(args, "OOd", "OOf", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3 &&
      *type != 100 && *type != 102 && *type != 103)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray2(array, varString, 
                              f, nil, njl, nkl, cnl, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: array is invalid.");
    return NULL;
  }

  E_Int posx, posy, posz;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res,array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyNodes: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int npts = f->getSize();
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);

  PyObject* ac = K_NUMPY::buildNumpyArray(npts, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Remplissage
#pragma omp parallel default(shared)
  {
  E_Float pt[3];
  E_Float xf,yf,zf,dx,dy,dz;
  E_Int ind;
#pragma omp for schedule(dynamic)
  for (E_Int i = 0; i < npts; i++)
  {
    xf = xp[i]; yf = yp[i]; zf = zp[i];
    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    ind = globalKdt->getClosest(pt); // closest pt
    dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
    if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
    else nptr[i] = -1;
  }
  }
  RELEASESHAREDB(res, array, f, cnl);
  return ac;
}

// ============================================================================
/* Identifie les centres des faces de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne la liste des faces. */
// ============================================================================
PyObject* K_CONVERTER::identifyFaces(PyObject* self, PyObject* args)
{ 
  PyObject* array; PyObject* hook;
  E_Float tol;
  if (!PYPARSETUPLEF(args, "OOd", "OOf", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3&&
      *type != 100 && *type != 102 && *type != 103)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray2(array, varString, 
                              f, nil, njl, nkl, cnl, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must be a NGON.");
    return NULL; 
  }
  if (strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must be a NGON.");
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
                    "identifyFaces: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int nfaces = (*cnl)[0];
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);
  
  PyObject* ac = K_NUMPY::buildNumpyArray(nfaces, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Remplissage
  FldArrayI posFace; K_CONNECT::getPosFaces(*cnl, posFace);
  E_Int* posFacep = posFace.begin();
  E_Int* cnp;
  if (cnl->isNGon() == 2) cnp = cnl->getNGon();
  else cnp = cnl->begin();

#pragma omp parallel default(shared)
  {
#pragma omp for schedule(dynamic)
  for (E_Int i = 0; i < nfaces; i++)
  {
    E_Int pos = posFacep[i];
    E_Int* ptr = &cnp[pos];
    E_Int nv = ptr[0];
    E_Float xf = 0.,yf = 0.,zf = 0.;

#ifdef QUADDOUBLE
    quad_double qxf, qyf, qzf;
    for (E_Int n = 1; n <= nv; n++)
    {
      E_Int ind = ptr[n]-1; 
      qxf = qxf+quad_double(xp[ind]); 
      qyf = qyf+quad_double(yp[ind]); 
      qzf = qzf+quad_double(zp[ind]); 
    }
    quad_double qinv = quad_double(nv); qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
    xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
#else
    {
    #ifdef __INTEL_COMPILER
    #pragma float_control(precise, on)
    #endif
    for (E_Int n = 1; n <= nv; n++)
    {
      E_Int ind = ptr[n]-1; 
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
    }
    E_Float inv = 1./E_Float(nv); xf *= inv; yf *= inv; zf *= inv;
    }
#endif
    E_Float pt[3];
    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    E_Int indk = globalKdt->getClosest(pt); // closest pt
    E_Float dx = xt[indk]-xf; E_Float dy = yt[indk]-yf; E_Float dz = zt[indk]-zf;
    // std::cout << "dx, dy, dz " << dx << " " << dy << " " << dz << std::endl;
    if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = indk+1; 
    else nptr[i] = -1;
  }
  }

  RELEASESHAREDU(array, f, cnl);
  return ac;
}

// ============================================================================
/* Identifie les centres des elements de a avec les points 
   stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne les indices des points du kdtree correspondant. */
// ============================================================================
PyObject* K_CONVERTER::identifyElements(PyObject* self, PyObject* args)
{ 
  PyObject* array; PyObject* hook;
  E_Float tol;
  if (!PyArg_ParseTuple(args, "OOd", &hook, &array, &tol)) return NULL;

  // recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 0 && *type != 2 && *type != 3&&
      *type != 100 && *type != 102 && *type != 103) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  //K_SEARCH::KdTree<FldArrayF>* coordAcc = 
  //  (K_SEARCH::KdTree<FldArrayF>*) packet[2];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray2(array, varString, 
                              f, nil, njl, nkl, cnl, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElts: array is invalid.");
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
                    "identifyElts: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;
  PyObject* ac = NULL;

  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);

  if (res == 1)
  {
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);
    E_Int nie = K_FUNC::E_max(nil-1,1);
    E_Int nje = K_FUNC::E_max(njl-1,1);
    E_Int nke = K_FUNC::E_max(nkl-1,1);
    E_Int nijl = nil*njl;
    
    ac = K_NUMPY::buildNumpyArray(nie*nje*nke, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);
    E_Int inci = K_FUNC::E_min(1,nil-1);
    E_Int incj = K_FUNC::E_min(1,njl-1);
    E_Int inck = K_FUNC::E_min(1,nkl-1);
#pragma omp parallel default(shared)
    {
    E_Float xf,yf,zf,dx,dy,dz;
    E_Int ind;
    E_Float pt[3];
    E_Int ic = 0;

    for (E_Int k = 0; k < nke; k++)
#pragma omp for
    for (E_Int j = 0; j < nje; j++)
    for (E_Int i = 0; i < nie; i++)
    {
      ic = i+j*nie+k*nie*nje;
      xf = 0.; yf = 0.; zf = 0.;
      ind = i+j*nil+k*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      ind = i+inci+j*nil+k*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      ind = i+(j+incj)*nil+k*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      ind = i+inci+(j+incj)*nil+k*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      ind = i+j*nil+(k+inck)*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      ind = i+inci+j*nil+(k+inck)*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      ind = i+(j+incj)*nil+(k+inck)*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
      ind = i+inci+(j+incj)*nil+(k+inck)*nijl;
      xf += xp[ind]; yf += yp[ind]; zf += zp[ind];

      xf = xf/8.; yf = yf/8.; zf = zf/8.;
      pt[0] = xf; pt[1] = yf; pt[2] = zf;

      ind = globalKdt->getClosest(pt); // closest pt
      dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
      //d = dx*dx + dy*dy + dz*dz;
      //if (d < tol) nptr[i] = ind+1;
      //else nptr[i] = -1;
      if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[ic] = ind+1;
      else nptr[ic] = -1;
    }
    }
  }
  else if (strcmp(eltType, "NGON") == 0)
  {  
    // Cree le numpy de sortie
    E_Int nelts = cnl->getNElts();
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);

    ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

    // Remplissage
    FldArrayI posFace; K_CONNECT::getPosFaces(*cnl, posFace);
    E_Int* posFacep = posFace.begin();
    FldArrayI posElts; K_CONNECT::getPosElts(*cnl, posElts);
    E_Int* posEltsp = posElts.begin();

    E_Int* ptrf; E_Int* ptre;
    if (cnl->isNGon() == 2) { ptrf = cnl->getNGon(); ptre = cnl->getNFace(); }
    else { ptrf = cnl->begin(); ptre = cnl->begin(); }
    
#pragma omp parallel default(shared)
    {
      E_Float pt[3];
#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int* ptrElt = &ptre[posEltsp[i]];
      E_Int nf = ptrElt[0];
      E_Float xf = 0.; E_Float yf = 0.; E_Float zf = 0.; E_Int c = 0;
#ifdef QUADDOUBLE
      quad_double qxf, qyf, qzf;
      for (E_Int n = 1; n <= nf; n++)
      { 
        E_Int ind = ptrElt[n]-1;
        E_Int pos = posFacep[ind];
        E_Int* ptrFace = &ptrf[pos];
        E_Int nv = ptrFace[0];
        for (E_Int p = 1; p <= nv; p++)
        {
          E_Int indp = ptrFace[p]-1; 
          qxf = qxf+quad_double(xp[indp]); 
          qyf = qyf+quad_double(yp[indp]);
          qzf = qzf+quad_double(zp[indp]); 
          c++;
        }
      }
      quad_double qinv = quad_double(c); qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
#else
      {
      #ifdef __INTEL_COMPILER
      #pragma float_control(precise, on)
      #endif
      for (E_Int n = 1; n <= nf; n++)
      { 
        E_Int ind = ptrElt[n]-1;
        E_Int pos = posFacep[ind];
        E_Int* ptrFace = &ptrf[pos];
        E_Int nv = ptrFace[0];
        for (E_Int p = 1; p <= nv; p++)
        {
          E_Int indp = ptrFace[p]-1; xf += xp[indp]; yf += yp[indp]; zf += zp[indp]; c++;
        }
      }
      E_Float inv = 1./E_Float(c); xf *= inv; yf *= inv; zf *= inv;
      }
#endif
      pt[0] = xf; pt[1] = yf; pt[2] = zf;
      E_Int ind = globalKdt->getClosest(pt); // closest pt
      E_Float dx = xt[ind]-xf; E_Float dy = yt[ind]-yf; E_Float dz = zt[ind]-zf;
      if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
      else nptr[i] = -1;
    }
    }
  }
  else if (strcmp(eltType,"BAR")==0 || strcmp(eltType,"TRI")==0 || 
           strcmp(eltType,"QUAD")==0 || strcmp(eltType,"TETRA")==0 || 
           strcmp(eltType,"HEXA")==0 || strcmp(eltType,"PENTA")==0 ||
           strcmp(eltType,"PYRA")==0)
  {
    E_Int nelts = cnl->getSize();
    E_Int nvert = cnl->getNfld();
    
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);
    ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);
    E_Float inv = 1./E_Float(nvert);
    quad_double qinv = quad_double(nvert);

#pragma omp parallel default(shared)
    {
      E_Float pt[3];
      E_Float xf,yf,zf,dx,dy,dz;
      E_Int ind;
#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < nelts; i++)
    {
      xf = 0.; yf = 0.; zf = 0.;

#ifdef QUADDOUBLE
      quad_double qxf, qyf, qzf;
      for (E_Int n = 1; n <= nvert; n++)
      {
        ind = (*cnl)(i,n)-1; 
        qxf = qxf+quad_double(xp[ind]);
        qyf = qyf+quad_double(yp[ind]); 
        qzf = qzf+quad_double(zp[ind]); 
      }
      qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      xf = E_Float(qxf);
      yf = E_Float(qyf);
      zf = E_Float(qzf);
#else
      {
      #ifdef __INTEL_COMPILER
      #pragma float_control(precise, on)
      #endif
      for (E_Int n = 1; n <= nvert; n++)
      {
        ind = (*cnl)(i,n)-1; 
        xf += xp[ind]; yf += yp[ind]; zf += zp[ind];        
      }
      xf*=inv; yf*=inv; zf*=inv;
      }
#endif
      pt[0] = xf; pt[1] = yf; pt[2] = zf;

      ind = globalKdt->getClosest(pt); // closest pt
      dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
      //d = dx*dx + dy*dy + dz*dz;
      //if (d < tol) nptr[i] = ind+1;
      //else nptr[i] = -1;
      if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
      else nptr[i] = -1;
    }
    }
  }
  else
  {
    RELEASESHAREDB(res, array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyElements: invalid type of array.");
    return NULL; 
  }
  RELEASESHAREDB(res, array, f, cnl);
  return ac;
}

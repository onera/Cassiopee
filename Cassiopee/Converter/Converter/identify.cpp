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

// Routines d'identification geometrique

# include "Connect/connect.h"
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
  PyObject* array; PyObject* hook;
  E_Float atol;  // absolute tolerance
  E_Float rtol;  // relative tolerance
  if (!PYPARSETUPLE_(args, OO_ RR_, &hook, &array, &atol, &rtol)) return NULL;

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
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, nil, njl, nkl, cn, eltType);

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
    RELEASESHAREDB(res, array, f, cn);
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

  E_Float etol = atol;  // effective tolerance
  if (!K_FUNC::fEqualZero(rtol, K_CONST::E_ZERO_MACHINE))
  {
    // Compute the maximum delta in x-, y-, z-coordinates
    // and add rtol*max(dx, dy, dz) to the absolute tolerance
    const E_Int nthreads = __NUMTHREADS__;
    E_Float dx, dy, dz;
    E_Float xmin = K_CONST::E_MAX_FLOAT; E_Float xmax = -K_CONST::E_MAX_FLOAT;
    E_Float ymin = K_CONST::E_MAX_FLOAT; E_Float ymax = -K_CONST::E_MAX_FLOAT;
    E_Float zmin = K_CONST::E_MAX_FLOAT; E_Float zmax = -K_CONST::E_MAX_FLOAT;
    E_Float* txmin = new E_Float [nthreads]; E_Float* txmax = new E_Float [nthreads];
    E_Float* tymin = new E_Float [nthreads]; E_Float* tymax = new E_Float [nthreads];
    E_Float* tzmin = new E_Float [nthreads]; E_Float* tzmax = new E_Float [nthreads];

    for (E_Int tid = 0; tid < nthreads; tid++)
    {
      txmin[tid] = K_CONST::E_MAX_FLOAT; txmax[tid] = -K_CONST::E_MAX_FLOAT;
      tymin[tid] = K_CONST::E_MAX_FLOAT; tymax[tid] = -K_CONST::E_MAX_FLOAT;
      tzmin[tid] = K_CONST::E_MAX_FLOAT; tzmax[tid] = -K_CONST::E_MAX_FLOAT;
    }

    #pragma omp parallel
    {
      const E_Int tid = __CURRENT_THREAD__;
      #pragma omp for
      for (E_Int i = 0; i < npts; i++)
      {
        txmin[tid] = K_FUNC::E_min(txmin[tid], xp[i]);
        txmax[tid] = K_FUNC::E_max(txmax[tid], xp[i]);
        tymin[tid] = K_FUNC::E_min(tymin[tid], yp[i]);
        tymax[tid] = K_FUNC::E_max(tymax[tid], yp[i]);
        tzmin[tid] = K_FUNC::E_min(tzmin[tid], zp[i]);
        tzmax[tid] = K_FUNC::E_max(tzmax[tid], zp[i]);
      }
    }

    for (E_Int tid = 0; tid < nthreads; tid++)
    {
      xmin = K_FUNC::E_min(txmin[tid], xmin);
      xmax = K_FUNC::E_max(txmax[tid], xmax);
      ymin = K_FUNC::E_min(tymin[tid], ymin);
      ymax = K_FUNC::E_max(tymax[tid], ymax);
      zmin = K_FUNC::E_min(tzmin[tid], zmin);
      zmax = K_FUNC::E_max(tzmax[tid], zmax);
    }

    dx = xmax - xmin; dy = ymax - ymin; dz = zmax - zmin;
    etol += rtol*K_FUNC::E_max(dx, K_FUNC::E_max(dy, dz));
    delete[] txmin; delete[] txmax;
    delete[] tymin; delete[] tymax;
    delete[] tzmin; delete[] tzmax;
  }

  // Remplissage
#pragma omp parallel
  {
    E_Float pt[3];
    E_Float xf, yf, zf, dx, dy, dz, dist;
    E_Int ind;

#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < npts; i++)
    {
      xf = xp[i]; yf = yp[i]; zf = zp[i];
      pt[0] = xf; pt[1] = yf; pt[2] = zf;
      ind = globalKdt->getClosest(pt); // closest pt
      dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
      dist = sqrt(dx*dx + dy*dy + dz*dz);
      if (dist < etol) nptr[i] = ind+1;
      else nptr[i] = -1;
    }
  }
  RELEASESHAREDB(res, array, f, cn);
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
  E_Float atol;  // absolute tolerance
  E_Float rtol;  // relative tolerance
  if (!PYPARSETUPLE_(args, OO_ RR_, &hook, &array, &atol, &rtol)) return NULL;

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
                    "identifyFaces: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, nil, njl, nkl, cn, eltType);

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
    RELEASESHAREDU(array, f, cn);
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
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "identifyFaces: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int nfaces = cn->getNFaces();
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);
  
  PyObject* ac = K_NUMPY::buildNumpyArray(nfaces, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Acces non universel sur les ptrs
  E_Int* ngon = cn->getNGon();
  E_Int* indPG = cn->getIndPG();

#pragma omp parallel default(shared)
  {
    E_Int v1, v2, nv, ind;
    E_Float inv, xf, yf, zf, dx, dy, dz, len, lmax, dist, etol;
    E_Float pt[3];
#ifdef QUADDOUBLE
    quad_double qxf, qyf, qzf, qdx, qdy, qdz, qlen;
    quad_double qinv;
#endif

#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < nfaces; i++)
    {
      lmax = -K_CONST::E_MAX_FLOAT;
          
      // Acces universel face i
      E_Int* face = cn->getFace(i, nv, ngon, indPG);

#ifdef QUADDOUBLE
      qxf = 0.; qyf = 0.; qzf = 0.;
      for (E_Int j = 0; j < nv; j++)
      {
        v1 = face[j]-1;
        qxf += quad_double(xp[v1]); 
        qyf += quad_double(yp[v1]); 
        qzf += quad_double(zp[v1]);

        // local max edge length
        if (j == nv-1) v2 = face[0]-1;
        else v2 = face[j+1]-1;
        qdx = quad_double(xp[v2]) - quad_double(xp[v1]);
        qdy = quad_double(yp[v2]) - quad_double(yp[v1]);
        qdz = quad_double(zp[v2]) - quad_double(zp[v1]);
        qlen = sqrt(qdx*qdx + qdy*qdy + qdz*qdz);
        len = E_Float(qlen);
        lmax = K_FUNC::E_max(lmax, len);
      }
      qinv = quad_double(nv); qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
#else
      {
        xf = 0.; yf = 0.; zf = 0.;
        #ifdef __INTEL_COMPILER
        #pragma float_control(precise, on)
        #endif
        for (E_Int j = 0; j < nv; j++)
        {
          v1 = face[j]-1;
          xf += xp[v1]; yf += yp[v1]; zf += zp[v1];

          // local max edge length
          if (j == nv-1) v2 = face[0]-1;
          else v2 = face[j+1]-1;
          dx = xp[v2]-xp[v1]; dy = yp[v2]-yp[v1]; dz = zp[v2]-zp[v1];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
        }
        inv = 1./E_Float(nv); xf *= inv; yf *= inv; zf *= inv;
      }
#endif
      pt[0] = xf; pt[1] = yf; pt[2] = zf;
      ind = globalKdt->getClosest(pt); // closest pt
      dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
      dist = sqrt(dx*dx + dy*dy + dz*dz);
      // local effective tolerance
      etol = atol + rtol*lmax;
      if (dist < etol) nptr[i] = ind+1; 
      else nptr[i] = -1;
    }
  }

  RELEASESHAREDU(array, f, cn);
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
  E_Float atol;  // absolute tolerance
  E_Float rtol;  // relative tolerance
  if (!PYPARSETUPLE_(args, OO_ RR_, &hook, &array, &atol, &rtol)) return NULL;

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
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString, 
                               f, nil, njl, nkl, cn, eltType);

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
    RELEASESHAREDU(array, f, cn);
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
      E_Float xf, yf, zf, dx, dy, dz, dist, lmax, len, etol;
      E_Int ind, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
      E_Float pt[3];
      E_Int ic = 0;

#pragma omp for schedule(dynamic) collapse(3)
      for (E_Int k = 0; k < nke; k++)
      for (E_Int j = 0; j < nje; j++)
      for (E_Int i = 0; i < nie; i++)
      {
        lmax = -K_CONST::E_MAX_FLOAT;
        etol = atol;

        ic = i+j*nie+k*nie*nje;
        xf = 0.; yf = 0.; zf = 0.;
        ind1 = i+j*nil+k*nijl;
        xf += xp[ind1]; yf += yp[ind1]; zf += zp[ind1];
        ind2 = i+inci+j*nil+k*nijl;
        xf += xp[ind2]; yf += yp[ind2]; zf += zp[ind2];
        ind3 = i+(j+incj)*nil+k*nijl;
        xf += xp[ind3]; yf += yp[ind3]; zf += zp[ind3];
        ind4 = i+inci+(j+incj)*nil+k*nijl;
        xf += xp[ind4]; yf += yp[ind4]; zf += zp[ind4];
        ind5 = i+j*nil+(k+inck)*nijl;
        xf += xp[ind5]; yf += yp[ind5]; zf += zp[ind5];
        ind6 = i+inci+j*nil+(k+inck)*nijl;
        xf += xp[ind6]; yf += yp[ind6]; zf += zp[ind6];
        ind7 = i+(j+incj)*nil+(k+inck)*nijl;
        xf += xp[ind7]; yf += yp[ind7]; zf += zp[ind7];
        ind8 = i+inci+(j+incj)*nil+(k+inck)*nijl;
        xf += xp[ind8]; yf += yp[ind8]; zf += zp[ind8];

        xf = xf/8.; yf = yf/8.; zf = zf/8.;
        pt[0] = xf; pt[1] = yf; pt[2] = zf;

        if (!K_FUNC::fEqualZero(rtol, K_CONST::E_ZERO_MACHINE))
        {
          // Edge (ind1, ind2)
          dx = xp[ind2]-xp[ind1]; dy = yp[ind2]-yp[ind1]; dz = zp[ind2]-zp[ind1];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind1, ind3)
          dx = xp[ind3]-xp[ind1]; dy = yp[ind3]-yp[ind1]; dz = zp[ind3]-zp[ind1];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind1, ind5)
          dx = xp[ind5]-xp[ind1]; dy = yp[ind5]-yp[ind1]; dz = zp[ind5]-zp[ind1];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind2, ind4)
          dx = xp[ind4]-xp[ind2]; dy = yp[ind4]-yp[ind2]; dz = zp[ind4]-zp[ind2];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind2, ind6)
          dx = xp[ind6]-xp[ind2]; dy = yp[ind6]-yp[ind2]; dz = zp[ind6]-zp[ind2];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind3, ind4)
          dx = xp[ind4]-xp[ind3]; dy = yp[ind4]-yp[ind3]; dz = zp[ind4]-zp[ind3];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind3, ind7)
          dx = xp[ind7]-xp[ind3]; dy = yp[ind7]-yp[ind3]; dz = zp[ind7]-zp[ind3];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind4, ind8)
          dx = xp[ind8]-xp[ind4]; dy = yp[ind8]-yp[ind4]; dz = zp[ind8]-zp[ind4];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind5, ind6)
          dx = xp[ind6]-xp[ind5]; dy = yp[ind6]-yp[ind5]; dz = zp[ind6]-zp[ind5];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind5, ind7)
          dx = xp[ind7]-xp[ind5]; dy = yp[ind7]-yp[ind5]; dz = zp[ind7]-zp[ind5];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind6, ind8)
          dx = xp[ind8]-xp[ind6]; dy = yp[ind8]-yp[ind6]; dz = zp[ind8]-zp[ind6];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Edge (ind7, ind8)
          dx = xp[ind8]-xp[ind7]; dy = yp[ind8]-yp[ind7]; dz = zp[ind8]-zp[ind7];
          len = sqrt(dx*dx + dy*dy + dz*dz);
          lmax = K_FUNC::E_max(lmax, len);
          // Effective tolerance
          etol += rtol*lmax;
        }

        ind = globalKdt->getClosest(pt); // closest pt
        dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
        dist = sqrt(dx*dx + dy*dy + dz*dz);
        if (dist < etol) nptr[ic] = ind+1;
        else nptr[ic] = -1;
      }
    }
  }
  else if (strcmp(eltType, "NGON") == 0)
  {  
    // Cree le numpy de sortie
    E_Int nelts = cn->getNElts();
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);

    ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

    // Acces non universel sur les ptrs
    E_Int* ngon = cn->getNGon();
    E_Int* nface = cn->getNFace();
    E_Int* indPG = cn->getIndPG();
    E_Int* indPH = cn->getIndPH();
    
#pragma omp parallel default(shared)
    {
      E_Float pt[3];
      E_Int nf, c, nv, ind, v1, v2;
      E_Float xf, yf, zf, inv, dx, dy, dz;
      E_Float lmax, len, dist, etol;

#ifdef QUADDOUBLE
      quad_double qxf, qyf, qzf, qdx, qdy, qdz, qlen;
      quad_double qinv;
#endif

#pragma omp for schedule(dynamic)
      for (E_Int i = 0; i < nelts; i++)
      {
        lmax = -K_CONST::E_MAX_FLOAT;

        // Acces universel element i
        E_Int* elem = cn->getElt(i, nf, nface, indPH);
        xf = 0.; yf = 0.; zf = 0.; c = 0;

#ifdef QUADDOUBLE
        qxf = 0.; qyf = 0.; qzf = 0.;
        for (E_Int n = 0; n < nf; n++)
        { 
          E_Int* face = cn->getFace(elem[n]-1, nv, ngon, indPG);
          for (E_Int j = 0; j < nv; j++)
          {
            v1 = face[j]-1;
            qxf += quad_double(xp[v1]);
            qyf += quad_double(yp[v1]);
            qzf += quad_double(zp[v1]);
            c++;

            // local max edge length
            if (j == nv-1) v2 = face[0]-1;
            else v2 = face[j+1]-1;
            qdx = quad_double(xp[v2]) - quad_double(xp[v1]);
            qdy = quad_double(yp[v2]) - quad_double(yp[v1]);
            qdz = quad_double(zp[v2]) - quad_double(zp[v1]);
            qlen = sqrt(qdx*qdx + qdy*qdy + qdz*qdz);
            len = E_Float(qlen);
            lmax = K_FUNC::E_max(lmax, len);
          }
        }
        qinv = quad_double(c); qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
        xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
#else
        {
          #ifdef __INTEL_COMPILER
          #pragma float_control(precise, on)
          #endif
          for (E_Int n = 0; n < nf; n++)
          {
            E_Int* face = cn->getFace(elem[n]-1, nv, ngon, indPG);
            for (E_Int j = 0; j < nv; j++)
            {
              v1 = face[j]-1;
              xf += xp[v1]; yf += yp[v1]; zf += zp[v1]; c++;

              // local max edge length
              if (j == nv-1) v2 = face[0]-1;
              else v2 = face[j+1]-1;
              dx = xp[v2]-xp[v1]; dy = yp[v2]-yp[v1]; dz = zp[v2]-zp[v1];
              len = sqrt(dx*dx + dy*dy + dz*dz);
              lmax = K_FUNC::E_max(lmax, len);
            }
          }
          inv = 1./E_Float(c); xf *= inv; yf *= inv; zf *= inv;
        }
#endif

        pt[0] = xf; pt[1] = yf; pt[2] = zf;
        ind = globalKdt->getClosest(pt); // closest pt
        dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
        dist = sqrt(dx*dx + dy*dy + dz*dz);
        // local effective tolerance
        etol = atol + rtol*lmax;
        if (dist < etol) nptr[i] = ind+1;
        else nptr[i] = -1;
      }
    }
  }
  else
  {
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);

    // Acces universel sur BE/ME
    E_Int nc = cn->getNConnect();
    // Acces universel aux eltTypes
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);
    std::vector<E_Int> nepc(nc+1, 0);

    // Number of faces per element for each connectivity
    std::vector<E_Int> nfpe;
    E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, true);
    if (ierr != 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "identifyElements: cannot determine the number of faces "
                      "per element.");
      return NULL;
    }

    // Boucle sur toutes les connectivites une premiere fois pour savoir si
    // elles sont valides et calculer le nombre total d'elements - permet
    // d'allouer ac
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      char* eltTypConn = eltTypes[ic];
      // Check that this connectivity is valid
      if (not(strcmp(eltTypConn,"BAR")==0 || strcmp(eltTypConn,"TRI")==0 || 
              strcmp(eltTypConn,"QUAD")==0 || strcmp(eltTypConn,"TETRA")==0 || 
              strcmp(eltTypConn,"HEXA")==0 || strcmp(eltTypConn,"PENTA")==0 ||
              strcmp(eltTypConn,"PYRA")==0))
      {
        RELEASESHAREDB(res, array, f, cn);
        PyErr_SetString(PyExc_TypeError,
                        "identifyElements: invalid type of array.");
        return NULL;
      }
      nepc[ic+1] += nepc[ic] + cm.getSize();
    }

    ac = K_NUMPY::buildNumpyArray(nepc[nc], 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

    // Boucle sur toutes les connectivites pour remplir ac
#pragma omp parallel default(shared)
    {
      E_Int nvpf;
      std::vector<std::vector<E_Int> > facets;

      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cn->getConnect(ic));
        E_Int nelts = cm.getSize();
        E_Int offset = nepc[ic];
        E_Int nvpe = cm.getNfld();
        E_Float inv = 1./E_Float(nvpe);
        K_CONNECT::getEVFacets(facets, eltTypes[ic], false);
  
#ifdef QUADDOUBLE
        quad_double qinv = quad_double(nvpe);
#endif

        E_Float pt[3];
        E_Float xf, yf, zf, dx, dy, dz, lmax, len, dist, etol;
        E_Int ind, v1, v2;
#ifdef QUADDOUBLE
        quad_double qxf, qyf, qzf, qdx, qdy, qdz, qlen;
#endif

#pragma omp for schedule(dynamic)
        for (E_Int i = 0; i < nelts; i++)
        {
          lmax = -K_CONST::E_MAX_FLOAT;

#ifdef QUADDOUBLE
          qxf = 0.; qyf = 0.; qzf = 0.;
          for (E_Int n = 1; n <= nvpe; n++)
          {
            v1 = cm(i,n)-1;
            qxf += quad_double(xp[v1]);
            qyf += quad_double(yp[v1]); 
            qzf += quad_double(zp[v1]); 
          }
          qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
          xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);

          // Loop over each facet of this element
          for (E_Int f = 0; f < nfpe[ic]; f++)
          {
            // Number of vertices per face
            nvpf = facets[f].size();
            for (E_Int j = 0; j < nvpf; j++)
            {
              v1 = cm(i, facets[f][j]) - 1;
              if (j == nvpf-1) v2 = cm(i, facets[f][0]) - 1;
              else v2 = cm(i, facets[f][j+1]) - 1;
              qdx = quad_double(xp[v2]) - quad_double(xp[v1]);
              qdy = quad_double(yp[v2]) - quad_double(yp[v1]);
              qdz = quad_double(zp[v2]) - quad_double(zp[v1]);
              qlen = sqrt(qdx*qdx + qdy*qdy + qdz*qdz);
              len = E_Float(qlen);
              lmax = K_FUNC::E_max(lmax, len);
            }
          }
#else
          xf = 0.; yf = 0.; zf = 0.;
          {
            #ifdef __INTEL_COMPILER
            #pragma float_control(precise, on)
            #endif
            for (E_Int n = 1; n <= nvpe; n++)
            {
              v1 = cm(i,n)-1; 
              xf += xp[v1]; yf += yp[v1]; zf += zp[v1];
            }
            xf *= inv; yf *= inv; zf *= inv;
          }

          // Loop over each facet of this element
          for (E_Int f = 0; f < nfpe[ic]; f++)
          {
            // Number of vertices per face
            nvpf = facets[f].size();
            for (E_Int j = 0; j < nvpf; j++)
            {
              v1 = cm(i, facets[f][j]) - 1;
              if (j == nvpf-1) v2 = cm(i, facets[f][0]) - 1;
              else v2 = cm(i, facets[f][j+1]) - 1;
              dx = xp[v2]-xp[v1]; dy = yp[v2]-yp[v1]; dz = zp[v2]-zp[v1];
              len = sqrt(dx*dx + dy*dy + dz*dz);
              lmax = K_FUNC::E_max(lmax, len);
            }
          }
#endif
          pt[0] = xf; pt[1] = yf; pt[2] = zf;
          ind = globalKdt->getClosest(pt); // closest pt
          dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
          dist = sqrt(dx*dx + dy*dy + dz*dz);
          // local effective tolerance
          etol = atol + rtol*lmax;
          if (dist < etol) nptr[offset+i] = ind+1;
          else nptr[offset+i] = -1;
        }
      }
    }

    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  }

  RELEASESHAREDB(res, array, f, cn);
  return ac;
}

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
/* Trouve pour les noeuds de a le point le plus proche parmi les points
   stockes dans un KdTree (hook).
   Retourne la liste des noeuds dans la numerotation du KDT. */
// ============================================================================
PyObject* K_CONVERTER::nearestNodes(PyObject* self, PyObject* args)
{ 
  PyObject* array; PyObject* hook;
  if (!PyArg_ParseTuple(args, "OO", &hook, &array)) return NULL;

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
                    "nearestNodes: this function requires a identify KDT hook.");
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
                    "nearestNodes: array is invalid.");
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
                    "nearestNodes: array must have coordinates.");
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
  E_Int* nptr1 = K_NUMPY::getNumpyPtrI(ac);
  PyObject* dist = K_NUMPY::buildNumpyArray(npts, 1, 0);
  E_Float* nptr2 = K_NUMPY::getNumpyPtrF(dist);

  // Remplissage
  E_Int ind;
  E_Float xf, yf, zf, d;
  E_Float pt[3];
  
  for (E_Int i = 0; i < npts; i++)
  {
    xf = xp[i]; yf = yp[i]; zf = zp[i];
    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    ind = globalKdt->getClosest(pt); // closest pt
    d = (xt[ind]-xf)*(xt[ind]-xf)+(yt[ind]-yf)*(yt[ind]-yf)+(zt[ind]-zf)*(zt[ind]-zf);
    nptr1[i] = ind+1;
    nptr2[i] = sqrt(d);
  }
  RELEASESHAREDB(res,array, f, cnl);
  PyObject* tpl;
  tpl = Py_BuildValue("[OO]",ac,dist);
  return tpl;
}

// ============================================================================
/* Trouve pour les  centres des faces de a le point le plus proche parmi 
   les points stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne la liste des faces. */
// ============================================================================
PyObject* K_CONVERTER::nearestFaces(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* hook;
  if (!PyArg_ParseTuple(args, "OO", &hook, &array)) return NULL;

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
                    "nearestFaces: this function requires a identify KDT hook.");
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
                    "nearestFaces: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "nearestFaces: array must be a NGON.");
    return NULL; 
  }
  if (strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "nearestFaces: array must be a NGON.");
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
                    "nearestFaces: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;

  // Cree le numpy de sortie
  E_Int* ptr = cnl->begin();
  E_Int nfaces = ptr[0];
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);

  PyObject* ac = K_NUMPY::buildNumpyArray(nfaces, 1, 1);
  E_Int* nptr1 = K_NUMPY::getNumpyPtrI(ac);
  PyObject* dist = K_NUMPY::buildNumpyArray(nfaces, 1, 0);
  E_Float* nptr2 = K_NUMPY::getNumpyPtrF(dist);


  FldArrayI posFace; K_CONNECT::getPosFaces(*cnl, posFace);

# pragma omp parallel for default(shared)
  for (E_Int i = 0; i < nfaces; i++)
  {
    E_Int posf = posFace[i];
    E_Int* ptrFace = &ptr[posf];
    E_Int nv = ptrFace[0];
    E_Float xf=0.,yf=0.,zf=0.;
    quad_double  qxf, qyf, qzf;
    for (E_Int n = 1; n <= nv; n++)
    {
      E_Int indv = ptrFace[n]-1;

#ifdef QUADDOUBLE
      qxf = qxf+quad_double(xp[indv]); 
      qyf = qyf+quad_double(yp[indv]); 
      qzf = qzf+quad_double(zp[indv]); 
    }//loop on vertices

    E_Float qinv = quad_double(nv); 
    qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
    xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
# else
    {
      #ifdef __INTEL_COMPILER
      #pragma float_control(precise, on)
      #endif
      xf += xp[indv]; yf += yp[indv]; zf += zp[indv];
    }
    }//loop on vertices
    E_Float inv = 1./E_Float(nv); xf *= inv; yf *= inv; zf *= inv;
#endif

    E_Float pt[3];  pt[0] = xf; pt[1] = yf; pt[2] = zf;
    E_Int ind = globalKdt->getClosest(pt); // closest pt
    E_Float d = (xt[ind]-xf)*(xt[ind]-xf)+(yt[ind]-yf)*(yt[ind]-yf)+(zt[ind]-zf)*(zt[ind]-zf);
    nptr1[i] = ind+1;
    nptr2[i] = sqrt(d);
  }

  RELEASESHAREDU(array, f, cnl);
  PyObject* tpl;
  tpl = Py_BuildValue("[OO]",ac,dist);
  return tpl;
}

// ============================================================================
/* Trouve pour les centres des elements de a le plus proche parmi les points 
   stockes dans un KdTree (hook).
   IN: a: NGON
   Retourne la liste des elements. */
// ============================================================================
PyObject* K_CONVERTER::nearestElements(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* hook;
  if (!PyArg_ParseTuple(args, "OO", &hook, &array)) return NULL;

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
                    "nearestElts: this function requires a identify KDT hook.");
    return NULL;
  }
  FldArrayF* centers = (FldArrayF*)packet[1];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    (K_SEARCH::KdTree<FldArrayF>*) packet[3];

  // Recupere l'array a identifier
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "nearestElts: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "nearestElts: array must be a NGON.");
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
                    "nearestElts: array must have coordinates.");
    return NULL; 
  }
  posx++; posy++; posz++;
  if (strcmp(eltType, "NGON") == 0)
  {  
    // Cree le numpy de sortie
    E_Int sizef = (*cnl)[1];
    E_Int nelts = (*cnl)[sizef+2];
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);

    PyObject* ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr1 = K_NUMPY::getNumpyPtrI(ac);
    PyObject* dist = K_NUMPY::buildNumpyArray(nelts, 1, 0);
    E_Float* nptr2 = K_NUMPY::getNumpyPtrF(dist);

    // Remplissage
    E_Int* ptr = cnl->begin();
    E_Int* ptrElt = cnl->begin(); ptrElt += ptrElt[1]+4;
    FldArrayI posFace; K_CONNECT::getPosFaces(*cnl, posFace);
    FldArrayI posElt; K_CONNECT::getPosElts(*cnl, posElt);

  # pragma omp parallel for default(shared)
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int pose = posElt[i];
      E_Int* ptrElt = &ptr[pose];
      E_Int nf = ptrElt[0]; 
      E_Float xf = 0., yf = 0., zf = 0.; 
      E_Int c = 0;
      quad_double qxf, qyf, qzf;

      for (E_Int n = 1; n <= nf; n++)
      { 
        E_Int ind = ptrElt[n]-1;
        E_Int pos = posFace[ind];
        E_Int* ptrFace = &ptr[pos];
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
      quad_double qinv = quad_double(c); qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
  #else
      {
      #ifdef __INTEL_COMPILER
      #pragma float_control(precise, on)
      #endif
      for (E_Int p = 1; p <= nv; p++)
      { 
        ind = ptrFace[p]-1; 
        xf += xp[ind]; yf += yp[ind]; zf += zp[ind]; c++; 
      }
      }
    }
    E_Float inv = 1./E_Float(c); xf *= inv; yf *= inv; zf *= inv;
  #endif

    E_Float pt[3];
    pt[0] = xf; pt[1] = yf; pt[2] = zf;
    E_Int ind = globalKdt->getClosest(pt); // closest pt
    E_Float d = (xt[ind]-xf)*(xt[ind]-xf)+(yt[ind]-yf)*(yt[ind]-yf)+(zt[ind]-zf)*(zt[ind]-zf);
    nptr1[i] = ind+1; nptr2[i] = sqrt(d);
  }
  RELEASESHAREDU(array, f, cnl);
  return Py_BuildValue("[OO]",ac,dist);
  }//NGON
  else if (strcmp(eltType,"BAR")==0 || strcmp(eltType,"TRI")==0 || 
           strcmp(eltType,"QUAD")==0 || strcmp(eltType,"TETRA")==0 || 
           strcmp(eltType,"HEXA")==0 || strcmp(eltType,"PENTA")==0 ||
           strcmp(eltType,"PYRA")==0)
  {
    E_Int nelts = cnl->getSize();
    E_Int nvert = cnl->getNfld();
    E_Float* xp = f->begin(posx);
    E_Float* yp = f->begin(posy);
    E_Float* zp = f->begin(posz);
    E_Float* xt = centers->begin(1);
    E_Float* yt = centers->begin(2);
    E_Float* zt = centers->begin(3);

    PyObject* ac = K_NUMPY::buildNumpyArray(nelts, 1, 1);
    E_Int* nptr1 = K_NUMPY::getNumpyPtrI(ac);
    PyObject* dist = K_NUMPY::buildNumpyArray(nelts, 1, 0);
    E_Float* nptr2 = K_NUMPY::getNumpyPtrF(dist);
    E_Float inv = 1./E_Float(nvert);
    quad_double qinv = quad_double(nvert);
    E_Float pt[3];

#pragma omp parallel for default(shared) private(pt) schedule(dynamic)
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Float dx,dy,dz;
      E_Int ind;
      E_Float xf=0.,yf=0.,zf=0.;
      quad_double qxf, qyf, qzf;

#ifdef QUADDOUBLE
      for (E_Int n = 1; n <= nvert; n++)
      {
        ind = (*cnl)(i,n)-1; 
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

      E_Float d = dx*dx+dy*dy*dz*dz;
      nptr1[i] = ind+1; nptr2[i] = sqrt(d);
    }
    RELEASESHAREDU(array, f, cnl);
    return Py_BuildValue("[OO]",ac,dist);
  }//EB
  else
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "nearestElts: invalid element type.");
    return NULL; 
  }
}


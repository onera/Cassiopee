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
  if (!PYPARSETUPLE_(args, OO_ R_, &hook, &array, &tol)) return NULL;

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
  res = K_ARRAY::getFromArray3(array, varString, 
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
  if (!PYPARSETUPLE_(args, OO_ R_, &hook, &array, &tol)) return NULL;

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
  res = K_ARRAY::getFromArray3(array, varString, 
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
  E_Int nfaces = cnl->getNFaces();
  E_Float* xp = f->begin(posx);
  E_Float* yp = f->begin(posy);
  E_Float* zp = f->begin(posz);
  E_Float* xt = centers->begin(1);
  E_Float* yt = centers->begin(2);
  E_Float* zt = centers->begin(3);
  
  PyObject* ac = K_NUMPY::buildNumpyArray(nfaces, 1, 1);
  E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

  // Acces non universel sur les ptrs
  E_Int* ngon = cnl->getNGon();
  E_Int* indPG = cnl->getIndPG();

#pragma omp parallel default(shared)
  {
    E_Int ind, nv, indk;
    E_Float inv, xf, yf, zf, dx, dy, dz;
    E_Float pt[3];
#ifdef QUADDOUBLE
    quad_double qxf, qyf, qzf;
    quad_double qinv;
#endif

#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < nfaces; i++)
    {
      // Acces universel face i
      E_Int* face = cnl->getFace(i, nv, ngon, indPG);

#ifdef QUADDOUBLE
      qxf = 0.; qyf = 0.; qzf = 0.;
      for (E_Int n = 0; n < nv; n++)
      {
        ind = face[n]-1;
        qxf = qxf+quad_double(xp[ind]); 
        qyf = qyf+quad_double(yp[ind]); 
        qzf = qzf+quad_double(zp[ind]); 
      }
      qinv = quad_double(nv); qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
      xf = E_Float(qxf); yf = E_Float(qyf); zf = E_Float(qzf);
#else
      {
        xf = 0.; yf = 0.; zf = 0.;
        #ifdef __INTEL_COMPILER
        #pragma float_control(precise, on)
        #endif
        for (E_Int n = 0; n < nv; n++)
        {
          ind = face[n]-1;
          xf += xp[ind]; yf += yp[ind]; zf += zp[ind];
        }
        inv = 1./E_Float(nv); xf *= inv; yf *= inv; zf *= inv;
      }
#endif
      pt[0] = xf; pt[1] = yf; pt[2] = zf;
      indk = globalKdt->getClosest(pt); // closest pt
      dx = xt[indk]-xf; dy = yt[indk]-yf; dz = zt[indk]-zf;
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
  if (!PYPARSETUPLE_(args, OO_ R_, &hook, &array, &tol)) return NULL;

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
  res = K_ARRAY::getFromArray3(array, varString, 
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

    // Acces non universel sur les ptrs
    E_Int* ngon = cnl->getNGon();
    E_Int* nface = cnl->getNFace();
    E_Int* indPG = cnl->getIndPG();
    E_Int* indPH = cnl->getIndPH();
    
#pragma omp parallel default(shared)
    {
      E_Float pt[3];
      E_Int nf, c, nv, ind, indp;
      E_Float xf, yf, zf, inv, dx, dy, dz;
#ifdef QUADDOUBLE
      quad_double qxf, qyf, qzf;
      quad_double qinv;
#endif

#pragma omp for schedule(dynamic)
      for (E_Int i = 0; i < nelts; i++)
      {
        // Acces universel element i
        E_Int* elem = cnl->getElt(i, nf, nface, indPH);
        xf = 0.; yf = 0.; zf = 0.; c = 0;

#ifdef QUADDOUBLE
        qxf = 0.; qyf = 0.; qzf = 0.;
        for (E_Int n = 0; n < nf; n++)
        { 
          // Acces universel face elem[n]-1
          E_Int* face = cnl->getFace(elem[n]-1, nv, ngon, indPG);
          for (E_Int p = 0; p < nv; p++)
          {
            indp = face[p]-1;
            qxf = qxf+quad_double(xp[indp]); 
            qyf = qyf+quad_double(yp[indp]);
            qzf = qzf+quad_double(zp[indp]); 
            c++;
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
            // Acces universel face elem[n]-1
            E_Int* face = cnl->getFace(elem[n]-1, nv, ngon, indPG);
            for (E_Int p = 0; p < nv; p++)
            {
              indp = face[p]-1; xf += xp[indp]; yf += yp[indp]; zf += zp[indp]; c++;
            }
          }
          inv = 1./E_Float(c); xf *= inv; yf *= inv; zf *= inv;
        }
#endif
        pt[0] = xf; pt[1] = yf; pt[2] = zf;
        ind = globalKdt->getClosest(pt); // closest pt
        dx = xt[ind]-xf; dy = yt[ind]-yf; dz = zt[ind]-zf;
        if (K_FUNC::E_abs(dx) < tol && K_FUNC::E_abs(dy) < tol && K_FUNC::E_abs(dz) < tol) nptr[i] = ind+1;
        else nptr[i] = -1;
      }
    }
  }
  else
  {
    // Acces universel sur BE/ME
    E_Int nc = cnl->getNConnect();
    // Acces universel aux eltTypes
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);
    E_Int nvert, nelts, ntotelts = 0;

    // Boucle sur toutes les connectivites une premiere fois pour savoir si
    // elles sont valides et calculer le nombre total d'elements - permet
    // d'allouer ac
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cnl->getConnect(ic));
      char* eltTypConn = eltTypes[ic];
      // Check that this connectivity is valid
      if (not(strcmp(eltTypConn,"BAR")==0 || strcmp(eltTypConn,"TRI")==0 || 
          strcmp(eltTypConn,"QUAD")==0 || strcmp(eltTypConn,"TETRA")==0 || 
          strcmp(eltTypConn,"HEXA")==0 || strcmp(eltTypConn,"PENTA")==0 ||
          strcmp(eltTypConn,"PYRA")==0))
      {
        RELEASESHAREDB(res, array, f, cnl);
        PyErr_SetString(PyExc_TypeError, 
                        "identifyElements: invalid type of array.");
        return NULL;
      }
      ntotelts += cm.getSize();
    }

    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

    ac = K_NUMPY::buildNumpyArray(ntotelts, 1, 1);
    E_Int* nptr = K_NUMPY::getNumpyPtrI(ac);

    // Boucle sur toutes les connectivites pour remplir ac
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cnl->getConnect(ic));
      nelts = cm.getSize(); // Number of elements of this connectivity
      nvert = cm.getNfld(); // Nombre de points par elements de cette connectivite
      E_Float inv = 1./E_Float(nvert);
#ifdef QUADDOUBLE
      quad_double qinv = quad_double(nvert);
#endif
      E_Float* xt = centers->begin(1);
      E_Float* yt = centers->begin(2);
      E_Float* zt = centers->begin(3);

#pragma omp parallel default(shared)
      {
        E_Float pt[3];
        E_Float xf,yf,zf,dx,dy,dz;
        E_Int ind;
#ifdef QUADDOUBLE
        quad_double qxf, qyf, qzf;
#endif

#pragma omp for schedule(dynamic)
        for (E_Int i = 0; i < nelts; i++)
        {
#ifdef QUADDOUBLE
          qxf = 0.; qyf =0.; qzf = 0.;
          for (E_Int n = 1; n <= nvert; n++)
          {
            ind = cm(i,n)-1;
            qxf = qxf+quad_double(xp[ind]);
            qyf = qyf+quad_double(yp[ind]); 
            qzf = qzf+quad_double(zp[ind]); 
          }
          qxf = qxf/qinv; qyf = qyf/qinv; qzf = qzf/qinv;
          xf = E_Float(qxf);
          yf = E_Float(qyf);
          zf = E_Float(qzf);
#else
          xf = 0.; yf = 0.; zf = 0.;
          {
            #ifdef __INTEL_COMPILER
            #pragma float_control(precise, on)
            #endif
            for (E_Int n = 1; n <= nvert; n++)
            {
              ind = cm(i,n)-1; 
              xf += xp[ind]; yf += yp[ind]; zf += zp[ind];        
            }
            xf *= inv; yf *= inv; zf *= inv;
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
  }
  RELEASESHAREDB(res, array, f, cnl);
  return ac;
}

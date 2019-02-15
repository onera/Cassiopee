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

// Global hook on a kdtree - option: indirection for zones can be returned also
 
# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Enregistre les centres des faces de a dans un KdTree 
   hook type=100 
   IN: a: array NGON */
// ============================================================================
PyObject* K_CONVERTER::registerAllFaces(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Int extended;
  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arrays, &extended)) return NULL;

  E_Int extendedI = E_Int(extended);
  
  // Extract infos from arrays
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = true;
  E_Boolean skipUnstructured = false; 
  E_Boolean skipDiffVars = false;
  vector<E_Int> rest;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, rest, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  if (isOk == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createGlobalHook: arrays is not valid.");
    return NULL;
  }
  E_Int nzones = unstrF.size();
  for (E_Int v = 0; v < nzones; v++)
  {
    if (strcmp(eltTypet[v], "NGON") != 0)   
    {
      for (E_Int nos = 0; nos < nzones; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      PyErr_SetString(PyExc_TypeError, 
                    "createGlobalHook: array must be a NGON.");
      return NULL; 
    }
  }
  E_Int posx1, posy1, posz1;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[noz]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[noz]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[noz]); posz1++;
    posxu.push_back(posx1); posyu.push_back(posy1); poszu.push_back(posz1); 
  }

  // Parcours les faces, calcule les centres, les enregistre dans le KdTree
  E_Int nfacesTot=0;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    FldArrayI& cnl = *cnt[noz];
    nfacesTot+= cnl[0];
  }
  
  FldArrayF* centers = new FldArrayF(nfacesTot, 3);
  E_Float* cx = centers->begin(1);
  E_Float* cy = centers->begin(2);
  E_Float* cz = centers->begin(3);
  E_Int nv, ind;
  E_Float xf, yf, zf, inv;
  E_Int ig = 0;
  PyObject* indirZones=NULL;

  if (extendedI == 0)
  {
    for (E_Int noz = 0; noz < nzones; noz++)
    {
      E_Float* xp = unstrF[noz]->begin(posxu[noz]);
      E_Float* yp = unstrF[noz]->begin(posyu[noz]);
      E_Float* zp = unstrF[noz]->begin(poszu[noz]);        
      E_Int* ptr = cnt[noz]->begin(); ptr += 2;
      E_Int nfaces = (*cnt[noz])[0];
      for (E_Int i = 0; i < nfaces; i++)
      {
        nv = ptr[0]; 
        xf = 0.; yf = 0.; zf = 0.;
        for (E_Int n = 1; n <= nv; n++)
        { ind = ptr[n]-1; xf += xp[ind]; yf += yp[ind]; zf += zp[ind]; }
        inv = 1./nv; xf *= inv; yf *= inv; zf *= inv;
        cx[ig] = xf; cy[ig] = yf; cz[ig] = zf;
        ptr += nv+1; ig++; 
      }
    }
  }
  else 
  {
    FldArrayI* numZones = new FldArrayI(nfacesTot);    
    E_Int* ptrNumZones = numZones->begin();

    for (E_Int noz = 0; noz < nzones; noz++)
    {
      E_Float* xp = unstrF[noz]->begin(posxu[noz]);
      E_Float* yp = unstrF[noz]->begin(posyu[noz]);
      E_Float* zp = unstrF[noz]->begin(poszu[noz]);        
      E_Int* ptr = cnt[noz]->begin(); ptr += 2;
      E_Int nfaces = (*cnt[noz])[0];
      for (E_Int i = 0; i < nfaces; i++)
      {
        nv = ptr[0]; 
        xf = 0.; yf = 0.; zf = 0.;
        for (E_Int n = 1; n <= nv; n++)
        { ind = ptr[n]-1; xf += xp[ind]; yf += yp[ind]; zf += zp[ind]; }
        inv = 1./nv; xf *= inv; yf *= inv; zf *= inv;
        cx[ig] = xf; cy[ig] = yf; cz[ig] = zf;
        ptrNumZones[ig] = noz;
        ptr += nv+1; ig++; 
      }
    }
    indirZones = K_NUMPY::buildNumpyArray(*numZones);
  }
  ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*centers, 1, 2, 3); // ref sur centers
  K_SEARCH::KdTree<FldArrayF>* globalKdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc); // ref sur coordAcc

  PyObject* hook;
  E_Int* type = new E_Int [1]; type[0] = 100;
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
  /* Nettoyage et sortie */
  for (E_Int noz = 0; noz < nzones; noz++)
    RELEASESHAREDU(objut[noz], unstrF[noz], cnt[noz]);
  
  if (extendedI == 0) return hook;
  else 
  {                   
    PyObject* res = PyList_New(0);
    PyList_Append(res, hook);
    PyList_Append(res, indirZones); Py_DECREF(indirZones);
    return res;
  }
}

// ============================================================================
/* Enregistre les noeuds a dans un KdTree 
   hook type=102
   IN: a: tout type d'array avec coordonnees */
// ============================================================================
PyObject* K_CONVERTER::registerAllNodes(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Int extended;
#ifdef E_DOUBLEINT 
  if (!PyArg_ParseTuple(args, "Ol", &arrays, &extended)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oi", &arrays, &extended)) return NULL;
#endif
  if (PyList_Check(arrays) == false)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createGlobalHook: arg must be a list of arrays.");
    return NULL;
  }

  // Extract infos from arrays
  vector<E_Int> resl;  vector<char*> varString;
  vector<FldArrayF*> fields;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  E_Boolean skipNoCoord = true; E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false; E_Boolean skipDiffVars = false;
  E_Int isOk = K_ARRAY::getFromArrays(arrays, resl, varString, fields, a2, a3, a4, objs,  
                                      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = resl.size();
  if (isOk == -1) 
  {
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);  
    PyErr_SetString(PyExc_TypeError, 
                    "createGlobalHook: arrays is not valid.");
    return NULL;
  }

  E_Int posx1, posy1, posz1;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(varString[noz]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(varString[noz]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(varString[noz]); posz1++;
    posxt.push_back(posx1); posyt.push_back(posy1); poszt.push_back(posz1); 
  }

  E_Int nptsTot=0;
  for (E_Int noz = 0; noz < nzones; noz++) nptsTot += fields[noz]->getSize();
  // Parcours noeuds, les enregistre dans le KdTree
  FldArrayF* coords = new FldArrayF(nptsTot,3);
  E_Float* cx = coords->begin(1);
  E_Float* cy = coords->begin(2);
  E_Float* cz = coords->begin(3);

  E_Int npts;
  E_Int ig = 0;
  PyObject* indirZones=NULL;

  if (extended == 0)
  {
    for (E_Int noz = 0; noz < nzones; noz++)
    {
      npts = fields[noz]->getSize(); 
      E_Float* xp = fields[noz]->begin(posxt[noz]);
      E_Float* yp = fields[noz]->begin(posyt[noz]);
      E_Float* zp = fields[noz]->begin(poszt[noz]);
      for (E_Int i = 0; i < npts; i++)
      {
        cx[ig] = xp[i]; cy[ig] = yp[i]; cz[ig] = zp[i];
        ig++;
      }
    }
  }
  else
  {
    FldArrayI* numZones = new FldArrayI(nptsTot);
    E_Int* ptrNumZones = numZones->begin();
    for (E_Int noz = 0; noz < nzones; noz++)
    {
      npts = fields[noz]->getSize(); 
      E_Float* xp = fields[noz]->begin(posxt[noz]);
      E_Float* yp = fields[noz]->begin(posyt[noz]);
      E_Float* zp = fields[noz]->begin(poszt[noz]);
      for (E_Int i = 0; i < npts; i++)
      {
        cx[ig] = xp[i]; cy[ig] = yp[i]; cz[ig] = zp[i];
        ptrNumZones[ig] = noz;
        ig++;
      }
    }
    indirZones = K_NUMPY::buildNumpyArray(*numZones);
  }
  ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*coords, 1, 2, 3); // ref sur coords
  K_SEARCH::KdTree<FldArrayF>* globalKdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc); // ref sur coordAcc

  PyObject* hook;
  E_Int* type = new E_Int [1]; type[0] = 102;
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
  for (E_Int no = 0; no < nzones; no++)
    RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);  
 
  if (extended == 0) return hook;
  else 
  {                   
    PyObject* res = PyList_New(0);
    PyList_Append(res, hook);
    PyList_Append(res, indirZones); Py_DECREF(indirZones);
    return res;
  }
}

// ============================================================================
/* Enregistre les centres des elements de a dans un KdTree 
   hook type=103 
   IN: a: array NGON */
// ============================================================================
PyObject* K_CONVERTER::registerAllElements(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Int extended;
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Ol", &arrays, &extended)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oi", &arrays, &extended)) return NULL;
#endif
  if (PyList_Check(arrays) == false)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createGlobalHook: arg must be a list of arrays.");
    return NULL;
  }
  
  // Extract infos from arrays
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = true;
  E_Boolean skipUnstructured = false; 
  E_Boolean skipDiffVars = false;
  vector<E_Int> rest;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, rest, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  if (isOk == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "createGlobalHook: arrays is not valid.");
    return NULL;
  }
  E_Int nzones = unstrF.size();
  for (E_Int v = 0; v < nzones; v++)
  {
    if (strcmp(eltTypet[v], "NGON") != 0)
    {
      for (E_Int nos = 0; nos < nzones; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      PyErr_SetString(PyExc_TypeError, 
                      "createGlobalHook: array must be a NGON.");
      return NULL; 
    }
  }
  E_Int posx1, posy1, posz1;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[noz]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[noz]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[noz]); posz1++;
    posxt.push_back(posx1); posyt.push_back(posy1); poszt.push_back(posz1); 
  }
  E_Int neltsTot = 0;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    E_Int* ptr0 = cnt[noz]->begin();
    neltsTot+= ptr0[ptr0[1]+2];
  }
  
  // Parcours les elements, calcule les centres, les enregistre dans le KdTree
  FldArrayF* centers = new FldArrayF(neltsTot,3);
  E_Float* cx = centers->begin(1);
  E_Float* cy = centers->begin(2);
  E_Float* cz = centers->begin(3);
  E_Int ig = 0;
  E_Int nv, ind, nf, pos, c;
  E_Float xf, yf, zf, inv;
  FldArrayI posFace;
  E_Int* ptrFace;
  PyObject* indirZones=NULL;

  if (extended == 0) 
  {
    for (E_Int noz = 0; noz < nzones; noz++)
    {    
      E_Int* ptr = cnt[noz]->begin();
      E_Int nelts = ptr[ptr[1]+2]; 
      E_Float* xp = unstrF[noz]->begin(posxt[noz]);
      E_Float* yp = unstrF[noz]->begin(posyt[noz]);
      E_Float* zp = unstrF[noz]->begin(poszt[noz]);
      E_Int* ptrElt = cnt[noz]->begin(); ptrElt += ptrElt[1]+4;
      K_CONNECT::getPosFaces(*cnt[noz], posFace);
    
      for (E_Int i = 0; i < nelts; i++)
      {
        nf = ptrElt[0];
        xf = 0.; yf = 0.; zf = 0.; c = 0;

        for (E_Int n = 1; n <= nf; n++)
        { 
          ind = ptrElt[n]-1;
          pos = posFace[ind];
          ptrFace = &ptr[pos];
          nv = ptrFace[0];
          for (E_Int p = 1; p <= nv; p++)
          {
            ind = ptrFace[p]-1; xf += xp[ind]; yf += yp[ind]; zf += zp[ind]; c++;
          }
        }
        inv = 1./c; xf *= inv; yf *= inv; zf *= inv;
        cx[ig] = xf; cy[ig] = yf; cz[ig] = zf;
        ptrElt += nf+1; ig++;
      }
    }
  }
  else 
  {
    FldArrayI* numZones = new FldArrayI(neltsTot);
    E_Int* ptrNumZones = numZones->begin();

    for (E_Int noz = 0; noz < nzones; noz++)
    {    
      E_Int* ptr = cnt[noz]->begin();
      E_Int nelts = ptr[ptr[1]+2]; 
      E_Float* xp = unstrF[noz]->begin(posxt[noz]);
      E_Float* yp = unstrF[noz]->begin(posyt[noz]);
      E_Float* zp = unstrF[noz]->begin(poszt[noz]);
      E_Int* ptrElt = cnt[noz]->begin(); ptrElt += ptrElt[1]+4;
      K_CONNECT::getPosFaces(*cnt[noz], posFace);
    
      for (E_Int i = 0; i < nelts; i++)
      {
        nf = ptrElt[0];
        xf = 0.; yf = 0.; zf = 0.; c = 0;

        for (E_Int n = 1; n <= nf; n++)
        { 
          ind = ptrElt[n]-1;
          pos = posFace[ind];
          ptrFace = &ptr[pos];
          nv = ptrFace[0];
          for (E_Int p = 1; p <= nv; p++)
          {
            ind = ptrFace[p]-1; xf += xp[ind]; yf += yp[ind]; zf += zp[ind]; c++;
          }
        }
        inv = 1./c; xf *= inv; yf *= inv; zf *= inv;
        cx[ig] = xf; cy[ig] = yf; cz[ig] = zf;
        ptrNumZones[ig] = noz;
        ptrElt += nf+1; ig++;
      }
    }
    indirZones = K_NUMPY::buildNumpyArray(*numZones);
  }

  ArrayAccessor<FldArrayF>* coordAcc = 
    new ArrayAccessor<FldArrayF>(*centers, 1, 2, 3); // ref sur centers
  K_SEARCH::KdTree<FldArrayF>* globalKdt = 
    new K_SEARCH::KdTree<FldArrayF>(*coordAcc); // ref sur coordAcc

  PyObject* hook;
  E_Int* type = new E_Int [1]; type[0] = 103;

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
  for (E_Int noz = 0; noz < nzones; noz++)
    RELEASESHAREDU(objut[noz], unstrF[noz], cnt[noz]);

  if (extended == 0) return hook;
  else
  {                   
    PyObject* res = PyList_New(0);
    PyList_Append(res,hook);
    PyList_Append(res, indirZones); Py_DECREF(indirZones);
    return res;
  }
}

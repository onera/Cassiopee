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

# include "kcore.h"
# include "Nuga/include/KdTree.h"
#include "parallel.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Correspondance d'indice entre les vertex d'un array structure et les 
   vertex d'un array non structure equivalent.
   IN: structArray: array structure dont on cherche la correspondance
   IN: unstrArray: array non structure equivalent
   IN: arrayOfStructIndices: numpy array des indices des vertex dans le 
   maillage structure dont on cherche la correspondance dans le maillage non
   structure
   OUT: numpy array des indices des vertex du maillage non structure 
   correspondants
 */
//=============================================================================
PyObject* K_KCORE::indiceStruct2Unstr(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject *arrayOfStructIndices, *structArray, *unstrArray;
  E_Float eps;
  if (!PyArg_ParseTuple(args, "OOOd", &structArray, &unstrArray, 
                        &arrayOfStructIndices, &eps)) return NULL;

  // Maillage structure
  FldArrayF* f1; FldArrayI* cn1;
  E_Int ni1, nj1, nk1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray(structArray, varString1, 
                                     f1, ni1, nj1, nk1, cn1, eltType1, true);
  if (res1 != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceStruct2Unstr: 1st arg must be structured.");
    if (res1 == 2) RELEASESHAREDU(structArray, f1, cn1);
    return NULL;
  }
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceStruct2Unstr: 1st arg must contain coordinates.");
    RELEASESHAREDS(structArray, f1); return NULL;
  }
  posx1++; posy1++; posz1++;

  // Maillage non structure 
  FldArrayF* f2; FldArrayI* cn2;
  E_Int ni2, nj2, nk2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray(unstrArray, varString2, 
                                     f2, ni2, nj2, nk2, cn2, eltType2, true);
  if (res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceStruct2Unstr: 2nd arg must be unstructured.");
    if (res2 == 1) RELEASESHAREDS(unstrArray, f2);
    RELEASESHAREDS(structArray, f1); return NULL;
  }
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceStruct2Unstr: 2nd arg must contain coordinates.");
    RELEASESHAREDS(structArray, f1); RELEASESHAREDU(unstrArray, f2, cn2); return NULL;
  }
  posx2++; posy2++; posz2++;

  // Indices des points du maillage structure a traiter
  E_Int nind, nf; E_Int* indices;
  E_Int ret = K_NUMPY::getFromNumpyArray(arrayOfStructIndices, indices, 
                                         nind, nf, true);

  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceStruct2Unstr: invalid indices.");
    RELEASESHAREDS(structArray, f1);
    RELEASESHAREDU(unstrArray, f2, cn2);
    return NULL;
  }

  ArrayAccessor<FldArrayF> coordAcc(*f2, posx2, posy2, posz2);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);
  E_Float* xt1 = f1->begin(posx1);
  E_Float* yt1 = f1->begin(posy1); 
  E_Float* zt1 = f1->begin(posz1); 
  E_Float* xt2 = f2->begin(posx2); 
  E_Float* yt2 = f2->begin(posy2);
  E_Float* zt2 = f2->begin(posz2); 

  // Sortie dans un numpy
  PyObject* a = K_NUMPY::buildNumpyArray(nind, 1, 1);
  E_Int* ptr = K_NUMPY::getNumpyPtrI(a);

#pragma omp parallel default(shared)
  {
    E_Float pt[3]; 
    E_Int ind, indu;
    E_Float dx, dy, dz;
#pragma omp for
    for (E_Int noind = 0; noind < nind; noind++)
    {
      ind = indices[noind];
      pt[0] = xt1[ind]; pt[1] = yt1[ind]; pt[2] = zt1[ind];
      indu = globalKdt.getClosest(pt);
      dx = xt2[indu]-pt[0];
      dy = yt2[indu]-pt[1];
      dz = zt2[indu]-pt[2];
      if (dx*dx + dy*dy + dz*dz > eps) indu = -1;
      ptr[noind] = indu;
    }
  }
  RELEASESHAREDS(structArray, f1); RELEASESHAREDU(unstrArray, f2, cn2);
  Py_DECREF(arrayOfStructIndices);
  return a;
}

//=============================================================================
/* Correspondance d'indice entre les vertex d'arrays structures et les 
   vertex d'un array non structure equivalent.
   IN: structArrays: arrays structures dont on cherche la correpondance
   IN: unstrArray: array non structure equivalent
   OUT: numpy array des indices des vertex du maillage non structure 
   correspondants
 */
//=============================================================================
PyObject* K_KCORE::indiceStruct2Unstr2(PyObject* self, PyObject* args)
{
  PyObject *structArrays, *unstrArray;
  E_Float eps;
  if (!PyArg_ParseTuple(args, "OOd", &structArrays, &unstrArray, &eps)) 
    return NULL;

  // Maillages structures
  vector<E_Int> res;
  vector<char*> structVarString, unstructVarString;
  vector<FldArrayF*> structF, unstructF;
  vector<E_Int> ni, nj, nk;
  vector<FldArrayI*> cn; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  K_ARRAY::getFromArrays(
    structArrays, res, structVarString, unstructVarString,
    structF, unstructF, ni, nj, nk, cn, eltType, objs, obju,
    false, true, false, true, true);

  // Maillage non structure 
  FldArrayF* f2; FldArrayI* cn2;
  E_Int ni2, nj2, nk2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray(unstrArray, varString2, 
                                     f2, ni2, nj2, nk2, cn2, eltType2, true);
  if (res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceStruct2Unstr2: 2nd arg must be unstructured.");
    if (res2 == 1) RELEASESHAREDS(unstrArray, f2);
    for (unsigned int z = 0; z < structF.size(); z++)
      RELEASESHAREDS(objs[z], structF[z]); 
    return NULL;
  }
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceStruct2Unstr2: 2nd arg must contain coordinates.");
    RELEASESHAREDU(unstrArray, f2, cn2);
    for (unsigned int z = 0; z < structF.size(); z++)
      RELEASESHAREDS(objs[z], structF[z]); 
    return NULL;
  }
  posx2++; posy2++; posz2++;

  ArrayAccessor<FldArrayF> coordAcc(*f2, posx2, posy2, posz2);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);
  E_Float* xt2 = f2->begin(posx2);
  E_Float* yt2 = f2->begin(posy2);
  E_Float* zt2 = f2->begin(posz2); 

  PyObject* out = PyList_New(0);
  //npy_intp dim[1];

  for (unsigned int z = 0; z < structF.size(); z++)
  {
    E_Int posx1 = K_ARRAY::isCoordinateXPresent(structVarString[z]);
    E_Int posy1 = K_ARRAY::isCoordinateYPresent(structVarString[z]);
    E_Int posz1 = K_ARRAY::isCoordinateZPresent(structVarString[z]);
    posx1++; posy1++; posz1++;
    E_Float* xt1 = structF[z]->begin(posx1); 
    E_Float* yt1 = structF[z]->begin(posy1); 
    E_Float* zt1 = structF[z]->begin(posz1);

    // Sortie dans un numpy
    E_Int nind = ni[z]*nj[z]*nk[z];
    PyObject* a = K_NUMPY::buildNumpyArray(nind, 1, 1);
    E_Int* ptr = K_NUMPY::getNumpyPtrI(a);

#pragma omp parallel default(shared)
    {
      E_Float pt[3];
      E_Int indu;
      E_Float dx, dy, dz;
#pragma omp for schedule(dynamic)
      for (E_Int ind = 0; ind < nind; ind++)
      {
        pt[0] = xt1[ind]; pt[1] = yt1[ind]; pt[2] = zt1[ind];
        indu = globalKdt.getClosest(pt);
        dx = xt2[indu]-pt[0];
        dy = yt2[indu]-pt[1];
        dz = zt2[indu]-pt[2];
        if (dx*dx + dy*dy + dz*dz > eps) indu = -1;
        ptr[ind] = indu;
      }
    }

    PyList_Append(out, a);
    Py_DECREF(a);
  }

  RELEASESHAREDU(unstrArray, f2, cn2);
  for (unsigned int z = 0; z < structF.size(); z++)
    RELEASESHAREDS(objs[z], structF[z]);

  return out;
}

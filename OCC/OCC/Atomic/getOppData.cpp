/*    
    Copyright 2013-2024 Onera.

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

#include "occ.h"
#include <set>

// ============================================================================
/* Return the opposite data for an edge on a face */
// ============================================================================
PyObject* K_OCC::getOppData(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* indArray;
  if (!PYPARSETUPLE_(args, OO_ , &array, &indArray)) return NULL;

  // input array is the face TRI mesh (without ghost cells)
  // inds is the indices in TRI mesh of edge

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, 
                                     cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: invalid array.");
    return NULL;
  }

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: only for TRI array.");
    return NULL;
  }

  if (strcmp(eltType, "TRI") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: only for TRI array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: coordinates missing in TRI array.");
    return NULL;
  }
  E_Float* fx = f->begin(posx+1);
  E_Float* fy = f->begin(posy+1);
  E_Float* fz = f->begin(posz+1);

  // Check inds
  E_Int* inds; E_Int size; E_Int nfld;
  E_Int res2 = K_NUMPY::getFromNumpyArray(indArray, inds, size, nfld, true);
  if (res2 == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: invalid index array.");
    return NULL;
  }

  // Export numpy of coordinates of inds (first)
  PyObject* tplx = K_NUMPY::buildNumpyArray(size, 1, 0, 1);
  E_Float* px = K_NUMPY::getNumpyPtrF(tplx);
  PyObject* tply = K_NUMPY::buildNumpyArray(size, 1, 0, 1);
  E_Float* py = K_NUMPY::getNumpyPtrF(tply);
  PyObject* tplz = K_NUMPY::buildNumpyArray(size, 1, 0, 1);
  E_Float* pz = K_NUMPY::getNumpyPtrF(tplz);
  
  // remplissage
  E_Int ind;
  for (E_Int i = 0; i < size; i++)
  {
    ind = inds[i];
    px[i] = fx[ind];
    py[i] = fy[ind];
    pz[i] = fz[ind];
  }

  // Recherche des vertex voisins
  E_Int nv = f->getSize();
  printf("nv=%d \n", nv);
  
  std::vector< std::vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(*cn, cVE);

  std::set<E_Int> myset;
  E_Int indn, inde;
  FldArrayI& cnp = *cn;
  for (E_Int i = 0; i < size; i++)
  {
    ind = inds[i];
    //printf("ind=%d size=%d\n", ind, cVE.size());
    std::vector<E_Int>& elts = cVE[ind];
    for (size_t e = 0; e < elts.size(); e++)
    {
      inde = elts[e];
      //printf("inde=%d size=\n", inde);
      indn = cnp(inde, 1)-1;
      myset.insert(indn);
      indn = cnp(inde, 2)-1;
      myset.insert(indn);
      indn = cnp(inde, 3)-1;
      myset.insert(indn);
    }
  }
  
  // Export set to numpy
  E_Int vsize = myset.size();

  // Export numpy of coordinates of inds (first)
  PyObject* tpl2x = K_NUMPY::buildNumpyArray(vsize, 1, 0, 1);
  E_Float* p2x = K_NUMPY::getNumpyPtrF(tpl2x);
  PyObject* tpl2y = K_NUMPY::buildNumpyArray(vsize, 1, 0, 1);
  E_Float* p2y = K_NUMPY::getNumpyPtrF(tpl2y);
  PyObject* tpl2z = K_NUMPY::buildNumpyArray(vsize, 1, 0, 1);
  E_Float* p2z = K_NUMPY::getNumpyPtrF(tpl2z);

  E_Int i = 0;
  for (auto it = myset.begin(); it != myset.end(); it++)
  {
    ind = *it;
    p2x[i] = fx[ind];
    p2y[i] = fy[ind];
    p2z[i] = fz[ind]; i++;
  }

  // Recherche des triangles supplementaires
  

  RELEASESHAREDU(array, f, cn);
  Py_DECREF(indArray);

  return Py_BuildValue("OOO", tpl2x, tpl2y, tpl2z);
} 
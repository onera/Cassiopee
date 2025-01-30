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

#include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Tag defined BC for a zone: definedBC variable is set to 1 */
//=============================================================================
PyObject* K_CONVERTER::tagDefinedBC(PyObject* self, PyObject* args)
{
  PyObject *allwins, *array;
  E_Int dimPb;
  if (!PYPARSETUPLE_(args, OO_ I_, &array, &allwins, &dimPb)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, 
                                    cn, eltType, true);
  
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "tagDefinedBC: 1st arg must be structured.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  E_Int posd = K_ARRAY::isNamePresent("definedBC", varString);
  if (posd == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "tagDefinedBC: definedBC variable does not exist.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posd++;
  E_Int dim = 3; if (km == 1 || dimPb == 2) dim = 2;
  
  // Ranges des CL deja definies
  if (PyList_Check(allwins) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "tagDefinedBC: allwins must be a list.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  E_Int nwins = PyList_Size(allwins);
  FldArrayI ranges(nwins, 6); ranges.setAllValuesAt(1);

  for (int i = 0; i < nwins; i++)
  {
    E_Int i0 = E_Int(i);
    // recuperation de chq win [i1,i2,j1,j2,k1,k2]
    PyObject* win = PyList_GetItem(allwins, i);
    if (PyList_Check(win) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "tagDefinedBC: each win must be a list.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    E_Int sizewin = PyList_Size(win);
    if (sizewin == 6  || sizewin == 4)
    {
      for (int j = 0; j < sizewin; j++)
      {    
        E_Int j0 = E_Int(j)+1;
        PyObject* iwin = PyList_GetItem(win, j);
        ranges(i0, j0) = PyInt_AsLong(iwin); 
      }
    }// dimensions of win OK
    else 
    {
      PyErr_SetString(PyExc_TypeError, 
                      "tagDefinedBC: wrong dimensions for a window.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
  }// for all wins in allwins

  E_Int nfld = f->getNfld();
  PyObject* tpl = K_ARRAY::buildArray(nfld, varString, im, jm, km);
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(f->getSize(), nfld, fnp, true); fn = *f;
  
  if (dim == 3) tagDefinedBC3D(posd, im, jm, km, fn, ranges);
  else tagDefinedBC2D(posd, im, jm, fn, ranges);
  
  RELEASESHAREDB(res, array, f, cn); 
  return tpl;
}
//=============================================================================
/* Tag points to 2 between indices in wins - 2D case */
//=============================================================================
void K_CONVERTER::tagDefinedBC2D(E_Int posd, E_Int ni, E_Int nj, 
                                 FldArrayF& f, FldArrayI& wins)
{
  E_Float* tag = f.begin(posd);
  E_Int* it1 = wins.begin(1);
  E_Int* it2 = wins.begin(2);
  E_Int* jt1 = wins.begin(3);
  E_Int* jt2 = wins.begin(4);
  E_Int nwins = wins.getSize();
  E_Int i1, i2, j1, j2, is1, js1;
  E_Int ind;
  for (E_Int now = 0; now < nwins; now++)
  {
    i1 = it1[now]-1; i2 = it2[now]-1; j1 = jt1[now]-1; j2 = jt2[now]-1; 
    is1 = K_FUNC::E_min(i1+1,i2); js1 = K_FUNC::E_min(j1+1,j2);
    if (i1 == i2) 
    {
      for (E_Int j = js1; j < j2; j++)
      {ind = i1 + j*ni; tag[ind] = 2.;}
      
      ind = i1 + j1*ni; tag[ind] = 1.;
      ind = i1 + j2*ni; tag[ind] = 1.;
    }
    else if (j1 == j2)
    {
      for (E_Int i = is1; i < i2; i++)
      {ind = i + j1*ni; tag[ind] = 2.;}
      
      ind = i1 + j1*ni; tag[ind] = 1.;
      ind = i2 + j1*ni; tag[ind] = 1.;
    }
  }
  return;
}

//=============================================================================
/* Tag points to 1 between indices in wins - 3D case */
//=============================================================================
void K_CONVERTER::tagDefinedBC3D(E_Int posd, E_Int ni, E_Int nj, E_Int nk, 
                                 FldArrayF& f, FldArrayI& wins)
{
//#define INCTAG(a) ((a > 0) ? 2 : 1)
#define INCTAG(a) 1

  E_Float* tag = f.begin(posd);
  E_Int ninj = ni*nj;
  E_Int* it1 = wins.begin(1);
  E_Int* it2 = wins.begin(2);
  E_Int* jt1 = wins.begin(3);
  E_Int* jt2 = wins.begin(4);
  E_Int* kt1 = wins.begin(5);
  E_Int* kt2 = wins.begin(6);
  E_Int nwins = wins.getSize();
  E_Int i1, i2, j1, j2, k1, k2, is1, js1,  ks1;
  E_Int ind;
  for (E_Int now = 0; now < nwins; now++)
  {
    i1 = it1[now]-1; j1 = jt1[now]-1; k1 = kt1[now]-1;
    i2 = it2[now]-1; j2 = jt2[now]-1; k2 = kt2[now]-1;
    
    ks1 = K_FUNC::E_min(k1+1,k2); js1 = K_FUNC::E_min(j1+1,j2); is1 = K_FUNC::E_min(i1+1,i2);
    if (i1 == i2) 
    {
      for (E_Int k = ks1; k < k2; k++)
        for (E_Int j = js1; j < j2; j++)
        {
          ind = i1 + j*ni + k*ninj;
          tag[ind] = 2.;
        }
      
      for (E_Int j = j1; j < j2+1; j++)
      {
        ind = i1 + j*ni + k1*ninj; tag[ind] = INCTAG(tag[ind]); 
        ind = i1 + j*ni + k2*ninj; tag[ind] = INCTAG(tag[ind]);
      }
      for (E_Int k = k1; k < k2+1; k++)
      {
        ind = i1 + j1*ni + k*ninj; tag[ind] = INCTAG(tag[ind]);
        ind = i1 + j2*ni + k*ninj; tag[ind] = INCTAG(tag[ind]);
      }
    }
    else if (j1 == j2) 
    {
      for (E_Int k = ks1; k < k2; k++)
        for (E_Int i = is1; i < i2; i++)
        {
          ind = i + j1*ni + k*ninj; tag[ind] = 2.;
        }
      for (E_Int i = i1; i < i2+1; i++)
      {
        ind = i + j1*ni + k1*ninj; tag[ind] = INCTAG(tag[ind]);
        ind = i + j1*ni + k2*ninj; tag[ind] = INCTAG(tag[ind]);
      }
      for (E_Int k = k1; k < k2+1; k++)
      {
        ind = i1 + j1*ni + k*ninj; tag[ind] = INCTAG(tag[ind]);
        ind = i2 + j1*ni + k*ninj; tag[ind] = INCTAG(tag[ind]);
      } 
    }
    else 
    {
      for (E_Int j = js1; j < j2; j++)
        for (E_Int i = is1; i < i2; i++)
        {
          ind = i + j*ni + k1*ninj; tag[ind] = 2.;    
        }
      for (E_Int i = i1; i < i2+1; i++)
      {
        ind = i + j1*ni + k1*ninj; tag[ind] = INCTAG(tag[ind]);
        ind = i + j2*ni + k1*ninj; tag[ind] = INCTAG(tag[ind]);
      }
      for (E_Int j = j1; j < j2+1; j++)
      {
        ind = i1 + j*ni + k1*ninj; tag[ind] = INCTAG(tag[ind]);
        ind = i2 + j*ni + k1*ninj; tag[ind] = INCTAG(tag[ind]);
      }
    }// k1 = k2
  }//now
}

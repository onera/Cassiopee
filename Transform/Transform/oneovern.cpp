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

# include "transform.h"

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Prend 1 point sur N 
   Si le nbre de points ne tombe pas juste, le dernier plan est ajoute.
*/
// ============================================================================
PyObject* K_TRANSFORM::oneovern(PyObject* self, PyObject* args)
{
  E_Int Ni, Nj, Nk, Add;
  PyObject* array;
  if (!PYPARSETUPLEI(args,
                    "O(lll)l", "O(iii)i",
                    &array, &Ni, &Nj, &Nk, &Add))
  {
      return NULL;
  }
  // Get E_Int value
  E_Int ni = Ni; E_Int nj = Nj; E_Int nk = Nk;

  if (Ni <= 0 || Nj <= 0 || Nk <= 0)
  {
    PyErr_SetString(PyExc_ValueError,
                    "oneovern: Ni, Nj, Nk must be >= 1.");
    return NULL;
  }
  E_Int add = Add;
  if ( add != 0 && add != 1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "oneovern: invalid value for add.");
    return NULL;
  }
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray2(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "oneovern: unknown type of array.");
    return NULL;
  }

  if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "oneovern: can not be used on an unstructured array.");
    return NULL;
  }
  
  E_Int addi=0, addj=0, addk=0;
  E_Int imjm = im*jm;
  E_Int in = E_Int(floor(1.*(im-1)/ni)+1);
  E_Int jn = E_Int(floor(1.*(jm-1)/nj)+1);
  E_Int kn = E_Int(floor(1.*(km-1)/nk)+1);

  if (add == 1)
  {
    if (im - in*ni != 1-ni) addi = 1; 
    if (jm - jn*nj != 1-nj) addj = 1; 
    if (km - kn*nk != 1-nk) addk = 1;
  }
  E_Int nfld = f->getNfld();
  //printf("size %d %d %d\n", in+addi, jn+addj, kn+addk);

  PyObject* tpl;
  tpl = K_ARRAY::buildArray2(nfld, varString, in+addi, jn+addj, kn+addk);
  E_Float* subzonep = K_ARRAY::getFieldPtr(tpl);
  FldArrayF subzone((in+addi)*(jn+addj)*(kn+addk), nfld, subzonep, true);
  E_Int ind, ind2;

  for (E_Int n = 1; n <= nfld; n++)
  {
    E_Float* sp = subzone.begin(n);
    E_Float* fp = f->begin(n);
    ind2 = 0;

    for (E_Int k = 0; k < km; k = k+nk)
    {
      for (E_Int j = 0; j < jm; j = j+nj)
      {
        for (E_Int i = 0; i < im; i = i+ni)
        {
          ind = i + j*im + k*imjm;
          sp[ind2] = fp[ind]; ind2++;
        }
        if (addi == 1) 
        {ind = im-1 + j*im + k*imjm; sp[ind2] = fp[ind]; ind2++;}
      }
      if (addj == 1)
      {
        for (E_Int i = 0; i < im; i = i+ni)
        {
          ind = i + (jm-1)*im + k*imjm;
          sp[ind2] = fp[ind]; ind2++;
        }
        if (addi == 1) 
        {ind = im-1 + (jm-1)*im + k*imjm; sp[ind2] = fp[ind]; ind2++;}
      }
    }
    if (addk == 1)
    {
      for (E_Int j = 0; j < jm; j = j+nj)
      {
        for (E_Int i = 0; i < im; i = i+ni)
        {
          ind = i + j*im + (km-1)*imjm;
          sp[ind2] = fp[ind]; ind2++;
        }
        if (addi == 1) 
        {ind = im-1 + j*im + (km-1)*imjm; sp[ind2] = fp[ind]; ind2++;}
      }
      if (addj == 1)
      {
        for (E_Int i = 0; i < im; i = i+ni)
        {
          ind = i + (jm-1)*im + (km-1)*imjm;
          sp[ind2] = fp[ind]; ind2++;
        }
        if (addi == 1) 
        {ind = im-1 + (jm-1)*im + (km-1)*imjm; sp[ind2] = fp[ind]; ind2++;}
      }
    }
  }
  
  RELEASESHAREDS(array, f);
  return tpl;
}

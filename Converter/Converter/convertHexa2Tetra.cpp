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

// Convert unstructured hexa array to tetra array

# include <stdlib.h>
# include <string.h>
# include <vector>

# include "converter.h"
# include "kcore.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert unstructured hexa array to a tetraedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertHexa2Tetra(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  E_Int res;

  E_Int elt = 1;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cnl, eltType, true);

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertHexa2Tetra: array must be unstructured Hexa or Quad.");
    return NULL; 
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertHexa2Tetra: array is invalid.");
    return NULL;
  }

  E_Int dim0 = -1;
  if (strcmp(eltType, "HEXA") != 0 && strcmp(eltType, "QUAD") != 0)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_ValueError,
                    "convertHexa2Tetra: array must be hexa.");
    return NULL;
  }
      
  // 2D or 3D ?
  dim0 = 3;
  if (strcmp(eltType, "QUAD") == 0) dim0 = 2;
      
  // build the unstructured mesh
  E_Int ncells = cnl->getSize(); // nb de cellules structurees
  FldArrayI& cl = *cnl;

  // build array (shared)
  E_Int nelts, nvert;
  if (dim0 == 2) {nelts = 2*ncells; elt = 2; nvert = 3;}// TRI
  else {nelts = 5*ncells; elt = 4; nvert = 4;}//TETRA
  PyObject* tpl = K_ARRAY::buildArray(f->getNfld(), varString, f->getSize(), nelts, elt, NULL);
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(f->getSize(), f->getNfld(), fnp, true); fn = *f;
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI cn(nelts, nvert, cnnp, true);

  if (dim0 == 2)
  {
    // define the type TRI: tetra 2d 
    E_Int* cl1 = cl.begin(1);
    E_Int* cl2 = cl.begin(2);
    E_Int* cl3 = cl.begin(3);
    E_Int* cl4 = cl.begin(4);
    E_Int* cn1 = cn.begin(1);
    E_Int* cn2 = cn.begin(2);
    E_Int* cn3 = cn.begin(3);

#pragma omp parallel default(shared) if (ncells > __MIN_SIZE_MEAN__)
    {
      E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        //starts from 1
        ind1 = cl1[i]; 
        ind2 = cl2[i];
        ind3 = cl3[i];
        ind4 = cl4[i];
        
        ind = 2*i;
        cn1[ind] = ind1;
        cn2[ind] = ind2;
        cn3[ind] = ind3;
        
        ind++;
        cn1[ind] = ind1;
        cn2[ind] = ind3;
        cn3[ind] = ind4;
      }
    }
  }//dim 2
  else 
  { 
    E_Int* cl1 = cl.begin(1);
    E_Int* cl2 = cl.begin(2);
    E_Int* cl3 = cl.begin(3);
    E_Int* cl4 = cl.begin(4);
    E_Int* cl5 = cl.begin(5);
    E_Int* cl6 = cl.begin(6);
    E_Int* cl7 = cl.begin(7);
    E_Int* cl8 = cl.begin(8);

    E_Int* cn1 = cn.begin(1);
    E_Int* cn2 = cn.begin(2);
    E_Int* cn3 = cn.begin(3);
    E_Int* cn4 = cn.begin(4);
#pragma omp parallel default(shared) if (ncells > __MIN_SIZE_MEAN__)
    {
      E_Int ind, sum, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
#pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        sum = i;
        ind1 = cl1[i];
        ind2 = cl2[i];
        ind3 = cl3[i];
        ind4 = cl4[i];
        ind5 = cl5[i];
        ind6 = cl6[i];
        ind7 = cl7[i];
        ind8 = cl8[i];
        
        if (sum%2 == 0) // pair 
        {
          //tetra ABDE
          ind = 5*i;
          cn1[ind] = ind1;
          cn2[ind] = ind2;
          cn3[ind] = ind4;
          cn4[ind] = ind5;
          
          //tetra BCDG
          ind++;
          cn1[ind] = ind2;
          cn2[ind] = ind3;
          cn3[ind] = ind4;
          cn4[ind] = ind7;
          
          //tetra DEGH
          ind++;
          cn1[ind] = ind4;
          cn2[ind] = ind5;
          cn3[ind] = ind7;
          cn4[ind] = ind8;
          
          //tetra BEFG
          ind++;
          cn1[ind] = ind2;
          cn2[ind] = ind5;
          cn3[ind] = ind6;
          cn4[ind] = ind7;
          
          //tetra BDEG
          ind++;
          cn1[ind] = ind2;
          cn2[ind] = ind4;
          cn3[ind] = ind5;
          cn4[ind] = ind7;
        }
        else // impair 
        {
          //tetra ACDH : 1348
          ind = 5*i;
          cn1[ind] = ind1;
          cn2[ind] = ind3;
          cn3[ind] = ind4;
          cn4[ind] = ind8;
          
          //tetra AFBC : 1623
          ind++;
          cn1[ind] = ind1;
          cn2[ind] = ind2;
          cn3[ind] = ind3;
          cn4[ind] = ind6;
          //cn1[ind] = ind1;
          //cn2[ind] = ind6;
          //cn3[ind] = ind2;
          //cn4[ind] = ind3;

          //tetra HFGC : 8673
          ind++;
          cn1[ind] = ind3;
          cn2[ind] = ind6;
          cn3[ind] = ind7;
          cn4[ind] = ind8;
          //cn1[ind] = ind8;
          //cn2[ind] = ind6;
          //cn3[ind] = ind7;
          //cn4[ind] = ind3;

          //tetra FHAE : 6815
          ind++;
          cn1[ind] = ind6;
          cn2[ind] = ind5;
          cn3[ind] = ind8;
          cn4[ind] = ind1;
          //cn1[ind] = ind6;
          //cn2[ind] = ind8;
          //cn3[ind] = ind1;
          //cn4[ind] = ind5;

          //tetra FHAC : 6813
          ind++;
          cn1[ind] = ind6;
          cn2[ind] = ind3;
          cn3[ind] = ind1;
          cn4[ind] = ind8;
          //cn1[ind] = ind6;
          //cn2[ind] = ind8;
          //cn3[ind] = ind1;
          //cn4[ind] = ind3;
        }
      }
    }
  } //dim = 3
  
  RELEASESHAREDU(array, f, cnl);
  return tpl;
}

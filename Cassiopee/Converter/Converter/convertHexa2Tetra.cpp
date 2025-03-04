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
  E_Int nil, njl, nkl, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType; string eltType2;

  res = K_ARRAY::getFromArray3(array, varString, f,
                               nil, njl, nkl, cnl, eltType);

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
  dim0 = 3; eltType2 = "TETRA";
  if (strcmp(eltType, "QUAD") == 0) 
  {
    dim0 = 2; eltType2 = "TRI";
  }
      
  // Build the unstructured mesh
  FldArrayI& cm = *(cnl->getConnect(0));
  E_Int ncells = cm.getSize();

  // build array (shared)
  E_Int nelts;
  if (dim0 == 2) nelts = 2*ncells; // TRI
  else {nelts = 5*ncells;} // TETRA
  E_Int npts = f->getSize();
  E_Int api = f->getApi();
  E_Int nfld = f->getNfld();

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts,
                                       eltType2.c_str(), false, api);
  FldArrayF* f2; FldArrayI* cnl2;
  K_ARRAY::getFromArray3(tpl, f2, cnl2);
  FldArrayI& cm2 = *(cnl2->getConnect(0));

#pragma omp parallel default(shared) if (ncells > __MIN_SIZE_MEAN__)
  {
    if (dim0 == 2)
    {
      // define the type TRI: tetra 2d 
      E_Int ind, ind1, ind2, ind3, ind4;
#pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        //starts from 1
        ind1 = cm(i,1); 
        ind2 = cm(i,2);
        ind3 = cm(i,3);
        ind4 = cm(i,4);
        
        ind = 2*i;
        cm2(ind,1) = ind1;
        cm2(ind,2) = ind2;
        cm2(ind,3) = ind3;
        
        ind++;
        cm2(ind,1) = ind1;
        cm2(ind,2) = ind3;
        cm2(ind,3) = ind4;
      }
    }//dim 2
    else 
    { 
      E_Int ind, sum, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
#pragma omp for
      for (E_Int i = 0; i < ncells; i++)
      {
        sum = i;
        ind1 = cm(i,1);
        ind2 = cm(i,2);
        ind3 = cm(i,3);
        ind4 = cm(i,4);
        ind5 = cm(i,5);
        ind6 = cm(i,6);
        ind7 = cm(i,7);
        ind8 = cm(i,8);
        
        if (sum%2 == 0) // pair 
        {
          //tetra ABDE
          ind = 5*i;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind4;
          cm2(ind,4) = ind5;
          
          //tetra BCDG
          ind++;
          cm2(ind,1) = ind2;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
          cm2(ind,4) = ind7;
          
          //tetra DEGH
          ind++;
          cm2(ind,1) = ind4;
          cm2(ind,2) = ind5;
          cm2(ind,3) = ind7;
          cm2(ind,4) = ind8;
          
          //tetra BEFG
          ind++;
          cm2(ind,1) = ind2;
          cm2(ind,2) = ind5;
          cm2(ind,3) = ind6;
          cm2(ind,4) = ind7;
          
          //tetra BDEG
          ind++;
          cm2(ind,1) = ind2;
          cm2(ind,2) = ind4;
          cm2(ind,3) = ind5;
          cm2(ind,4) = ind7;
        }
        else // impair 
        {
          //tetra ACDH : 1348
          ind = 5*i;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
          cm2(ind,4) = ind8;
          
          //tetra AFBC : 1623
          ind++;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind3;
          cm2(ind,4) = ind6;
          //cm2(ind,1) = ind1;
          //cm2(ind,2) = ind6;
          //cm2(ind,3) = ind2;
          //cm2(ind,4) = ind3;

          //tetra HFGC : 8673
          ind++;
          cm2(ind,1) = ind3;
          cm2(ind,2) = ind6;
          cm2(ind,3) = ind7;
          cm2(ind,4) = ind8;
          //cm2(ind,1) = ind8;
          //cm2(ind,2) = ind6;
          //cm2(ind,3) = ind7;
          //cm2(ind,4) = ind3;

          //tetra FHAE : 6815
          ind++;
          cm2(ind,1) = ind6;
          cm2(ind,2) = ind5;
          cm2(ind,3) = ind8;
          cm2(ind,4) = ind1;
          //cm2(ind,1) = ind6;
          //cm2(ind,2) = ind8;
          //cm2(ind,3) = ind1;
          //cm2(ind,4) = ind5;

          //tetra FHAC : 6813
          ind++;
          cm2(ind,1) = ind6;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind1;
          cm2(ind,4) = ind8;
          //cm2(ind,1) = ind6;
          //cm2(ind,2) = ind8;
          //cm2(ind,3) = ind1;
          //cm2(ind,4) = ind3;
        }
      }//dim = 3
    }

    // Copy fields to f2
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
#pragma omp for
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }
  
  RELEASESHAREDU(array, f, cnl);
  RELEASESHAREDU(tpl, f2, cnl2);
  return tpl;
}
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
// Convert structured array to tetra array

# include <stdlib.h>
# include "converter.h"
# include "kcore.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert  structured array to a tetraedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertStruct2Tetra(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  PyObject* tpl = NULL; // array de sortie en noeuds
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType; string eltType2;
  E_Int res, dim0 = -1;
  res = K_ARRAY::getFromArray3(array, varString, f,
                               nil, njl, nkl, cnl, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2Tetra: array is invalid.");
    return NULL;
  }
  if (res == 2)
  {
    RELEASESHAREDU(array, f, cnl);
    PyErr_SetString(PyExc_TypeError, 
                    "convertStruct2Tetra: array must be structured.");
    return NULL; 
  }
      
  // 1D, 2D or 3D ?
  dim0 = 3; eltType2 = "TETRA";
  if (nil == 1)
  {
    if (njl == 1 || nkl == 1 ) dim0 = 1;
    else dim0 = 2;
  }
  else if (njl == 1)
  {
    if (nil == 1 || nkl == 1 ) dim0 = 1;
    else dim0 = 2;
  }
  else if (nkl == 1)
  {
    if (nil == 1 || njl == 1 ) dim0 = 1;
    else dim0 = 2;
  }

  // Build the unstructured mesh
  E_Int ni1 = E_max(1, E_Int(nil)-1);
  E_Int nj1 = E_max(1, E_Int(njl)-1);
  E_Int nk1 = E_max(1, E_Int(nkl)-1);
  E_Int ninj = nil*njl;
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
  E_Int npts = f->getSize(), api = f->getApi(), nfld = f->getNfld();

  E_Int nelts = 5*ncells; // nb d'elements non structures
  if (dim0 == 1) {eltType2 = "BAR"; nelts = ncells;}
  else if (dim0 == 2) {eltType2 = "TRI"; nelts = 2*ncells;}

  tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts,
                             eltType2.c_str(), false, api);
  FldArrayF* f2; FldArrayI* cnl2;
  K_ARRAY::getFromArray3(tpl, f2, cnl2);
  FldArrayI& cm = *(cnl2->getConnect(0));
                
#pragma omp parallel
  {
    if (dim0 == 1)
    {
      E_Int ind1, ind2;
      if (nk1 == 1 && nj1 == 1)
      {
#pragma omp for
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = i + 1;
          ind2 = ind1 + 1;
          cm(i,1) = ind1;
          cm(i,2) = ind2;
        }
      }
      else if (ni1 == 1 && nj1 == 1)
      {
#pragma omp for
        for (E_Int k = 0; k < nk1; k++)
        {
          ind1 = k*nil*njl + 1;
          ind2 = ind1 + nil*njl;
          cm(k,1) = ind1;
          cm(k,2) = ind2;
        }
      }
      else if (ni1 == 1 && nk1 == 1)
      {
#pragma omp for
        for (E_Int j = 0; j < nj1; j++)
        {
          ind1 = j*nil + 1;
          ind2 = ind1 + nil;
          cm(j,1) = ind1;
          cm(j,2) = ind2;
        }
      }
    }
    else if (dim0 == 2)
    {
      E_Int ind, ind1, ind2, ind3, ind4;
      if (nk1 == 1)
      {
        for (E_Int j = 0; j < nj1; j++)
#pragma omp for
          for (E_Int i = 0; i < ni1; i++)
          {
            //starts from 1
            ind1 = i + j*nil + 1; //(i,j,1)
            ind2 = ind1 + 1;      //(i+1,j,1)
            ind3 = ind2 + nil;    //(i+1,j+1,1)
            ind4 = ind3 - 1;      //(i,j+1,1)
            
            ind = 2*(i + j*ni1);
            cm(ind,1) = ind1;
            cm(ind,2) = ind2;
            cm(ind,3) = ind3;
            
            ind++;
            cm(ind,1) = ind1;
            cm(ind,2) = ind3;
            cm(ind,3) = ind4;
          }
      }
      else if (nj1 == 1)
      {
        for (E_Int k = 0; k < nk1; k++)
#pragma omp for
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + k*ninj + 1;  //(i,1,k)
            ind2 = ind1 + ninj;     //(i,1,k+1)
            ind3 = ind2 + 1;        //(i+1,1,k+1)
            ind4 = ind1 + 1;        //(i+1,1,k)
            
            ind = 2*(i + k*ni1);
            cm(ind,1) = ind1;
            cm(ind,2) = ind2;
            cm(ind,3) = ind3;
            
            ind++;
            cm(ind,1) = ind1;
            cm(ind,2) = ind3;
            cm(ind,3) = ind4;
          }
      }
      else // i1 = 1 
      {
        for (E_Int k = 0; k < nk1; k++)
#pragma omp for
          for (E_Int j = 0; j < nj1; j++)
          {
            ind1 = 1 + j*nil + k*ninj; //(1,j,k)
            ind2 = ind1 + nil;         //(1,j+1,k)
            ind3 = ind2 + ninj;        //(1,j+1,k+1)
            ind4 = ind3 - nil;         //(1,j,k+1)
            
            ind = 2*(j + k*nj1);
            cm(ind,1) = ind1;
            cm(ind,2) = ind2;
            cm(ind,3) = ind3;
            
            ind++;
            cm(ind,1) = ind1;
            cm(ind,2) = ind3;
            cm(ind,3) = ind4;
          }
      }// i1 = 1
    }
    else
    {
      E_Int ni1nj1 = ni1*nj1;
      E_Int sum, ind, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
      for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
#pragma omp for
          for (E_Int i = 0; i < ni1; i++)
          {
            sum = i+j+k;
            ind1 = 1 + i + j*nil + k*ninj; //A(  i,  j,k)
            ind2 = ind1 + 1;               //B(i+1,  j,k)
            ind3 = ind2 + nil;             //C(i+1,j+1,k)
            ind4 = ind3 - 1;               //D(  i,j+1,k)
            ind5 = ind1 + ninj;            //E(  i,  j,k+1)
            ind6 = ind2 + ninj;            //F(i+1,  j,k+1)
            ind7 = ind3 + ninj;            //G(i+1,j+1,k+1)
            ind8 = ind4 + ninj;            //H(  i,j+1,k+1) 
            
            if (sum%2 == 0) // pair 
            {
              // tetra ABDE
              ind = 5*(i+j*ni1+k*ni1nj1);
              cm(ind,1) = ind1;
              cm(ind,2) = ind2;
              cm(ind,3) = ind4;
              cm(ind,4) = ind5;
              
              // tetra BCDG
              ind++;
              cm(ind,1) = ind2;
              cm(ind,2) = ind3;
              cm(ind,3) = ind4;
              cm(ind,4) = ind7;
              
              // tetra DEGH
              ind++;
              cm(ind,1) = ind4;
              cm(ind,2) = ind5;
              cm(ind,3) = ind7;
              cm(ind,4) = ind8;
              
              // tetra BEFG
              ind++;
              cm(ind,1) = ind2;
              cm(ind,2) = ind5;
              cm(ind,3) = ind6;
              cm(ind,4) = ind7;
              
              // tetra BDEG
              ind++;
              cm(ind,1) = ind2;
              cm(ind,2) = ind4;
              cm(ind,3) = ind5;
              cm(ind,4) = ind7;
            }
            else // impair 
            {
              // tetra ACDH: 1348
              ind = 5*(i+j*ni1+k*ni1nj1);
              cm(ind,1) = ind1;
              cm(ind,2) = ind3;
              cm(ind,3) = ind4;
              cm(ind,4) = ind8;
              
              // tetra AFBC: 1623
              ind++;
              cm(ind,1) = ind1;
              cm(ind,2) = ind2;
              cm(ind,3) = ind3;
              cm(ind,4) = ind6;
              //cm(ind,1) = ind1;
              //cm(ind,2) = ind6;
              //cm(ind,3) = ind2;
              //cm(ind,4) = ind3;

              // tetra HFGC: 8673
              ind++;
              cm(ind,1) = ind3;
              cm(ind,2) = ind6;
              cm(ind,3) = ind7;
              cm(ind,4) = ind8;
              //cm(ind,1) = ind8;
              //cm(ind,2) = ind6;
              //cm(ind,3) = ind7;
              //cm(ind,4) = ind3;

              // tetra FHAE: 6815
              ind++;
              cm(ind,1) = ind6;
              cm(ind,2) = ind5;
              cm(ind,3) = ind8;
              cm(ind,4) = ind1;
              //cm(ind,1) = ind6;
              //cm(ind,2) = ind8;
              //cm(ind,3) = ind1;
              //cm(ind,4) = ind5;
            
              // tetra FHAC: 6813
              ind++;
              cm(ind,1) = ind6;
              cm(ind,2) = ind3;
              cm(ind,3) = ind1;
              cm(ind,4) = ind8;
              //cm(ind,1) = ind6;
              //cm(ind,2) = ind8;
              //cm(ind,3) = ind1;
              //cm(ind,4) = ind3;
            }
          }
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
  
  // Building numpy array
  RELEASESHAREDB(res, array, f, cnl);
  RELEASESHAREDU(tpl, f2, cnl2);
  return tpl;
}
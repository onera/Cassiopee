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

# include "converter.h"

using namespace K_FLD;

// ============================================================================
/* Convert  structured array to a tetraedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertStruct2Tetra(const char* varString, FldArrayF* f,
                                           E_Int ni, E_Int nj, E_Int nk)
{
  // Get dimensionality
  E_Int dim = 3;
  if (ni == 1)
  {
    if (nj == 1 || nk == 1) dim = 1;
    else dim = 2;
  }
  else if (nj == 1)
  {
    if (ni == 1 || nk == 1) dim = 1;
    else dim = 2;
  }
  else if (nk == 1)
  {
    if (ni == 1 || nj == 1) dim = 1;
    else dim = 2;
  }

  // Build the unstructured mesh
  E_Int ni1 = K_FUNC::E_max(1, ni-1);
  E_Int nj1 = K_FUNC::E_max(1, nj-1);
  E_Int nk1 = K_FUNC::E_max(1, nk-1);
  E_Int ninj = ni*nj;
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
  E_Int npts = f->getSize(), api = f->getApi(), nfld = f->getNfld();

  const char* eltType2;
  E_Int nelts; // nb d'elements non structures

  if (dim == 1) { eltType2 = "BAR"; nelts = ncells; }
  else if (dim == 2) { eltType2 = "TRI"; nelts = 2*ncells; }
  else { eltType2 = "TETRA"; nelts = 5*ncells; }

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts,
                                       eltType2, false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);
  FldArrayI& cm2 = *(cn2->getConnect(0));
                
  #pragma omp parallel
  {
    if (dim == 1)
    {
      E_Int ind1, ind2;
      if (nk1 == 1 && nj1 == 1)
      {
        #pragma omp for
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = i + 1;
          ind2 = ind1 + 1;
          cm2(i,1) = ind1;
          cm2(i,2) = ind2;
        }
      }
      else if (ni1 == 1 && nj1 == 1)
      {
        #pragma omp for
        for (E_Int k = 0; k < nk1; k++)
        {
          ind1 = k*ni*nj + 1;
          ind2 = ind1 + ni*nj;
          cm2(k,1) = ind1;
          cm2(k,2) = ind2;
        }
      }
      else if (ni1 == 1 && nk1 == 1)
      {
        #pragma omp for
        for (E_Int j = 0; j < nj1; j++)
        {
          ind1 = j*ni + 1;
          ind2 = ind1 + ni;
          cm2(j,1) = ind1;
          cm2(j,2) = ind2;
        }
      }
    }
    else if (dim == 2)
    {
      E_Int ind, ind1, ind2, ind3, ind4;
      if (nk1 == 1)
      {
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          //starts from 1
          ind1 = i + j*ni + 1;  // (i,j,1)
          ind2 = ind1 + 1;      // (i+1,j,1)
          ind3 = ind2 + ni;     // (i+1,j+1,1)
          ind4 = ind3 - 1;      // (i,j+1,1)
          
          ind = 2*(i + j*ni1);
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind3;
          
          ind++;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
        }
      }
      else if (nj1 == 1)
      {
        #pragma omp for collapse(2)
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = i + k*ninj + 1;  // (i,1,k)
          ind2 = ind1 + ninj;     // (i,1,k+1)
          ind3 = ind2 + 1;        // (i+1,1,k+1)
          ind4 = ind1 + 1;        // (i+1,1,k)
          
          ind = 2*(i + k*ni1);
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind3;
          
          ind++;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
        }
      }
      else // i1 = 1 
      {
        #pragma omp for collapse(2)
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        {
          ind1 = 1 + j*ni + k*ninj;  // (1,j,k)
          ind2 = ind1 + ni;          // (1,j+1,k)
          ind3 = ind2 + ninj;        // (1,j+1,k+1)
          ind4 = ind3 - ni;          // (1,j,k+1)
          
          ind = 2*(j + k*nj1);
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind3;
          
          ind++;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
        }
      }// i1 = 1
    }
    else
    {
      E_Int ni1nj1 = ni1*nj1;
      E_Int sum, ind, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
      #pragma omp for collapse(3)
      for (E_Int k = 0; k < nk1; k++)
      for (E_Int j = 0; j < nj1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        sum = i + j + k;
        ind1 = 1 + i + j*ni + k*ninj;  // A(  i,  j,k)
        ind2 = ind1 + 1;               // B(i+1,  j,k)
        ind3 = ind2 + ni;              // C(i+1,j+1,k)
        ind4 = ind3 - 1;               // D(  i,j+1,k)
        ind5 = ind1 + ninj;            // E(  i,  j,k+1)
        ind6 = ind2 + ninj;            // F(i+1,  j,k+1)
        ind7 = ind3 + ninj;            // G(i+1,j+1,k+1)
        ind8 = ind4 + ninj;            // H(  i,j+1,k+1) 
        
        if (sum%2 == 0) // pair 
        {
          // tetra ABDE
          ind = 5*(i+j*ni1+k*ni1nj1);
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind4;
          cm2(ind,4) = ind5;
          
          // tetra BCDG
          ind++;
          cm2(ind,1) = ind2;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
          cm2(ind,4) = ind7;
          
          // tetra DEGH
          ind++;
          cm2(ind,1) = ind4;
          cm2(ind,2) = ind5;
          cm2(ind,3) = ind7;
          cm2(ind,4) = ind8;
          
          // tetra BEFG
          ind++;
          cm2(ind,1) = ind2;
          cm2(ind,2) = ind5;
          cm2(ind,3) = ind6;
          cm2(ind,4) = ind7;
          
          // tetra BDEG
          ind++;
          cm2(ind,1) = ind2;
          cm2(ind,2) = ind4;
          cm2(ind,3) = ind5;
          cm2(ind,4) = ind7;
        }
        else // impair 
        {
          // tetra ACDH: 1348
          ind = 5*(i+j*ni1+k*ni1nj1);
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind3;
          cm2(ind,3) = ind4;
          cm2(ind,4) = ind8;
          
          // tetra AFBC: 1623
          ind++;
          cm2(ind,1) = ind1;
          cm2(ind,2) = ind2;
          cm2(ind,3) = ind3;
          cm2(ind,4) = ind6;
          //cm2(ind,1) = ind1;
          //cm2(ind,2) = ind6;
          //cm2(ind,3) = ind2;
          //cm2(ind,4) = ind3;

          // tetra HFGC: 8673
          ind++;
          cm2(ind,1) = ind3;
          cm2(ind,2) = ind6;
          cm2(ind,3) = ind7;
          cm2(ind,4) = ind8;
          //cm2(ind,1) = ind8;
          //cm2(ind,2) = ind6;
          //cm2(ind,3) = ind7;
          //cm2(ind,4) = ind3;

          // tetra FHAE: 6815
          ind++;
          cm2(ind,1) = ind6;
          cm2(ind,2) = ind5;
          cm2(ind,3) = ind8;
          cm2(ind,4) = ind1;
          //cm2(ind,1) = ind6;
          //cm2(ind,2) = ind8;
          //cm2(ind,3) = ind1;
          //cm2(ind,4) = ind5;
        
          // tetra FHAC: 6813
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
      }
    }

    // Copy fields to f2
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }
  
  RELEASESHAREDU(tpl, f2, cn2);
  return tpl;
}
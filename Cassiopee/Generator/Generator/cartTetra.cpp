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

// create cartesian unstructured mesh
// 1D: BAR, 2D: TRI, 3D: TETRA
# include "generator.h"

#include <string.h>
using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/* Generateur de grille non structuree tetraedrique 
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille 
   IN: ni, nj, nk: nombre de points
   OUT: array definissant le maillage cree.
   Attention : pour le maillage tetraedrique, il faut que les facettes soient
   identiques pour 2 cellules adjacentes : test si premier noeud de la cellule
   est pair ou impair -> 
   impair: ACDH, AFBC, CFGH, AEFH, ACFH
   pair: ABDE, BCDG, DEGH, BEFG,BDGE
   Les tetraedres sont construits de telle sorte que le tetra ABCD forme 
   un triedre direct */
//=============================================================================
PyObject* K_GENERATOR::cartTetra(PyObject* self, PyObject* args)
{
  E_Int ni, nj, nk;
  E_Float xo, yo, zo;
  E_Float hi, hj, hk;
  E_Int api=1;
  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ TIII_ I_, 
                    &xo, &yo, &zo, &hi, &hj, &hk, &ni, &nj, &nk, &api))
  {
    return NULL;
  }

  // check ni, nj, nk
  if (ni < 1 || nj < 1 || nk < 1)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cartTetra: ni, nj, nk must be >=1.");
    return NULL;
  }

  // 1D, 2D or 3D ?
  E_Int dim0 = 3;
  if (ni == 1)
  {
    if (nj == 1 || nk == 1 ) dim0 = 1;
    else dim0 = 2;
  }
  else if (nj == 1)
  {
    if (ni == 1 || nk == 1 ) dim0 = 1;
    else dim0 = 2;
  }
  else if (nk == 1)
  {
    if (ni == 1 || nj == 1) dim0 = 1;
    else dim0 = 2;
  }

  // Create cartesian mesh
  E_Int ninj = ni*nj; E_Int npts = ninj*nk;
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
  char eltType[8]; 
  E_Int nelts = 0;
  if (dim0 == 1) { strcpy(eltType, "BAR"); nelts = ncells; }
  else if (dim0 == 2) { strcpy(eltType, "TRI"); nelts = 2*ncells; }
  else { strcpy(eltType, "TETRA"); nelts = 5*ncells; }

  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", npts, nelts, eltType, false, api);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  K_ARRAY::getFromArray3(tpl, f, cn);

  K_FLD::FldArrayI& cm = *(cn->getConnect(0));
  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  // Build the unstructured mesh (BE connectivity and fields)
#pragma omp parallel if (nelts > __MIN_SIZE_MEAN__)
  {
    E_Int i, j, k, c;

#pragma omp for 
    for (E_Int ind = 0; ind < npts; ind++)
    {
      k = ind/ninj;
      j = (ind-k*ninj)/ni;
      i = ind-j*ni-k*ninj;
      xt[ind] = xo + i * hi;
      yt[ind] = yo + j * hj;
      zt[ind] = zo + k * hk;
    } 
    
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
          ind1 = k*ni*nj + 1;
          ind2 = ind1 + ni*nj;
          cm(k,1) = ind1;
          cm(k,2) = ind2;
        }
      }
      else if (ni1 == 1 && nk1 == 1)
      {
#pragma omp for
        for (E_Int j = 0; j < nj1; j++)
        {
          ind1 = j*ni + 1;
          ind2 = ind1 + ni;
          cm(j,1) = ind1;
          cm(j,2) = ind2;
        }
      }
      // type BAR 
    }
    else if (dim0 == 2)
    {
      E_Int ind1, ind2, ind3, ind4;
      if (nk1 == 1)
      {
#pragma omp for
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            // starts from 1
            ind1 = i + j*ni + 1; //(i,j,1)
            ind2 = ind1 + 1;     //(i+1,j,1)
            ind3 = ind2 + ni;    //(i+1,j+1,1)
            ind4 = ind3 - 1;     //(i,j+1,1)
            c = 2*(j*ni1 + i);

            cm(c,1) = ind1;
            cm(c,2) = ind2;
            cm(c,3) = ind3;
            
            cm(c+1,1) = ind1;
            cm(c+1,2) = ind3;
            cm(c+1,3) = ind4;
          }
      }
      else if (nj1 == 1)
      {
#pragma omp for
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + k*ninj + 1;  //(i,1,k)
            ind2 = ind1 + ninj;     //(i,1,k+1)
            ind3 = ind2 + 1;        //(i+1,1,k+1)
            ind4 = ind3 - 1;        //(i,1,k+1)
            c = 2*(k*ni1 + i);
                
            cm(c,1) = ind1;
            cm(c,2) = ind2;
            cm(c,3) = ind3;
            
            cm(c+1,1) = ind1;
            cm(c+1,2) = ind3;
            cm(c+1,3) = ind4;
          }
      }
      else // ni1 = 1 
      {
#pragma omp for
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int k = 0; k < nk1; k++)
          {
            ind1 = 1 + j*ni + k*ninj; //(1,j,k)
            ind2 = ind1 + ni;         //(1,j+1,k)
            ind3 = ind2 + ninj;       //(1,j+1,k+1)
            ind4 = ind3 - ni;         //(1,j,k+1)
            c = 2*(j*nk1 + k);
                
            cm(c,1) = ind1;
            cm(c,2) = ind2;
            cm(c,3) = ind3;
            
            cm(c+1,1) = ind1;
            cm(c+1,2) = ind3;
            cm(c+1,3) = ind4;
          }
      }// ni1 = 1
    }
    else
    {
      E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, sum;
#pragma omp for
      for(E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            sum = i+j+k;
            ind1 = 1 + i + j*ni + k*ninj; //A(  i,  j,k)
            ind2 = ind1 + 1;              //B(i+1,  j,k)
            ind3 = ind2 + ni;             //C(i+1,j+1,k)
            ind4 = ind3 - 1;              //D(  i,j+1,k)
            ind5 = ind1 + ninj;           //E(  i,  j,k+1)
            ind6 = ind2 + ninj;           //F(i+1,  j,k+1)
            ind7 = ind3 + ninj;           //G(i+1,j+1,k+1)
            ind8 = ind4 + ninj;           //H(  i,j+1,k+1) 
            c = 5*((k*nj1 + j)*ni1 + i);
            
            if (sum%2 == 0) // pair 
            {
              //tetra ABDE
              cm(c,1) = ind1;
              cm(c,2) = ind2;
              cm(c,3) = ind4;
              cm(c,4) = ind5;
              
              //tetra BCDG
              cm(c+1,1) = ind2;
              cm(c+1,2) = ind3;
              cm(c+1,3) = ind4;
              cm(c+1,4) = ind7;
              
              //tetra DEGH
              cm(c+2,1) = ind4;
              cm(c+2,2) = ind5;
              cm(c+2,3) = ind7;
              cm(c+2,4) = ind8;
              
              //tetra BEFG
              cm(c+3,1) = ind2;
              cm(c+3,2) = ind5;
              cm(c+3,3) = ind6;
              cm(c+3,4) = ind7;
              
              //tetra BDEG
              cm(c+4,1) = ind2;
              cm(c+4,2) = ind4;
              cm(c+4,3) = ind5;
              cm(c+4,4) = ind7;   
            }
            else // impair 
            {
              //tetra ACDH : 1348
              cm(c,1) = ind1;
              cm(c,2) = ind3;
              cm(c,3) = ind4;
              cm(c,4) = ind8;
              
              //tetra AFBC : 1623
              cm(c+1,1) = ind1;
              cm(c+1,2) = ind6;
              cm(c+1,3) = ind2;
              cm(c+1,4) = ind3;
              
              //tetra HFGC : 8763
              cm(c+2,1) = ind8;
              cm(c+2,2) = ind7;
              cm(c+2,3) = ind6;
              cm(c+2,4) = ind3;
              
              //tetra FHAE : 6815
              cm(c+3,1) = ind6;
              cm(c+3,2) = ind8;
              cm(c+3,3) = ind1;
              cm(c+3,4) = ind5;
              
              //tetra FHAC : 6183
              cm(c+4,1) = ind6;
              cm(c+4,2) = ind1;
              cm(c+4,3) = ind8;
              cm(c+4,4) = ind3;
            }
          }
    }
  }

  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}

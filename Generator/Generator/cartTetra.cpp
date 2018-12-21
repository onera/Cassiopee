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
  if (!PYPARSETUPLE(args, 
                    "(ddd)(ddd)(lll)l", "(ddd)(ddd)(iii)i",
                    "(fff)(fff)(lll)l", "(fff)(fff)(iii)i", 
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
  E_Int i, j, k, ind;
  E_Int ninj = ni*nj; E_Int npts = ninj*nk;
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
  char eltType[8]; 
  E_Int nelts=0;
  if (dim0 == 1) { strcpy(eltType, "BAR"); nelts = ncells; }
  else if (dim0 == 2) { strcpy(eltType, "TRI"); nelts = 2*ncells; }
  else { strcpy(eltType, "TETRA"); nelts = 5*ncells; }

  PyObject* tpl = K_ARRAY::buildArray2(3, "x,y,z", npts, nelts, -1, eltType, 0, 0, 0, 0, api);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char*eltType2;
  K_ARRAY::getFromArray2(tpl, varString, f, ni, nj, nk, cn, eltType2);
  E_Int stride = cn->getStride();

  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

#pragma omp parallel for default(shared) private(k,j,i,ind)
  for (ind = 0; ind < npts; ind++)
  {
    k = ind/ninj;
    j = (ind-k*ninj)/ni;
    i = ind-j*ni-k*ninj;
    xt[ind] = xo + i * hi;
    yt[ind] = yo + j * hj;
    zt[ind] = zo + k * hk;
  } 

  // Build the unstructured mesh
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
  E_Int c = 0;
  E_Int* cn1 = NULL;
  E_Int* cn2 = NULL;
  E_Int* cn3 = NULL;
  E_Int* cn4 = NULL;
  switch (dim0)
  {
    case 1:
      cn1 = cn->begin(1);
      cn2 = cn->begin(2);
      if (nk1 == 1 && nj1 == 1)
      {
        for (E_Int i = 0; i < ni1; i++)
        {
          ind1 = i + 1;
          ind2 = ind1 + 1;
          cn1[c] = ind1;
          cn2[c]= ind2;
          c += stride;
        }
      }
      else if (ni1 == 1 && nj1 == 1)
      {
        for (E_Int k = 0; k < nk1; k++)
        {
          ind1 = k*ni*nj + 1;
          ind2 = ind1 + ni*nj;
          cn1[c] = ind1;
          cn2[c] = ind2;
          c += stride;
        }
      }
      else if (ni1 == 1 && nk1 == 1)
      {
        for (E_Int j = 0; j < nj1; j++)
        {
          ind1 = j*ni + 1;
          ind2 = ind1 + ni;
          cn1[c] = ind1;
          cn2[c] = ind2;
          c += stride;
        }
      }
      // type BAR 
      break;

    case 2:
      cn1 = cn->begin(1);
      cn2 = cn->begin(2);
      cn3 = cn->begin(3);
      if (nk1 == 1)
      {
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            // starts from 1
            ind1 = i + j*ni + 1; //(i,j,1)
            ind2 = ind1 + 1;  //(i+1,j,1)
            ind3 = ind2 + ni;//(i+1,j+1,1)
            ind4 = ind3 - 1; //(i,j+1,1)
            
            cn1[c] = ind1;
            cn2[c] = ind2;
            cn3[c] = ind3;
            c += stride;
            
            cn1[c] = ind1;
            cn2[c] = ind3;
            cn3[c] = ind4;
            c += stride;
          }
      }
      else if (nj1 == 1)
      {
        for (E_Int k = 0; k < nk1; k++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = i + k*ninj + 1;  //(i,1,k)
            ind2 = ind1 + ninj; //(i,1,k+1)
            ind3 = ind2 + 1;    //(i+1,1,k+1)
            ind4 = ind3 - 1;    //(i,1,k+1)
                
            cn1[c] = ind1;
            cn2[c] = ind2;
            cn3[c] = ind3;
            c += stride;
            
            cn1[c] = ind1;
            cn2[c] = ind3;
            cn3[c] = ind4;
            c += stride;
          }
      }
      else // ni1 = 1 
      {
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int k = 0; k < nk1; k++)
          {
            ind1 = 1 + j*ni + k*ninj; //(1,j,k)
            ind2 = ind1 + ni;  //(1,j+1,k)
            ind3 = ind2 + ninj;//(1,j+1,k+1)
            ind4 = ind3 - ni;   //(1,j,k+1)
                
            cn1[c] = ind1;
            cn2[c] = ind2;
            cn3[c] = ind3;
            c += stride;
            
            cn1[c] = ind1;
            cn2[c] = ind3;
            cn3[c] = ind4;
            c += stride;
          }
      }// ni1 = 1
      break;
  
    case 3:
      cn1 = cn->begin(1);
      cn2 = cn->begin(2);
      cn3 = cn->begin(3);
      cn4 = cn->begin(4);
    
      E_Int sum;
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
            
            if (sum%2 == 0) // pair 
            {
              //tetra ABDE
              cn1[c] = ind1;
              cn2[c] = ind2;
              cn3[c] = ind4;
              cn4[c] = ind5;
              c += stride;
              
              //tetra BCDG
              cn1[c] = ind2;
              cn2[c] = ind3;
              cn3[c] = ind4;
              cn4[c] = ind7;
              c += stride;
              
              //tetra DEGH
              cn1[c] = ind4;
              cn2[c] = ind5;
              cn3[c] = ind7;
              cn4[c] = ind8;
              c += stride;
              
              //tetra BEFG
              cn1[c] = ind2;
              cn2[c] = ind5;
              cn3[c] = ind6;
              cn4[c] = ind7;
              c += stride;
              
              //tetra BDEG
              cn1[c] = ind2;
              cn2[c] = ind4;
              cn3[c] = ind5;
              cn4[c] = ind7;
              c += stride;    
            }
            else // impair 
            {
              //tetra ACDH : 1348
              cn1[c] = ind1;
              cn2[c] = ind3;
              cn3[c] = ind4;
              cn4[c] = ind8;
              c += stride;
              
              //tetra AFBC : 1623
              cn1[c] = ind1;
              cn2[c] = ind6;
              cn3[c] = ind2;
              cn4[c] = ind3;
              c += stride;
              
              //tetra HFGC : 8763
              cn1[c] = ind8;
              cn2[c] = ind7;
              cn3[c] = ind6;
              cn4[c] = ind3;
              c += stride;
              
              //tetra FHAE : 6815
              cn1[c] = ind6;
              cn2[c] = ind8;
              cn3[c] = ind1;
              cn4[c] = ind5;
              c += stride;
              
              //tetra FHAC : 6183
              cn1[c] = ind6;
              cn2[c] = ind1;
              cn3[c] = ind8;
              cn4[c] = ind3;
              c += stride;
            }
          }
      break;
  }
  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}

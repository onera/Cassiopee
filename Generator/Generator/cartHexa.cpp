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

// Generateur non structure hexaedrique a partir d'une grille cartesienne
// 1D: BAR, 2D: QUAD, 3D: HEXA

#include "generator.h"

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/* Generateur de grille non structuree hexaedrique 
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille 
   IN: ni, nj, nk: nombre de points
   OUT: array definissant le maillage cree. */
//=============================================================================
PyObject* K_GENERATOR::cartHexa(PyObject* self,  PyObject* args)
{
  E_Int ni, nj, nk;
  E_Float xo, yo, zo;
  E_Float hi, hj, hk;
  E_Int api = 1;
  if (!PYPARSETUPLE(args, 
                    "(ddd)(ddd)(lll)l", "(ddd)(ddd)(iii)i", 
                    "(fff)(fff)(lll)l", "(fff)(fff)(iii)i",
                    &xo, &yo, &zo, &hi, &hj, &hk, &ni, &nj, &nk, &api))
  {
    return NULL;
  }
  
  // Check ni, nj, nk
  if (ni < 1 || nj < 1 || nk < 1)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cartHexa: ni, nj, nk must be >=1.");
    return NULL;
  }
  
  // 1D, 2D or 3D ?
  E_Int dim0 = 3;
  if (ni == 1)
  {
    if (nj == 1 || nk == 1) dim0 = 1;
    else dim0 = 2;
  }
  else if (nj == 1)
  {
    if (ni == 1 || nk == 1) dim0 = 1;
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
  if (dim0 == 1) { strcpy(eltType, "BAR"); }
  else if (dim0 == 2) { strcpy(eltType, "QUAD"); }
  else { strcpy(eltType, "HEXA"); }

  PyObject* tpl = K_ARRAY::buildArray2(3, "x,y,z", npts, ncells, -1, eltType, 0, 0, 0, 0, api);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  K_ARRAY::getFromArray2(tpl, f, cn);
  
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
  E_Int ind1, ind2, ind3, ind4;
  E_Int c = 0;
  E_Int* cn1 = NULL;
  E_Int* cn2 = NULL;
  E_Int* cn3 = NULL;
  E_Int* cn4 = NULL;
  E_Int* cn5 = NULL;
  E_Int* cn6 = NULL;
  E_Int* cn7 = NULL;
  E_Int* cn8 = NULL;

  switch (dim0)
  {
    case 1:
    cn1 = cn->begin(1);
    cn2 = cn->begin(2);
   
    if (nk1 == 1 && nj1 == 1)
    {
      for (E_Int i = 0; i < ni1; i++)
      {
        cn1[c] = i+1;
        cn2[c] = i+2;
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
    break;
    
    case 2:
      cn1 = cn->begin(1);
      cn2 = cn->begin(2);
      cn3 = cn->begin(3);
      cn4 = cn->begin(4);
      cn->setAllValuesAtNull();
      if (nk1 == 1)
      {
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            //starts from 1
            ind1 = i + j*ni + 1; //(i,j,1)
            ind2 = ind1 + 1;  //(i+1,j,1)
            ind3 = ind2 + ni; //(i+1,j+1,1)
            ind4 = ind3 - 1; //(i,j+1,1)
            
            cn1[c] = ind1;
            cn2[c] = ind2;
            cn3[c] = ind3;
            cn4[c] = ind4;
            c += stride;
          }
          //for (E_Int i = 0; i < ni1*nj1; i=i+4) printf("%d %d %d %d\n", cn1[i], cn1[i+1], cn1[i+2], cn1[i+3]);
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
            cn4[c] = ind4;
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
            cn4[c] = ind4;
            c += stride;
          }
      }// ni1 = 1
      break;

    default: //3
      cn1 = cn->begin(1);
      cn2 = cn->begin(2);
      cn3 = cn->begin(3);
      cn4 = cn->begin(4);
      cn5 = cn->begin(5);
      cn6 = cn->begin(6);
      cn7 = cn->begin(7);
      cn8 = cn->begin(8);
      for(E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = 1 + i + j*ni + k*ninj; //(i,j,k)
            ind2 = ind1 + 1; // (i+1,j,k)
            ind3 = ind2 + ni; // (i+1,j+1,k)
            ind4 = ind3 -1; // (i,j+1,k)
            cn1[c] = ind1;
            cn2[c] = ind2;
            cn3[c] = ind3;
            cn4[c] = ind4;
            cn5[c] = ind1 + ninj;
            cn6[c] = ind2 + ninj;
            cn7[c] = ind3 + ninj;
            cn8[c] = ind4 + ninj;
            c += stride;
          }
      break;
  }
  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}

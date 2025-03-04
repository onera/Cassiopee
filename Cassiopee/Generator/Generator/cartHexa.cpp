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
  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ TIII_ I_, 
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
  E_Int ninj = ni*nj; E_Int npts = ninj*nk;
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees

  char eltType[8];
  if (dim0 == 1) { strcpy(eltType, "BAR"); }
  else if (dim0 == 2) { strcpy(eltType, "QUAD"); }
  else { strcpy(eltType, "HEXA"); }

  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", npts, ncells, eltType, false, api);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  K_ARRAY::getFromArray3(tpl, f, cn);

  K_FLD::FldArrayI& cm = *(cn->getConnect(0));
  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  // Build the unstructured mesh (BE connectivity and fields)
#pragma omp parallel if (ncells > __MIN_SIZE_MEAN__)
  {
    E_Int i, j, k, c;
    E_Int ind1, ind2, ind3, ind4;
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
      if (nk1 == 1 && nj1 == 1)
      {
#pragma omp for
        for (E_Int i = 0; i < ni1; i++)
        {
          cm(i,1) = i+1;
          cm(i,2) = i+2;
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
    }
    else if (dim0 == 2)
    {
      if (nk1 == 1)
      {
#pragma omp for
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            //starts from 1
            ind1 = i + j*ni + 1; //(i,j,1)
            ind2 = ind1 + 1;     //(i+1,j,1)
            ind3 = ind2 + ni;    //(i+1,j+1,1)
            ind4 = ind3 - 1;     //(i,j+1,1)
            c = j*ni1 + i;
            cm(c,1) = ind1;
            cm(c,2) = ind2;
            cm(c,3) = ind3;
            cm(c,4) = ind4;
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
            c = k*ni1 + i;
            cm(c,1) = ind1;
            cm(c,2) = ind2;
            cm(c,3) = ind3;
            cm(c,4) = ind4;
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
            c = j*nk1 + k;
            cm(c,1) = ind1;
            cm(c,2) = ind2;
            cm(c,3) = ind3;
            cm(c,4) = ind4;
          }
      }// ni1 = 1
    }
    else
    {
#pragma omp for
      for(E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
          for (E_Int i = 0; i < ni1; i++)
          {
            ind1 = 1 + i + j*ni + k*ninj; //(i,j,k)
            ind2 = ind1 + 1;              // (i+1,j,k)
            ind3 = ind2 + ni;             // (i+1,j+1,k)
            ind4 = ind3 -1;               // (i,j+1,k)
            c = (k*nj1 + j)*ni1 + i;
            cm(c,1) = ind1;
            cm(c,2) = ind2;
            cm(c,3) = ind3;
            cm(c,4) = ind4;
            cm(c,5) = ind1 + ninj;
            cm(c,6) = ind2 + ninj;
            cm(c,7) = ind3 + ninj;
            cm(c,8) = ind4 + ninj;
          }
    }
  }

  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}
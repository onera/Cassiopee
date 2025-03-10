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

// Grille cartesienne avec facteurs d'expansion

#include "generator.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC; 

// ============================================================================
/* Create a cartesian mesh of nixnjxnk points with expension factors
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille. Enter negative values for extrusion in -x,-y or -z direction
   IN: ni, nj, nk: nombre de points
   IN: ri, rj, rk: facteur d'expansion dans chaque direction
   IN: doubleLeft=(1,1,1), force un double pas h a gauche
   IN: doubleRigh=(1,1,1), force un double pas h a droite
   OUT: array definissant le maillage cree. */
// ============================================================================
PyObject* K_GENERATOR::cartr1(PyObject* self, PyObject* args)
{
  E_Int ni, nj, nk;
  E_Float xo, yo, zo;
  E_Float hi, hj, hk;
  E_Float ri, rj, rk;
  E_Int doubleLefti, doubleRighti, doubleLeftj, doubleRightj, doubleLeftk, doubleRightk; 
  E_Int api = 1;
  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ TRRR_ TIII_ TIII_ TIII_ I_, 
                    &xo, &yo, &zo, &hi, &hj, &hk, &ri, &rj, &rk, &ni, &nj, &nk,  
                    &doubleLefti, &doubleLeftj, &doubleLeftk, &doubleRighti, 
                    &doubleRightj, &doubleRightk, &api))
  {
    return NULL;
  }
  if (ni < 1 || nj < 1 || nk < 1)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cart: ni, nj, nk must be >= 1.");
    return NULL;
  }

  E_Int i, j, k, ind;
  // Create cartesian mesh
  PyObject* tpl;
  tpl = K_ARRAY::buildArray2(3, "x,y,z", ni, nj, nk, api);
  
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  K_ARRAY::getFromArray2(tpl, varString, f, ni, nj, nk, c, eltType);

  E_Int nij = ni*nj;
  E_Int nijk = ni*nj*nk;

  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);
  
  if (K_FUNC::fEqual(ri, 1.0) == true)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      xt[ind] = xo + hi * i; 
    }   
  }   
  else if (doubleLefti == 0 && doubleRighti == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      xt[ind] = xo + hi * ( ( -1. + pow(ri, i) ) / (-1. + ri) ); 
    }
  }
  else if (doubleLefti == 1 && doubleRighti == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (i == 0) xt[ind] = xo;
      else xt[ind] = xo + hi + hi * ((-1. + pow(ri, (i-1))) / (-1. + ri));
    }
  }
  else if (doubleLefti == 0 && doubleRighti == 1)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (i == ni-1) xt[ind] = xo + hi*pow(ri,ni-3) + hi*((-1.+pow(ri,ni-2))/(-1.+ri));
      else xt[ind] = xo + hi*((-1.+pow(ri,i)) / (-1.+ri));
    }
  }
  else // all double 
  {
     #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (i == 0) xt[ind] = xo;
      else if (i == ni-1) xt[ind] = xo + hi + hi*pow(ri,ni-4) + hi*((-1.+pow(ri,ni-3))/(-1.+ri));
      else xt[ind] = xo + hi + hi*((-1.+pow(ri,(i-1))) / (-1.+ri));
    } 
  }

  if (K_FUNC::fEqual(rj, 1.0) == true)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      yt[ind] = yo + hj * j; 
    }   
  }   
  else if (doubleLeftj == 0 && doubleRightj == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      yt[ind] = yo + hj * ( ( -1. + pow(rj, j) ) / (-1. + rj) ); 
    }
  }
  else if (doubleLeftj == 1 && doubleRightj == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (j == 0) yt[ind] = yo;
      else yt[ind] = yo + hj + hj * ((-1. + pow(rj, (j-1))) / (-1. + rj));
    }
  }
  else if (doubleLeftj == 0 && doubleRightj == 1)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (j == nj-1) yt[ind] = yo + hj*pow(rj,nj-3) + hj*((-1.+pow(rj,nj-2))/(-1.+rj));
      else yt[ind] = yo + hj*((-1.+pow(rj,j)) / (-1.+rj));
    }
  }
  else // all double 
  {
     #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (j == 0) yt[ind] = yo;
      else if (j == nj-1) yt[ind] = yo + hj + hj*pow(rj,nj-4) + hj*((-1.+pow(rj,nj-3))/(-1.+rj));
      else yt[ind] = yo + hj + hj*((-1.+pow(rj,(j-1))) / (-1.+rj));
    } 
  }

  if (K_FUNC::fEqual(rk, 1.0) == true)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      zt[ind] = zo + hk * k; 
    }   
  }   
  else if (doubleLeftk == 0 && doubleRightk == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      zt[ind] = zo + hk * ( ( -1. + pow(rk, k) ) / (-1. + rk) ); 
    }
  }
  else if (doubleLeftk == 1 && doubleRightk == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (k == 0) zt[ind] = zo;
      else zt[ind] = zo + hk + hk * ((-1. + pow(rk, (k-1))) / (-1. + rk));
    }
  }
  else if (doubleLeftk == 0 && doubleRightk == 1)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (k == nk-1) zt[ind] = zo + hk*pow(rk,nk-3) + hk*((-1.+pow(rk,nk-2))/(-1.+rk));
      else zt[ind] = zo + hk*((-1.+pow(rk,k)) / (-1.+rk));
    }
  }
  else // all double 
  {
     #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (k == 0) zt[ind] = zo;
      else if (k == nk-1) zt[ind] = zo + hk + hk*pow(rk,nk-4) + hk*((-1.+pow(rk,nk-3))/(-1.+rk));
      else zt[ind] = zo + hk + hk*((-1.+pow(rk,(k-1))) / (-1.+rk));
    } 
  }

  // Return array
  RELEASESHAREDS(tpl, f);
  return tpl;
}

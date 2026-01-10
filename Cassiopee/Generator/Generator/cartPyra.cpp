/*    
    Copyright 2013-2026 ONERA.

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

// Generateur non structure pyra a partir d'une grille cartesienne 3D
# include "generator.h"

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/* Generateur de grille non structuree PYRA
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille 
   IN: ni, nj, nk: nombre de points
   OUT: array definissant le maillage cree. */
//=============================================================================
PyObject* K_GENERATOR::cartPyra(PyObject* self, PyObject* args)
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
  
  // Check ni, nj, nk
  if (ni < 2 || nj < 2 || nk < 2)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cartPyra: ni, nj, nk must be > 1.");
    return NULL;
  }
  
  E_Int ninj = ni*nj;
  E_Int ni1 = E_max(1, E_Int(ni)-1);
  E_Int nj1 = E_max(1, E_Int(nj)-1);
  E_Int nk1 = E_max(1, E_Int(nk)-1);  
  E_Int ncells = ni1*nj1*nk1; // nb de cellules structurees
  E_Int npts = ninj*nk+ncells;
  E_Int nelts = 6*ncells; // 6 pyra par cellule structuree
  E_Int shift = ninj*nk;

  PyObject* tpl = K_ARRAY::buildArray3(3, "x,y,z", npts, nelts, "PYRA", false, api);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  K_ARRAY::getFromArray3(tpl, f, cn);

  K_FLD::FldArrayI& cm = *(cn->getConnect(0));
  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  // Build the unstructured mesh (BE connectivity and fields)
#pragma omp parallel if (nelts > __MIN_SIZE_MEAN__)
  {
    E_Int ind, indp, indcell, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8, c;
#pragma omp for
  for (E_Int k = 0; k < nk; k++)
    for (E_Int j = 0; j < nj; j++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + j*ni + k*ninj;
        xt[ind] = xo + i*hi;
        yt[ind] = yo + j*hj;
        zt[ind] = zo + k*hk;        
      }

#pragma omp for
  for (E_Int kc = 0; kc < nk1; kc++)
    for (E_Int jc = 0; jc < nj1; jc++)
      for (E_Int ic = 0; ic < ni1; ic++)
      {
        indcell = (kc*nj1 + jc)*ni1 + ic + shift;
        xt[indcell] = xo + (ic+0.5)*hi;
        yt[indcell] = yo + (jc+0.5)*hj;
        zt[indcell] = zo + (kc+0.5)*hk;        
      } 

#pragma omp for
  for (E_Int k = 0; k < nk1; k++)
    for (E_Int j = 0; j < nj1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        indcell = (k*nj1 + j)*ni1 + i;
        ind1 = i+j*ni+k*ninj;
        ind2 = ind1+1;
        ind3 = ind2+ni;
        ind4 = ind1+ni;
        ind5 = ind1+ninj;
        ind6 = ind2+ninj;
        ind7 = ind3+ninj;
        ind8 = ind4+ninj;
        indp = indcell+shift;
        c = 6*indcell;

        //pyra a base 1234
        cm(c,1) = ind1+1;
        cm(c,2) = ind2+1;
        cm(c,3) = ind3+1;
        cm(c,4) = ind4+1;
        cm(c,5) = indp+1;

        //pyra a base 5876
        cm(c+1,1) = ind5+1;
        cm(c+1,2) = ind8+1;
        cm(c+1,3) = ind7+1;
        cm(c+1,4) = ind6+1;
        cm(c+1,5) = indp+1;

        //pyra a base 1485
        cm(c+2,1) = ind1+1;
        cm(c+2,2) = ind4+1;
        cm(c+2,3) = ind8+1;
        cm(c+2,4) = ind5+1;
        cm(c+2,5) = indp+1;

        //pyra a base 2673
        cm(c+3,1) = ind2+1;
        cm(c+3,2) = ind6+1;
        cm(c+3,3) = ind7+1;
        cm(c+3,4) = ind3+1;
        cm(c+3,5) = indp+1;
        
        // pyra a base 1652
        cm(c+4,1) = ind1+1;
        cm(c+4,2) = ind5+1;
        cm(c+4,3) = ind6+1;
        cm(c+4,4) = ind2+1;
        cm(c+4,5) = indp+1;

        //pyra a base4378
        cm(c+5,1) = ind4+1;
        cm(c+5,2) = ind3+1;
        cm(c+5,3) = ind7+1;
        cm(c+5,4) = ind8+1;
        cm(c+5,5) = indp+1;
      }
  }

  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}
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
  if (!PYPARSETUPLE(args, 
                    "(ddd)(ddd)(lll)l", "(ddd)(ddd)(iii)i",
                    "(fff)(fff)(lll)l", "(fff)(fff)(iii)i",
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

  PyObject* tpl = K_ARRAY::buildArray2(3, "x,y,z", npts, nelts, -1, "PYRA", 0, 0, 0, 0, api);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  K_ARRAY::getFromArray2(tpl, f, cn);
  E_Int stride = cn->getStride();

  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);
  E_Int ind, indp, indcell, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
  E_Int shift = ninj*nk;

#pragma omp parallel for default(shared) private(ind) if (nk > 50)
  for (E_Int k = 0; k < nk; k++)
    for (E_Int j = 0; j < nj; j++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + j*ni + k*ninj;
        xt[ind] = xo + i*hi;
        yt[ind] = yo + j*hj;
        zt[ind] = zo + k*hk;        
      }
  for (E_Int kc = 0; kc < nk1; kc++)
    for (E_Int jc = 0; jc < nj1; jc++)
      for (E_Int ic = 0; ic < ni1; ic++)
      {
        indcell = ic + jc*ni1 + kc*ni1*nj1 + shift;
        xt[indcell] = xo + (ic+1./2.)*hi;
        yt[indcell] = yo + (jc+1./2.)*hj;
        zt[indcell] = zo + (kc+1./2.)*hk;        
      } 

  // Build the unstructured mesh
  E_Int c = 0;
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  E_Int* cn4 = cn->begin(4);
  E_Int* cn5 = cn->begin(5);
  for (E_Int k = 0; k < nk1; k++)
    for (E_Int j = 0; j < nj1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        indcell = i + j*ni1 + k*ni1*nj1;
        ind1 = i+j*ni+k*ninj;
        ind2 = ind1+1;
        ind3 = ind2+ni;
        ind4 = ind1+ni;
        ind5 = ind1+ninj;
        ind6 = ind2+ninj;
        ind7 = ind3+ninj;
        ind8 = ind4+ninj;
        indp = indcell+shift;

        //pyra a base 1234
        cn1[c] = ind1+1;
        cn2[c] = ind2+1;
        cn3[c] = ind3+1;
        cn4[c] = ind4+1;
        cn5[c] = indp+1;
        c += stride;

        //pyra a base 5876
        cn1[c] = ind5+1;
        cn2[c] = ind8+1;
        cn3[c] = ind7+1;
        cn4[c] = ind6+1;
        cn5[c] = indp+1;
        c += stride;

        //pyra a base 1485
        cn1[c] = ind1+1;
        cn2[c] = ind4+1;
        cn3[c] = ind8+1;
        cn4[c] = ind5+1;
        cn5[c] = indp+1;
        c += stride;

        //pyra a base 2673
        cn1[c] = ind2+1;
        cn2[c] = ind6+1;
        cn3[c] = ind7+1;
        cn4[c] = ind3+1;
        cn5[c] = indp+1;
        c += stride;
        
        // pyra a base 1652
        cn1[c] = ind1+1;
        cn2[c] = ind5+1;
        cn3[c] = ind6+1;
        cn4[c] = ind2+1;
        cn5[c] = indp+1;
        c += stride;

        //pyra a base4378
        cn1[c] = ind4+1;
        cn2[c] = ind3+1;
        cn3[c] = ind7+1;
        cn4[c] = ind8+1;
        cn5[c] = indp+1;
        c += stride;
      }
  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}

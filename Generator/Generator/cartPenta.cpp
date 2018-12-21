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
// 3D: PENTA: maillage prismatique

# include "generator.h"

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/* Generateur de grille non structuree prismatique
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille 
   IN: ni, nj, nk: nombre de points
   OUT: array definissant le maillage cree.
   On decoupe chaque hexa en deux prismes */
//=============================================================================
PyObject* K_GENERATOR::cartPenta(PyObject* self, PyObject* args)
{
  E_Int ni, nj, nk;
  E_Float xo, yo, zo;
  E_Float hi, hj, hk;
  E_Int api=1;
  if (!PYPARSETUPLE(args, 
                    "(ddd)(ddd)(lll)l", "(ddd)(ddd)(iii)i",
                    "(fff)(fff)(lll)l", "(fff)(fff)(lll)l",
                    &xo, &yo, &zo, &hi, &hj, &hk, &ni, &nj, &nk, &api))
  {
    return NULL;
  }  

  // check ni, nj, nk
  if (ni < 2 || nj < 2 || nk < 2)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cartPenta: array must be 3D.");
    return NULL;
  }
  
  E_Int ninj = ni*nj; E_Int npts = ninj*nk;
  E_Int ni1 = ni-1;
  E_Int nj1 = nj-1;
  E_Int nk1 = nk-1;
  E_Int nhexas = ni1*nj1*nk1;
  E_Int nelts = 2 * nhexas;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
  E_Int i, j, k, ind;
  E_Int c = 0;
  PyObject* tpl = K_ARRAY::buildArray2(3, "x,y,z", npts, nelts, -1, "PENTA", 0, 0, 0, 0, api);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char*eltType2;
  K_ARRAY::getFromArray2(tpl, varString, f, ni, nj, nk, cn, eltType2);
  E_Int stride = cn->getStride();

  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  E_Int* cn4 = cn->begin(4);
  E_Int* cn5 = cn->begin(5);
  E_Int* cn6 = cn->begin(6);

  // Create cartesian mesh
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

  // build the unstructured mesh
  for(E_Int k = 0; k < nk1; k++)
    for (E_Int j = 0; j < nj1; j++)
      for (E_Int i = 0; i < ni1; i++)
      {
        ind1 = 1 + i + j*ni + k*ninj;
        ind2 = ind1 + 1;
        ind3 = ind2 + ni;
        ind4 = ind3 - 1;
        ind5 = ind1 + ninj;
        ind6 = ind2 + ninj;
        ind7 = ind3 + ninj;
        ind8 = ind4 + ninj;

        // prisme1
        cn1[c] = ind1;
        cn2[c] = ind3;
        cn3[c] = ind4;
        cn4[c] = ind5;
        cn5[c] = ind7;
        cn6[c] = ind8;
        c += stride;
        
        // prisme2
        cn1[c] = ind1;
        cn2[c] = ind2;
        cn3[c] = ind3;
        cn4[c] = ind5;
        cn5[c] = ind6;
        cn6[c] = ind7;
        c += stride;
      }  

  RELEASESHAREDU(tpl, f, cn);
  return tpl;
}

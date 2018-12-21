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

// Grille cartesienne reguliere

#include "generator.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC; 

// ============================================================================
/* Create a cartesian mesh of nixnjxnk points 
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas de la grille 
   IN: ni, nj, nk: nombre de points
   OUT: array definissant le maillage cree. */
// ============================================================================
PyObject* K_GENERATOR::cartStruct(PyObject* self, PyObject* args)
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

#pragma omp parallel for default(shared) private(k,j,i,ind)
  for (ind = 0; ind < nijk; ind++)
  {
    k = ind/nij;
    j = (ind-k*nij)/ni;
    i = ind-j*ni-k*nij;
    xt[ind] = xo + i * hi;
    yt[ind] = yo + j * hj;
    zt[ind] = zo + k * hk;
  }

  // Return array
  RELEASESHAREDS(tpl, f);
  return tpl;
}

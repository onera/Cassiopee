/*    
    Copyright 2013-2024 Onera.

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

#include "generator.h"
using namespace K_FLD;

//============================================================================
/* Maillage en chapeau */
//============================================================================
PyObject* K_GENERATOR::pointedHat(PyObject* self, PyObject* args)
{
  PyObject* array;
  double x, y, z;

  if (!PYPARSETUPLE_(args, O_ TRRR_, &array, &x, &y, &z)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1)
  {
    if (res == 2) {delete f; delete cn;}
    PyErr_SetString(PyExc_TypeError,
                    "pointedHat: array must be a structured array.");
    return NULL;
  }
  
  // traitement
  if (im < 1 || jm < 1 || km != 1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                      "pointedHat: array must be a i or a ij-array.");
    return NULL;
  }
  if (f->getNfld() != 3) 
  {
    printf("Warning: pointedHat: only coordinates are considered.\n");
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);    
  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "pointedHat: can't find coordinates in array.");
    return NULL;
  }
  
  posx++; posy++; posz++;
  E_Int npts = im*jm;
  FldArrayF* sol = new FldArrayF(2*npts,3);
  E_Float* xt = sol->begin(1);
  E_Float* yt = sol->begin(2);
  E_Float* zt = sol->begin(3);
  E_Float* xt0 = f->begin(posx);
  E_Float* yt0 = f->begin(posy);
  E_Float* zt0 = f->begin(posz);  

  for (E_Int i = 0; i < npts; i++)
  {
    xt[i] = xt0[i]; yt[i] = yt0[i]; zt[i] = zt0[i];
    xt[i+npts] = x; yt[i+npts] = y; zt[i+npts] = z;
  }
  
  E_Int im2 = im; E_Int jm2 = jm; E_Int km2 = 2;
  if (jm == 1) {jm2 = 2; km2 = 1;}
  PyObject* tpl = K_ARRAY::buildArray(*sol, "x,y,z", im2, jm2, km2);
  delete f; delete sol;
  return tpl;
}

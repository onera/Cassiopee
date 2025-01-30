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

#include "geom.h"
using namespace K_FLD;

// ============================================================================
/* Create a profile line of N points defined by Bernstein polynomials */
// ============================================================================
PyObject* K_GEOM::bezier(PyObject* self, PyObject* args)
{
  E_Int N, M;
  E_Float density;
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_ II_ R_,
                    &array, &N, &M, &density))
  {
      return NULL;
  }
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* et;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, et);  

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "bezier: array is not valid.");
    return NULL;
  }
  if (res == 2)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "bezier: array must be an i-array or a i,j-array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "bezier: coordinates not found in array.");
    return NULL;
  }
  posx++; posy++; posz++;
  
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  if (im != 1)
  {
    if (jm == 1) 
    {
      K_FLD::FldArrayF PF;
      K_COMPGEOM::regularBezier(im, N, density, xt, yt, zt, PF);
      delete f;
      PyObject* tpl = K_ARRAY::buildArray(PF, "x,y,z", PF.getSize(), 1, 1);
      return tpl;
    }
    else 
    {
      K_FLD::FldArrayF PF;
      E_Int niout, njout;
      K_COMPGEOM::regularBezier2D(im, jm, N, M, density, xt, yt, zt, PF,
                                  niout, njout);
      delete f;
      PyObject* tpl = K_ARRAY::buildArray(PF, "x,y,z", niout, njout, 1);
      return tpl;
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "bezier: control points array must be i-array or i,j-array.");
    return  NULL;
  }
}

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

// Information on geometries

# include "geom.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Get the curvature radius of all points */
// ============================================================================
PyObject* K_GEOM::getCurvatureRadius(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getCurvatureRadius: array must be a BAR or an i-array.");
    return NULL;  
  }
  if ( res == 1 && ( jm != 1 || km != 1 || im == 1 ) )
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "getCurvatureRadius: structured array must be an i-array.");
    return NULL;
  }
  else if ( res == 2  )
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "getCurvatureRadius: not for unstructured arrays.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "getCurvatureRadius: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;  
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  E_Int npts = f->getSize();

  // calcul du rayon de courbure
  FldArrayF* radp = new FldArrayF(npts);
  FldArrayF& rad = *radp;
  
  FldArrayF curv(npts);
  K_COMPGEOM::compCurvature(npts, xt, yt, zt, curv);
  E_Float c;
  for (E_Int i = 0; i < npts; i++)
  {
    c = curv[i];
    if ( K_FUNC::fEqualZero(c) && c >= 0)
      rad[i] = K_CONST::E_MAX_FLOAT;
    else if ( K_FUNC::fEqualZero(c) && c <= 0)
      rad[i] = -K_CONST::E_MAX_FLOAT;
    else rad[i] = 1./c;
  }
  PyObject* tpl = K_ARRAY::buildArray(rad, "radius", npts, 1, 1);
  delete f; delete radp;
  return tpl;
}     

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

// Information on geometries

# include <string.h>
# include "geom.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Get the curvature angle of all points */
// ============================================================================
PyObject* K_GEOM::getCurvatureAngle(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float dirVect[3];
  dirVect[0] = 0.;  dirVect[1] = 0.;  dirVect[2] = 1.;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  E_Int posx, posy, posz;

  if (res == 1 || res == 2)
  {
    // data check
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                      "getCurvatureAngle: can't find coordinates in array.");
      return NULL;  
    }
    posx++; posy++; posz++;  

    if (res == 1 && (jm != 1 || km != 1))
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                      "getCurvatureAngle: array must be a TRI, a BAR or an i-array.");
      return NULL;  
    }
    if ( (res == 2 && (strcmp(eltType, "BAR") != 0)))// || 
      if ((res == 2 && (strcmp(eltType, "TRI") != 0))) 
      {
        delete f; delete cn;
        PyErr_SetString(PyExc_TypeError,
                        "getCurvatureAngle: array must be a TRI, a BAR or an i-array.");
        return NULL;  
      }
    
    E_Int sizef = f->getSize();
    E_Int api = f->getApi();
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    PyObject* tpl = NULL;
    
    if (res == 1)
    {
      FldArrayF* an = new FldArrayF(sizef);
      FldArrayF& angle = *an;
      angle.setAllValuesAtNull();
      K_COMPGEOM::compCurvatureAngle(im, xt, yt, zt, dirVect, angle);
      tpl = K_ARRAY::buildArray3(*an, "angle", im, jm, km);
      delete an;
    }
    else if (res == 2 && strcmp(eltType, "BAR") == 0)
    {
      FldArrayF* an = new FldArrayF(sizef);
      FldArrayF& angle = *an;
      angle.setAllValuesAtNull();
      K_COMPGEOM::compCurvatureAngleForBar(sizef, xt, yt, zt, *cn, dirVect, angle);
      tpl = K_ARRAY::buildArray3(*an, "angle", *cn, eltType, api);
      delete an; delete cn;
    }
    else if (res == 2 && strcmp(eltType, "TRI") == 0)
    {
      FldArrayF* an = new FldArrayF(cn->getSize(), cn->getNfld());
      FldArrayF& angle = *an;
      angle.setAllValuesAtNull();
      E_Int ret = K_COMPGEOM::compCurvatureAngleForTri(sizef,
                                                       xt, yt, zt, *cn, angle);
      if (ret == 1)
      {
        tpl = K_ARRAY::buildArray3(*an, "angle1,angle2,angle3", *cn, eltType, api);
        delete an; delete cn;
      }
    }

    delete f;
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getCurvatureAngle: invalid array.");
    return NULL;
  }
}     

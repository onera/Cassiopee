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

# include "transform.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Split a i-array following curvature angle */
// ============================================================================
PyObject* K_TRANSFORM::splitCurvatureAngle( PyObject* self,
                                            PyObject* args )
{
  E_Float tol = 45.;
  E_Float dirVect[3];
  dirVect[0] = 0;  dirVect[1] = 0;  dirVect[2] = 1; 
  PyObject* array;

  if (!PYPARSETUPLEF(args,
                    "Od", "Of",
                    &array, &tol) && 
      !PYPARSETUPLEI(args,
                    "O", "O",
                    &array))
  {
      return NULL;
  }

  // Check arrays
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int is, ic;

  E_Int res = K_ARRAY::getFromArray(array, varString, 
                                    f, im, jm, km, cn, eltType);

  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                      "splitCurvatureAngle: can't find coordinates in array.");
      return NULL;        
    }
    if (im < 2 || jm != 1 || km != 1 )
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                      "splitCurvatureAngle: structured array must be an i-array.");
      return NULL;         
    }
    posx++; posy++; posz++;

    // Compute angle
    FldArrayF angle(im);
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    K_COMPGEOM::compCurvatureAngle(im, xt, yt, zt, dirVect, angle);

    // Split line
    is = 0; ic = 1;
    while (ic < im)
    {
      //printf("%f ", angle[ic]);
      if (angle[ic] < 180.-tol || angle[ic] > 180.+tol)
      {
        // split
        is = ic+1;
        break;
      }
      ic++;
    }
    
    delete f;
    PyObject* tpl = Py_BuildValue("l", long(is));
    return tpl;
  }
  else if (res == 2)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "splitCurvatureAngle: not for unstructured arrays.");
    return NULL;   
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitCurvatureAngle: invalid array.");
    return NULL;
  }
}

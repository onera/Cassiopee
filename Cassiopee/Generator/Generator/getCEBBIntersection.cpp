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

# include "generator.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Retourne 1 si les CEBB des 2 blocs s intersectent*/
//=============================================================================
PyObject* K_GENERATOR::getCEBBIntersectionOfArrays(PyObject* self, PyObject* args)
{
  PyObject* array1;
  PyObject* array2;
  E_Float tol;
  if (!PYPARSETUPLE_(args, OO_ R_, &array1, &array2, &tol)) return NULL;

  // Check array1
  E_Int ni1, nj1, nk1;
  FldArrayF* f1;
  FldArrayI* cn1;
  char* varString1;
  char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray(
    array1, varString1, f1, ni1, nj1, nk1, cn1, eltType1, true);
  
  // Check array2
  E_Int ni2, nj2, nk2;
  FldArrayF* f2;
  FldArrayI* cn2;
  char* varString2;
  char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray(
    array2, varString2, f2, ni2, nj2, nk2, cn2, eltType2, true);
  E_Int isIntersect = 0;

  // cas structure seulement
  if ( res1 == 1 && res2 == 1 )
  {
    // determination de x,y et z ds les arrays
    E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
    E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
    E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
    E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
    E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
    E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
    if (posx1 == -1 || posy1 == -1 || posz1 == -1 || 
        posx2 == -1 || posy2 == -1 || posz2 == -1)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_TypeError,
                      "CEBBIntersection: can't find coordinates in array.");
      return NULL;
    }
    posx1++; posy1++; posz1++; posx2++; posy2++; posz2++;

    // test si les elements cartesiens des deux arrays s'intersectent
    isIntersect = K_COMPGEOM::compCEBBIntersection(
      ni1, nj1, nk1, posx1, posy1, posz1, *f1,
      ni2, nj2, nk2, posx2, posy2, posz2, *f2, tol);

    RELEASESHAREDS(array1, f1); RELEASESHAREDS(array2, f2);
    return Py_BuildValue(I_, isIntersect);
  }
  else if (res1 == 2 || res2 == 2)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);

    PyErr_SetString(PyExc_TypeError,
                    "CEBBIntersection: not for unstructured arrays.");
    return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "CEBBIntersection: unknown type of arrays.");
    return NULL; 
  }
}

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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Conformisation topologique d'un NGON
   Le jeu de faces correspondant a une plus grande face la remplace */
//=============================================================================
PyObject* K_CONVERTER::conformizeNGon(PyObject* self, PyObject* args)
{
  PyObject* array; E_Float tol;
  if (!PYPARSETUPLE_(args, O_ R_, &array, &tol)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, ni, nj, nk, cnl, eltType, true);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "conformizeNGon: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType, "NGON") != 0)
  { RELEASESHAREDU(array, f, cnl); return array; }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (posx == 0 || posy == 0 || posz == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "conformizeNGon: coordinates not found.");
    RELEASESHAREDU(array, f, cnl); return NULL;
  }

  // nouveau code
  FldArrayI* cn;
  conformizeNGon(*f, posx, posy, posz, *cnl, tol, cn);
  
  // Construction de l'array de sortie
  PyObject* tpl = K_ARRAY::buildArray(*f, varString, *cn, 8, "NGON");
  delete cn;
  RELEASESHAREDU(array, f, cnl); 

  return tpl;
}

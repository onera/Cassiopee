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

// Information on grids

# include "generator.h"

using namespace K_FUNC;
using namespace K_FLD;

extern "C"
{
  void k6checkmesh_(const E_Int& ni, const E_Int& nj, const E_Int& nk,
                    const E_Float* x, const E_Float* y, const E_Float* z);
}

// ============================================================================
/* Check mesh for regularity, etc... */
// ============================================================================
PyObject* K_GENERATOR::checkMesh( PyObject* self,
                                  PyObject* args )
{ 
  PyObject* array;
  if ( !PyArg_ParseTuple(args, "O", &array ) ) return NULL;
 
  // Check array
  E_Int im, jm, km;
  FldArrayF* f;
  char* varString;
  char* eltType;
  FldArrayI* cn;

  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                      "check: can't find coordinates in array.");
      Py_INCREF(Py_None);
      return Py_None;        
    }
    posx++; posy++; posz++;
      
    // Check mesh
    k6checkmesh_(im, jm, km, f->begin(posx), f->begin(posy), f->begin(posz));
    delete f;
  }
  else if (res == 2)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "check: not used for unstructured arrays.");
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "check: unrecognized type of array.");
  }
  Py_INCREF(Py_None);
  return Py_None;
}

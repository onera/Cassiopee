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

//=============================================================================
/* Returns the barycenter of an array */
//=============================================================================
PyObject* K_GENERATOR::barycenter(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* weight;
  if (!PyArg_ParseTuple(args, "OO", &array, &weight)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, 
                                     ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "barycenter: invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "barycenter: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float xb, yb, zb;
  if (weight == Py_None)
    K_COMPGEOM::barycenter(f->getSize(), f->begin(posx), f->begin(posy),
                           f->begin(posz), xb, yb, zb);
  else 
  {
    E_Int niw, njw, nkw;
    FldArrayF* fw; FldArrayI* cnw;
    char* varStringw; char* eltTypew;
    E_Int resw = K_ARRAY::getFromArray3(weight, varStringw, fw, 
                                        niw, njw, nkw, cnw, eltTypew);
  
    if (resw != 1 && resw != 2) 
    {
      RELEASESHAREDB(res, array, f, cn);
      PyErr_SetString(PyExc_TypeError, 
                      "barycenter: invalid weight.");
      return NULL;
    }
    if (fw->getSize() != f->getSize())
    {
      RELEASESHAREDB(res, array, f, cn);
      PyErr_SetString(PyExc_TypeError, 
                      "barycenter: size of weight is different of size of array.");
      return NULL;
    }
    K_COMPGEOM::weightedBarycenter(
      f->getSize(), f->begin(posx), f->begin(posy),
      f->begin(posz), fw->begin(), xb, yb, zb);
    RELEASESHAREDB(resw, weight, fw, cnw);
  }

  RELEASESHAREDB(res, array, f, cn);

  PyObject* tpl;
  tpl = Py_BuildValue("[f,f,f]", xb, yb, zb);
  return tpl;
}

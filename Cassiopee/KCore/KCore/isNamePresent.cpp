/*    
    Copyright 2013-2026 ONERA.

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
#include "kcore.h"

//=============================================================================
// Interface python de isNamePresent
//=============================================================================
PyObject* K_KCORE::isNamePresent(PyObject* self, PyObject* args)
{
  PyObject* array; char* varName;
  if (!PYPARSETUPLE_(args, O_ S_, &array, &varName)) return NULL;

  char* varString;
  E_Int res = K_ARRAY::getVarStringFromArray(array, varString);
  
  E_Int posvar = -1;
  if (res == 1)
  {
    posvar = K_ARRAY::isNamePresent(varName, varString);
    return Py_BuildValue("l", long(posvar));
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "isNamePresent: argument is not an array.");
    return NULL;
  }
}
//=============================================================================
PyObject* K_KCORE::isCoordinateXPresent(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array))
  {
    PyErr_SetString(PyExc_TypeError,
                    "isCoordinateXPresent: wrong arguments.");
    return NULL;
  }

  char* varString;
  E_Int res = K_ARRAY::getVarStringFromArray(array, varString);
  
  E_Int posvar = -1;
  if (res == 1)
  {
    posvar = K_ARRAY::isCoordinateXPresent(varString);
    return Py_BuildValue("l", long(posvar));
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "isCoordinateXPresent: argument is not an array.");
    return NULL;
  }
}
//=============================================================================
PyObject* K_KCORE::isCoordinateYPresent(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array))
  {
    PyErr_SetString(PyExc_TypeError,
                    "isCoordinateYPresent: wrong arguments.");
    return NULL;
  }

  char* varString;
  E_Int res = K_ARRAY::getVarStringFromArray(array, varString);
  
  E_Int posvar = -1;
  if (res == 1)
  {
    posvar = K_ARRAY::isCoordinateYPresent(varString);
    return Py_BuildValue("l", long(posvar));
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "isCoordinateYPresent: argument is not an array.");
    return NULL;
  }
}
//=============================================================================
PyObject* K_KCORE::isCoordinateZPresent(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array))
  {
    PyErr_SetString(PyExc_TypeError,
                    "isCoordinateZPresent: wrong arguments.");
    return NULL;
  }

  char* varString;
  E_Int res = K_ARRAY::getVarStringFromArray(array, varString);
  
  E_Int posvar = -1;
  if (res == 1)
  {
    posvar = K_ARRAY::isCoordinateZPresent(varString);
    return Py_BuildValue("l", long(posvar));
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "isCoordinateZPresent: argument is not an array.");
    return NULL;
  }
}

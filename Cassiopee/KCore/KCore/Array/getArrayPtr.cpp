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

# include "Array/Array.h"

//=============================================================================
E_Float* K_ARRAY::getFieldPtr(PyObject* o)
{
  PyObject* a = PyList_GetItem(o, 1);
  return (E_Float*)PyArray_DATA((PyArrayObject*)a);
}

//=============================================================================
E_Int* K_ARRAY::getConnectPtr(PyObject* o)
{
  PyObject* a = PyList_GetItem(o, 2);
  return (E_Int*)PyArray_DATA((PyArrayObject*)a);
}

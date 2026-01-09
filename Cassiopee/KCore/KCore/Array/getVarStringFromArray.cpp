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

# include "Array/Array.h"

//=============================================================================
// Retourne la chaine des variables de l'array
//==============================================================================
E_Int K_ARRAY::getVarStringFromArray(PyObject* o, char*& varString)
{
  if (PyList_Check(o) == 0)
  {
    PyErr_Warn(PyExc_Warning,
               "getVarStringFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check list.");
    return -1;
  }
  if (PyList_Size(o) != 4 && PyList_Size(o) != 5)
  {
    PyErr_Warn(PyExc_Warning, 
               "getVarStringFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. Check number of elements in list.");
    return -2;
  }

  // -- varString --
  PyObject* l = PyList_GetItem(o,0);
  
  if (PyString_Check(l)) varString = PyString_AsString(l);  
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(l)) varString = (char*)PyUnicode_AsUTF8(l); 
#endif
  else
  {
    PyErr_Warn(PyExc_Warning,
               "getVarStringFromArray: an array must be a list of type ['vars', a, ni, nj, nk] or ['vars', a, c, 'ELTTYPE']. First element must be a string.");
    return -3;
  }
  return 1;
}

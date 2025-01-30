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

#include "cplot.h" 

//=============================================================================
// Extrait des strings d'un pyObj representant une chaine de strings
// retourne une list<char*>
//=============================================================================
int getStringsFromPyObj(PyObject* obj, std::vector<char*>& strings)
{
  if (PyList_Check(obj) == false) return 0; // failed
  E_Int l = PyList_Size(obj);
  for (E_Int i = 0; i < l; i++)
  {
    PyObject* o = PyList_GetItem(obj, i);
    if (PyString_Check(o))
    {
      char* s = new char[128];
      strcpy(s, PyString_AsString(o));
      strings.push_back(s);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(o))
    {
      char* s = new char[128];
      strcpy(s, PyUnicode_AsUTF8(o));
      strings.push_back(s);
    }
#endif
  }
  return 1;
}

//=============================================================================
// Extrait une string d'un pyObj representant string
// Retourne un char* alloue avec new
//=============================================================================
int getStringFromPyObj(PyObject* obj, char*& string)
{
  if (PyString_Check(obj))
  {
    char* s = new char[128];
    strcpy(s, PyString_AsString(obj));
    string = s;
    return 1;
  }
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(obj))
  {
    char* s = new char[128];
    strcpy(s, PyUnicode_AsUTF8(obj));
    string = s;
    return 1;
  }
#endif
  else
  {
    string = NULL;
    return 0; // fail
  } 
}

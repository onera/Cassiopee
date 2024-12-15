/*    
    Copyright 2013-2024 Onera.

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

#include "occ.h"
#include <map>

// ============================================================================
/* Identify int in __tag__ array */
// ============================================================================
PyObject* K_OCC::identifyTags(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, 
                               cn, eltType); 
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "identifyTags: invalid array.");
    return NULL;
  }

  E_Int post = K_ARRAY::isNamePresent(varString, "__tag__");
  if (post < 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "identifyTags: __tag__ not found.");
    return NULL;
  }
  post++;
  E_Float* ft = f->begin(post);

  std::map<E_Int, E_Int> kmap;
  E_Int n = f->getSize();
  E_Int val;
  for (E_Int i = 0; i < n; i++)
  {
    val = round(ft[i]);
    kmap[val] = 1;
  }

  // Export a list
  PyObject* l = PyList_New(0);
  for (const auto& pair : kmap) 
  {
    PyObject* o = PyLong_FromLong(pair.first);
    PyList_Append(l, o);
  }
  RELEASESHAREDB(res, array, f, cn);
  return l;
}
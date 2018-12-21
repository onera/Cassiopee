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
# include <stdlib.h>
# include <string.h>
# include <vector>

# include "converter.h"
# include "kcore.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert array to a node array */
// ============================================================================
PyObject* K_CONVERTER::convertArray2Node(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  PyObject* tpl;
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray(array, varString, 
                              f, nil, njl, nkl, cn, eltType, true);
  
  // Building numpy arrays and return it.
  if (res == 1)      // Structured
  {
    FldArrayI* cnl = new FldArrayI();
    tpl = K_ARRAY::buildArray(*f, varString, *cnl, 0);
    delete cnl;
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2) // Unstructured
  {
    FldArrayI* cnl = new FldArrayI();
    tpl = K_ARRAY::buildArray(*f, varString, *cnl, 0);
    delete cnl;
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "convertArray2Node: array is invalid.");
    return NULL;
  }
}

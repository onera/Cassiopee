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

# include "transform.h"

//=============================================================================
PyObject* K_TRANSFORM::breakElements(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, "O", &array)) return NULL;

  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString,
                               f, ni, nj, nk, cn, eltType);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "breakElements: array is invalid.");
    return NULL;
  }

  PyObject* l = NULL;
  if (K_STRING::cmp(eltType, "NGON") == 0)
    l = breakNGonElements(*f, *cn, varString);
  else if (K_STRING::cmp(eltType, "MIXED") == 0)
    l = breakMixedElements(*f, *cn, varString);
  else  // BE/ME
  {
    // BE: return a copy, ME: break into a list of BEs
    E_Int api = f->getApi();
    E_Int nc = cn->getNConnect();
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);

    l = PyList_New(0);
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
      PyObject* tpl = K_ARRAY::buildArray3(*f, varString, cm, eltTypes[ic], api);
      PyList_Append(l, tpl); Py_DECREF(tpl);
    }

    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  }

  RELEASESHAREDU(array, f, cn);
  return l;
}


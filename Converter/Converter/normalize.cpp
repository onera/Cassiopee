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

// normalize a list of vars

# include "converter.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
// normalize a list of vars
//=============================================================================
PyObject* K_CONVERTER::normalize(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* varsO;
  if (!PyArg_ParseTuple(args, "OO", &array, &varsO )) return NULL;

  // Check vars
  if (PyList_Check(varsO) == 0)
  {
    PyErr_SetString(PyExc_TypeError, "normalize: vars must be a list.");
    return NULL;
  }
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int res = K_ARRAY::getFromArray2(array, varString, 
                                     f, ni, nj, nk, cn, eltType);
  
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "normalize: array is invalid.");
    return NULL;
  }

  // Extraction des variables 
  vector<E_Int> pos;
  E_Int n = 0;
  char* var;
  E_Int m;
  for (E_Int v  = 0 ; v < PyList_Size(varsO); v++)
  {
    PyObject* l = PyList_GetItem(varsO, v);
    if (PyString_Check(l))
    {
      var = PyString_AsString(l);
      m = K_ARRAY::isNamePresent(var, varString);
      if (m == -1)
        printf("Warning: normalize: variable %d not present in array. Skipped...\n", v);
      else {m++; pos.push_back(m);}
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l))
    {
      var = PyBytes_AsString(PyUnicode_AsUTF8String(l));
      m = K_ARRAY::isNamePresent(var, varString);
      if (m == -1)
        printf("Warning: normalize: variable %d not present in array. Skipped...\n", v);
      else {m++; pos.push_back(m);}
    }
#endif
    else
    {
      printf("Warning: normalize: invalid string for variable %d. Skipped...\n", v);
    }  
  }

  n = pos.size();
  if (n == 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    printf("Warning: normalize: no variable in result.\n");;
    Py_INCREF(Py_None); return Py_None;
  }
  E_Int npts = f->getSize();

  FldArrayF fp(npts);
  E_Float* fpt0 = fp.begin();

#pragma omp parallel default(shared)
  {
#pragma omp for nowait
    for (E_Int i = 0; i < npts; i++) fpt0[i] = 0.;

    for (E_Int v = 0; v < n; v++)
    {
      E_Float* ft = f->begin(pos[v]);
#pragma omp for nowait
      for (E_Int i = 0; i < npts; i++) fpt0[i] += ft[i]*ft[i];
    }
#pragma omp for nowait
    for (E_Int i = 0; i < npts; i++)
      fpt0[i] = 1./K_FUNC::E_max(sqrt(fpt0[i]), 1.e-12);

    for (E_Int v  = n-1; v >= 0; v--)
    {
      E_Float* ft = f->begin(pos[v]);
#pragma omp for
      for (E_Int i = 0; i < npts; i++)
      {
        ft[i] = ft[i] * fpt0[i];
      }
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
} 
    

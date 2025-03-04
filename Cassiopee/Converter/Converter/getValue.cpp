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

// getValue
// Cette fonction est obsolete et a ete avantageusement remplacee par du 
// code en python

# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/* Get the values in an array for point of index ind or (i,j,k) */
// ============================================================================
PyObject* K_CONVERTER::getValueOfArray(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int i0, j0, k0;
  E_Int ind0;
  E_Int i = -1; 
  E_Int j = -1; 
  E_Int k = -1; 
  E_Int ind = -1;

  if (PYPARSETUPLE_(args, O_ TIII_, &array, &i0, &j0, &k0))
  { i = i0; j = j0; k = k0; }
  else if (PYPARSETUPLE(args, O_ I_, &array, &ind0)) ind = ind0;
  else return NULL;

  // Check array
  E_Int im, jm, km;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;  
  E_Int res = K_ARRAY::getFromArray(array, varString, f, 
                                    im, jm, km, cn, eltType, true);
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getValue: array is invalid."); 
    return NULL;
  }
  if (nfld <= 0 || npts <= 0) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "getValue: no field defined in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  E_Float* vals = new E_Float[nfld];
  if (res == 1)
  { 
    if (i != -1 && j != -1 && k != -1)
    {
      ind = (i-1) + (j-1) * im + (k-1)* im*jm; 
      if (ind < 0 || ind >= npts)
      {
        char error[256];
        sprintf(error, "getValue: index (%d) is out of bounds (0-%d).", 
                ind, npts-1);
        PyErr_SetString(PyExc_ValueError, error);
        RELEASESHAREDS(array, f); return NULL;
      }
    }
    else 
    {
      if (ind < 0 || ind >= npts)
      {
        char error[256];
        sprintf(error, "getValue: index (%d) is out of bounds (0-%d).", 
                ind, npts-1);
        PyErr_SetString(PyExc_ValueError, error);
        RELEASESHAREDS(array, f); return NULL;
      }
    }
    
    for (E_Int eq = 1; eq <= nfld; eq++) vals[eq-1] = (*f)(ind, eq);
    RELEASESHAREDS(array, f);
  }
  else if (res == 2)
  {
    if (ind == -1)
    {
      PyErr_SetString(PyExc_ValueError,
                      "getValue: a single index must be used for unstructured arrays.");
      RELEASESHAREDU(array, f, cn); return NULL;
    }  

    if (ind < 0 || ind >= f->getSize())
    {
      char error[256];
      sprintf(error, "getValue: index (%d) is out of bounds (0-%d).", 
              ind, f->getSize()-1);
      PyErr_SetString(PyExc_ValueError, error);
      RELEASESHAREDU(array, f, cn);  return NULL;
    }
    for (E_Int eq = 1; eq <= nfld; eq++) vals[eq-1] = (*f)(ind, eq);
    RELEASESHAREDU(array, f, cn); 
  }
  
  // build list of values 
  PyObject* tpl;
  PyObject* l = PyList_New(0);
#ifdef E_DOUBLEREAL
  for (E_Int i = 0; i < nfld; i++)
  {
    tpl = Py_BuildValue("d", vals[i]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
#else
  for (E_Int i = 0; i < nfld; i++)
  {
    tpl = Py_BuildValue("f", vals[i]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
#endif
  delete [] vals;
  return l;
}

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
#include "Array/Array.h"
#include "String/kstring.h"
#include <stdio.h>
#include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Add a field in an array1 or array2 (in place) 
   IN: array: input array
   IN: varName: variable name to add
   OUT: array (modified) */
//=============================================================================
void K_ARRAY::addFieldInArray(PyObject* array, char* varName)
{
  IMPORTNUMPY;

  // Le champ existe deja?
  PyObject* o = PyList_GetItem(array, 0);
  char* varString = NULL; 
  if (PyString_Check(o)) varString = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(o)) varString = PyBytes_AsString(PyUnicode_AsUTF8String(o)); 
#endif
  E_Int pos = K_ARRAY::isNamePresent(varString, varName);
  if (pos != -1) return;

  // Ajout du champ dans varString
  E_Int l = strlen(varString)+strlen(varName)+2;
  char newVarstring[l];
  strcpy(newVarstring, varString);
  strcat(newVarstring, ",");
  strcat(newVarstring, varName);
#if PY_VERSION_HEX >= 0x03000000
  PyList_SetItem(array, 0, PyUnicode_FromString(newVarstring));
#else
  PyList_SetItem(array, 0, PyString_FromString(newVarstring));
#endif
  
  // Recuperation des champs
  PyObject* a = PyList_GetItem(array, 1);

  // Recuperation api
  E_Int api = 1;
  if (PyList_Check(a) == true) api = 2;

  if (api == 1) // Array1
  { 
    PyArrayObject* ar = (PyArrayObject*)a;
    E_Int s = PyArray_DIMS(ar)[1];
    E_Int nfld = PyArray_DIMS(ar)[0];
    E_Float* pto = (E_Float*)PyArray_DATA(ar);
    npy_intp dim[2];  
    dim[1] = s;
    dim[0] = nfld+1;
    PyArrayObject* an = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    E_Float* pt = (E_Float*)PyArray_DATA(an);
#pragma omp parallel
    {
    E_Float* ptl; E_Float* ptlo;
    for (E_Int n = 0; n < nfld; n++) // recopie
    {
      ptl = pt+n*s; ptlo = pto+n*s;
#pragma omp for
      for (E_Int i = 0; i < s; i++) ptl[i] = ptlo[i];   
    }
    }
    PyList_SetItem(array, 1, (PyObject*)an);
  }
  else // Array2
  {
    PyArrayObject* ar = (PyArrayObject*)PyList_GetItem(a, 0);
    E_Int s = PyArray_SIZE(ar);
    npy_intp dim[1];
    dim[0] = s;
    PyArrayObject* an = (PyArrayObject*)PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    PyList_Append(a, (PyObject*)an);
    Py_DECREF(an);  
  }
  return;
}

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

// get the magnitude (norm) of 3 variables

# include "converter.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
// normalize a list of vars
//=============================================================================
PyObject* K_CONVERTER::magnitude(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* varsO;
  if (!PYPARSETUPLE_(args, OO_, &array, &varsO)) return NULL;

  FldArrayF* f; FldArrayI* cn;
  char* eltType; char* varString;
  E_Int res;
  E_Int ni, nj, nk; // number of points of array
  res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);
  
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "magnitude: array is invalid.");
    return NULL;
  }

  // Check vars
  if (PyList_Check(varsO) == 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, "magnitude: vars must be a list.");
    return NULL;
  }

  // Extraction des variables 
  vector<E_Int> pos;
  PyObject* tpl;
  E_Int n = 0;
  char* var;
  E_Int m;
  E_Int varStringL = strlen(varString);
  char* varStringOut = new char [varStringL+10];
  varStringOut[0] = '\0';
  char* commonVariable = new char [varStringL];
  commonVariable[0] = '\0';
  char* start = new char [varStringL];
  E_Bool common = false;

  for (E_Int v  = 0 ; v < PyList_Size(varsO); v++)
  {
    PyObject* l = PyList_GetItem(varsO, v);
    if (PyString_Check(l))
    {
      var = PyString_AsString(l);
      m = K_ARRAY::isNamePresent(var, varString);
      if (m == -1)
        printf("Warning: magnitude: variable " SF_D_ " not present in array. Skipped...\n", v);
      else 
      {
        m++; pos.push_back(m); strcat(varStringOut, var);
        char v = var[strlen(var)-1];
        strcpy(start, var); start[strlen(var)-1] = '\0';
        if (v == 'x' || v == 'y' || v == 'z' || v == 'X' || v == 'Y' || v == 'Z')
        { 
          if (strlen(commonVariable) == 0) 
          { strcpy(commonVariable, start); common = true; }
          else if (strcmp(commonVariable, start) != 0) common = false;
        }
      }
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l))
    {
      var = (char*)PyUnicode_AsUTF8(l);
      m = K_ARRAY::isNamePresent(var, varString);
      if (m == -1)
        printf("Warning: magnitude: variable " SF_D_ " not present in array. Skipped...\n", v);
      else 
      {
        m++; pos.push_back(m); strcat(varStringOut, var);
        char v = var[strlen(var)-1];
        strcpy(start, var); start[strlen(var)-1] = '\0';
        if (v == 'x' || v == 'y' || v == 'z' || v == 'X' || v == 'Y' || v == 'Z')
        { 
          if (strlen(commonVariable) == 0) 
          { strcpy(commonVariable, start); common = true; }
          else if (strcmp(commonVariable, start) != 0) common = false;
        }
      }
    }
#endif
    else
    {
      printf("Warning: magnitude: invalid string for variable " SF_D_ ". Skipped...\n", v);
    }
      
  }

  if (common == true) strcpy(varStringOut, commonVariable);
  strcat(varStringOut, "Magnitude");

  n = pos.size();
  if (n == 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    printf("Warning: magnitude: no variable in result.\n");
    Py_INCREF(Py_None); return Py_None;
  }

  E_Int npts = f->getSize();
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray3(1, varStringOut, ni, nj, nk);
  } 
  else //unstructured 
  {
    E_Int api = f->getApi();
    tpl = K_ARRAY::buildArray3(1, varStringOut, npts, *cn,
                              eltType, -1, api, true); // copy connect
  }
  FldArrayF* norm;
  K_ARRAY::getFromArray3(tpl, norm);
  E_Float* normp = norm->begin();

  FldArrayF& fp = *f;

#pragma omp parallel
  {
    E_Float norm;
    #pragma omp for
    for (E_Int i = 0; i < npts; i++)
    {
      norm = 0.;
      for (E_Int v  = 0; v < n; v++)
      {
        norm += fp(i, pos[v]) * fp(i, pos[v]);
      }
      normp[i] = sqrt(norm);
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, norm);
  delete [] varStringOut;
  delete [] commonVariable;
  delete [] start;
  return tpl;
} 

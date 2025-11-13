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

// getArgMin, getMax, getMean values of an array
// getNegativeValues, getSupValues of an array
# include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
/* Get array value where one of its variables is minimum */
//=============================================================================
PyObject* K_CONVERTER::getArgMin(PyObject* self, PyObject* args)
{
  PyObject* array;
  char* varName;
  if (!PYPARSETUPLE_(args, O_ S_, &array, &varName)) return NULL;

  // Check array
  E_Int im, jm, km;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);
  E_Int posvar = -1;
  if (res == 1 || res == 2)
  {
    posvar = K_ARRAY::isNamePresent(varName, varString);
    if (posvar == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getArgMin: variable %s is not present in array.", 
              varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }
    posvar++;

    E_Int fsize = f->getSize();
    if (fsize == 0)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getArgMin: array for variable %s is empty.", varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }
   
    E_Float* fp = f->begin(posvar);
    E_Int nfld = f->getNfld();
    E_Int ind = 0;
    E_Float minVal = K_CONST::E_MAX_FLOAT;

#pragma omp parallel default(shared)
    {
      E_Float minValLoc = K_CONST::E_MAX_FLOAT; E_Int indLoc = 0;
#pragma omp for nowait
      for (E_Int i = 0; i < fsize; i++)
      {
        if (fp[i] < minValLoc) { minValLoc = fp[i]; indLoc = i; }
      }
#pragma omp critical
      {
        if (minValLoc < minVal) { minVal = minValLoc; ind = indLoc; }
      }
    }

    // Build return
    FldArrayF field(1, nfld);
    field.setAllValuesAtNull();

    for (E_Int n = 1; n <= nfld; n++) field(0, n) = (*f)(ind, n);

    RELEASESHAREDB(res, array, f, cn);

    PyObject* l = PyList_New(0);

    for (E_Int i = 0; i < nfld; i++)
    {
      PyObject* tpl = Py_BuildValue(R_, field(0, i+1));
      PyList_Append(l, tpl); Py_DECREF(tpl);
    }
    return l;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getArgMin: argument is not an array.");
    return NULL;
  }
}

//=============================================================================
/* Get array value where one of its variables is maximum */
//=============================================================================
PyObject* K_CONVERTER::getArgMax(PyObject* self, PyObject* args)
{
  PyObject* array;
  char* varName;  
  if (!PYPARSETUPLE_(args, O_ S_, &array, &varName)) return NULL;

  // Check array
  E_Int im, jm, km;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);

  E_Int posvar = -1;
  if (res == 1 || res == 2)
  {
    posvar = K_ARRAY::isNamePresent(varName, varString);
    if (posvar == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getArgMax: variable %s is not present in array.", 
              varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }
    posvar++;

    E_Int fsize = f->getSize();
    if (fsize == 0)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getArgMax: array for variable %s is empty.", varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }
   
    E_Float* fp = f->begin(posvar);
    E_Int nfld = f->getNfld();
    E_Int ind = 0;
    E_Float maxVal = -K_CONST::E_MAX_FLOAT;

#pragma omp parallel default(shared)
    {
      E_Float maxValLoc = -K_CONST::E_MAX_FLOAT; E_Int indLoc = 0;
#pragma omp for nowait
      for (E_Int i = 0; i < fsize; i++)
      {
        if (fp[i] > maxValLoc) { maxValLoc = fp[i]; indLoc = i; }
      }
#pragma omp critical
      {
        if (maxValLoc > maxVal) { maxVal = maxValLoc; ind = indLoc; }
      }
    }

    // Build return
    FldArrayF field(1, nfld);
    field.setAllValuesAtNull();

    for (E_Int n = 1; n <= nfld; n++) field(0, n) = (*f)(ind, n);

    RELEASESHAREDB(res, array, f, cn);

    PyObject* l = PyList_New(0);

    for (E_Int i = 0; i < nfld; i++)
    {
      PyObject* tpl = Py_BuildValue(R_, field(0, i+1));
      PyList_Append(l, tpl); Py_DECREF(tpl);
    }
    return l;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getArgMax: argument is not an array.");
    return NULL;
  }
}

//=============================================================================
/* Get mean value of an array */
//=============================================================================
PyObject* K_CONVERTER::getMeanValue(PyObject* self, PyObject* args)
{
  PyObject* array; char* varName;
  if (!PYPARSETUPLE_(args, O_ S_, &array, &varName)) return NULL;

  // Check array
  E_Int im, jm, km;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);

  E_Int posvar = -1;
  if (res == 1 || res == 2)
  {
    posvar = K_ARRAY::isNamePresent(varName, varString);
    if (posvar == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getMeanValue: variable %s is not present in array.", 
              varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }
    posvar++;

    E_Int fsize = f->getSize();
    if (fsize == 0)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getArgMean: array for variable %s is empty.", varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }

    E_Float* fp = f->begin(posvar);
    E_Float meanVal = 0.;

#pragma omp parallel for default(shared) reduction(+:meanVal) 
    for (E_Int i = 0; i < fsize; i++) meanVal += fp[i];

    meanVal = meanVal / E_Float(fsize);
    RELEASESHAREDB(res, array, f, cn);
    return Py_BuildValue(R_, meanVal);
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getMeanValue: argument is not an array.");
    return NULL;
  }
}

//=============================================================================
/* Get mean value of a range in an array */
//=============================================================================
PyObject* K_CONVERTER::getMeanRangeValue(PyObject* self, PyObject* args)
{
  PyObject* array; char* varName;
  E_Float rmin, rmax;
  if (!PYPARSETUPLE_(args, O_ S_ RR_, &array, &varName, &rmin, &rmax)) return NULL;

  // Check array
  E_Int im, jm, km;
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);

  E_Int posvar = -1;
  if (res == 1 || res == 2)
  {
    posvar = K_ARRAY::isNamePresent(varName, varString);
    if (posvar == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getMeanRangeValue: variable %s is not present in array.", 
              varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }
    posvar++;

    E_Int fsize = f->getSize();
    if (fsize == 0)
    {
      RELEASESHAREDB(res, array, f, cn);
      char error[256];
      sprintf(error, "getMeanRangeValue: array for variable %s is empty.", varName);
      PyErr_SetString(PyExc_TypeError, error);
      return NULL;
    }

    E_Float* fp = f->begin(posvar);
    E_Float meanVal = 0.;
    
    std::vector<E_Float> vals;
    vals.insert(vals.end(), fp, fp+fsize);
    
    std::sort(vals.begin(), vals.end());
    
    size_t sz1 = size_t(E_Float(fsize)*rmin);
    size_t sz2 = size_t(E_Float(fsize)*rmax);
    
    for (size_t i = sz1; i < sz2; ++i) meanVal += vals[i];

    meanVal = meanVal / E_Float(1+sz2-sz1);

    RELEASESHAREDB(res, array, f, cn);
    return Py_BuildValue(R_, meanVal);
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "getMeanRangeValue: argument is not an array.");
    return NULL;
  }
}


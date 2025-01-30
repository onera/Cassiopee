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

// Computes the L0 and L2 norms of an array

# include "converter.h"
# include "kcore.h"
using namespace K_FUNC;
using namespace K_FLD;

//=============================================================================
// L0 norm of an array
//=============================================================================
PyObject* K_CONVERTER::normL0(PyObject* self, PyObject* args)
{
  PyObject* array;
  char* varName;
  if (!PyArg_ParseTuple(args, "Os", &array, &varName)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; char* varString; 
  char* eltType; FldArrayI* cn;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);
  
  E_Int npts = f->getSize();
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "normL0: invalid array." );
    return NULL;
  }

  E_Float L0err = 0.; 
  E_Int posv = K_ARRAY::isNamePresent(varName, varString); posv++;
  E_Int posc = K_ARRAY::isCellNatureField2Present(varString); posc++;
  if (posv == 0)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "normL0: variable not found in array.");
    return NULL;
  }

  E_Float* fv = f->begin(posv);

  E_Int nthreads = __NUMTHREADS__;
  E_Float delta = npts/nthreads;
  E_Float* fmax = new E_Float [nthreads];
  for (E_Int i = 0; i < nthreads; i++) fmax[i] = -1.;
  
  if (posc != 0) // celln exists
  {
    E_Float* fc = f->begin(posc);
#pragma omp parallel
    {
      E_Float celln;
      E_Int ithread = __CURRENT_THREAD__;
      E_Int istart = ithread*delta;
      E_Int iend = (ithread+1)*delta;
      if (ithread == nthreads-1) iend = npts;
      for (E_Int i = istart; i < iend; i++)
      {
        celln = (fc[i] == 0.) ? 0. : 1.; 
        fmax[ithread] = E_max(fmax[ithread], fv[i]) * celln;
      }
    }
  }
  else
  {
#pragma omp parallel
    {
      E_Int ithread = __CURRENT_THREAD__;
      E_Int istart = ithread*delta;
      E_Int iend = (ithread+1)*delta;
      if (ithread == nthreads-1) iend = npts;
      for (E_Int i = istart; i < iend; i++)
      {
        fmax[ithread] = E_max(fmax[ithread], fv[i]);
      }
    }
  }

  for (E_Int i = 0; i < nthreads; i++) L0err = E_max(L0err, fmax[i]);
  delete [] fmax;

  RELEASESHAREDB(res, array, f, cn);
  
#ifdef E_DOUBLEREAL
  return Py_BuildValue("d", L0err);
#else
  return Py_BuildValue("f", L0err);
#endif    
}

//=============================================================================
// L2 norm of an array
//=============================================================================
PyObject* K_CONVERTER::normL2(PyObject* self, PyObject* args)
{
  PyObject* array;
  char* varName;
  if (!PyArg_ParseTuple(args, "Os", &array, &varName)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);

  E_Int npts = f->getSize();

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "normL2: invalid array." );
    return NULL;
  }
  E_Float L2err = 0.; 
  E_Int posv = K_ARRAY::isNamePresent(varName, varString); posv++;
  E_Int posc = K_ARRAY::isCellNatureField2Present(varString); posc++;
  if (posv == 0)
  {
    PyErr_SetString( PyExc_TypeError, 
                     "normL2: variable not found in array.") ;
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  E_Float* fv = f->begin(posv);
  if (posc != 0) //celln exists
  {
    E_Float* fc = f->begin(posc);
#pragma omp parallel for reduction(+:L2err)
    for (E_Int ind = 0; ind < npts; ind++)
    {
      E_Float celln = (fc[ind] == 0.) ? 0. : 1.; 
      E_Float val = fv[ind] * celln;
      L2err += val*val;
    }
  }
  else
  {
#pragma omp parallel for reduction(+:L2err)
    for (E_Int ind = 0; ind < npts; ind++)
    {
      E_Float val = fv[ind];
      L2err += val*val;
    }
  }

  if (npts != 0) L2err = sqrt(L2err / npts);
  
  RELEASESHAREDB(res, array, f, cn);
  
#ifdef E_DOUBLEREAL
  return Py_BuildValue("d", L2err);
#else
  return Py_BuildValue("f", L2err);
#endif
}

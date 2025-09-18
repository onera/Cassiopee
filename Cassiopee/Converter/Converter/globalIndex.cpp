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
# include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* create global index as a field
   IN: array in vertex with 'globalIndex' field
   IN: start: starting numbering
 */
//=============================================================================
PyObject* K_CONVERTER::createGlobalIndex(PyObject* self, PyObject* args)
{
  PyObject *array;
  E_Int start;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &start)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "createGlobalIndex: invalid array.");
    return NULL;
  }

  E_Int pos = K_ARRAY::isNamePresent("globalIndex", varString);
  if (pos == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "createGlobalIndex: no globalIndex field found in array.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }
  pos += 1; 
  E_Float* fp = f->begin(pos);
  E_LONG* gi = (E_LONG*)fp;

#pragma omp parallel
  {
#pragma omp for
    for (E_Int ind = 0; ind < f->getSize(); ind++)
    {
      gi[ind] = (E_LONG)ind;
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* recover global index
   IN: array in vertex with 'globalIndex' field
   OUT: arrayo: output array with recovered field
 */
//=============================================================================
PyObject* K_CONVERTER::recoverGlobalIndex(PyObject* self, PyObject* args)
{
  PyObject *array;
  PyObject *arrayo;
  if (!PYPARSETUPLE_(args, OO_, &array, &arrayo)) return NULL;

  // Check input array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "recoverGlobalIndex: invalid array.");
    return NULL;
  }

  E_Int pos = K_ARRAY::isNamePresent("globalIndex", varString);
  if (pos == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "recoverGlobalIndex: no globalIndex field found in array.");
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }
  pos += 1; 
  E_Float* find = f->begin(pos);
  E_LONG* indp = (E_LONG*)find;

  // Get var list of input array
  vector<char*> vars;
  K_ARRAY::extractVars(varString, vars);

  // Check output array
  E_Int nio, njo, nko;
  FldArrayF* fo; FldArrayI* cno;
  char* varStringo; char* eltTypeo;
  E_Int reso = K_ARRAY::getFromArray3(arrayo, varStringo, fo, nio, njo, nko, 
                                      cno, eltTypeo);
  if (reso != 1 && reso != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "recoverGlobalIndex: invalid output array.");
    return NULL;
  }
  E_Int sizeo = fo->getSize();
  
  // Recovering
#pragma omp parallel
  {
    E_Int ind;
    for (E_Int n = 1; n <= f->getNfld(); n++)
    {
      if (n != pos)
      {
        char* v = vars[n-1];
        E_Int pv = K_ARRAY::isNamePresent(v, varStringo);
        if (pv != -1)
        {
          E_Float* fp = f->begin(n);
          E_Float* fop = fo->begin(pv+1);
#pragma omp for
          for (E_Int i = 0; i < f->getSize(); i++)
          {
            ind = indp[i];
            if (ind < sizeo) fop[ind] = fp[i]; 
          }
        }
      }
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
}

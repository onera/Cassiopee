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

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "converter.h"
# include "kcore.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
// extract vars from an array
//=============================================================================
PyObject* K_CONVERTER::extractVars(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* vars; PyObject* tpl;
  if (!PYPARSETUPLE_(args, OO_, &array, &vars)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res, nv, n;
  res = K_ARRAY::getFromArray3(array, varString, f,
                               nil, njl, nkl, cn, eltType);

  if (res == 1) nv = f->getNfld();
  else if (res == 2) nv = f->getNfld();
  else return NULL; // errors are alread set

  // Check vars
  if (PyList_Check(vars) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "extractVars: vars must be a list."); 
    RELEASESHAREDB(res, array, f, cn);
    return NULL;
  }

  vector<E_Int> nvar;
  for (E_Int i = 0; i < PyList_Size(vars); i++)
  {
    tpl = PyList_GetItem(vars, i);
    if (PyInt_Check(tpl) == 0)
    {
      printf("Warning: extractVars: Vars must be integers. Skipped...\n");
    }
    else
    {
      n = PyLong_AsLong(tpl);
      if (n <= nv && n >=1) nvar.push_back(n);
      else
        printf("Warning: extractVars: Variable is not in array. Skipped...\n");
    }
  }

  // Sortie
  if (nvar.size() == 0)
  {
    printf("Warning: extractVars: no variable in result.\n");
    RELEASESHAREDB(res, array, f, cn);
    Py_INCREF(Py_None);
    return Py_None;
  }

  E_Int nt = nvar.size();
  E_Int fSize = f->getSize();

  // Building varString
  vector<char*> vvars;
  K_ARRAY::extractVars(varString, vvars);
  char* fstring = new char [strlen(varString)+1];
  strcpy(fstring, vvars[nvar[0]-1]);
  if (nt > 1) strcat(fstring, ",");
  for (E_Int i = 1; i < nt; i++)
  {
    strcat(fstring, vvars[nvar[i]-1]);
    if (i != nt-1) strcat(fstring, ",");
  }
  E_Int sizevvars = vvars.size();
  for (E_Int i = 0; i < sizevvars; i++) delete [] vvars[i];

  // Build array here
  E_Int api = f->getApi();
  if (res == 1)
  {
    tpl = K_ARRAY::buildArray3(nt, fstring, nil, njl, nkl, api);
  }
  else
  {
    E_Int api = f->getApi();
    E_Bool compact = false;
    if (api == 1) compact = true;
    FldArrayF f2(fSize, nt, compact);
    tpl = K_ARRAY::buildArray3(f2, fstring, *cn, eltType, api);
  }
  FldArrayF* f2; K_ARRAY::getFromArray3(tpl, f2);

#pragma omp parallel
  {
    for (E_Int i = 0; i < nt; i++)
    {
      E_Float* f2p = f2->begin(i+1);
      E_Float* fp = f->begin(nvar[i]);
#pragma omp for
      for (E_Int j = 0; j < fSize; j++) f2p[j] = fp[j];
    }
  }
  delete [] fstring;
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

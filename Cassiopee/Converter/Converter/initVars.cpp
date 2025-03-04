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
/* Init fields given in nameArray to the constant value val */
//=============================================================================
PyObject* K_CONVERTER::initVars(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float val; char* varName;

  if (!PYPARSETUPLE_(args, O_ S_ R_, &array, &varName, &val)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, 
                                     ni, nj, nk, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "initVars: invalid array definition.");
    return NULL;
  }

  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  E_Int npts = f->getSize();
  if (res == 1) //structured
  {
    E_Int nfld = f->getNfld();
    tpl = K_ARRAY::buildArray3(nfld, varString, ni, nj, nk);
  } 
  else //unstructured 
  {
    tpl = K_ARRAY::buildArray3(*f, varString, *cn, eltType);
  }
  FldArrayF* f2; K_ARRAY::getFromArray3(tpl, f2);

  E_Int posvar = K_ARRAY::isNamePresent(varName, varString)+1;
  if (posvar == 0)
  {
    printf("Warning: initVars: variable name %s is not in array. Skipped...\n",varName);
  }
  else 
  {
    E_Float* f2p = f2->begin(posvar);
#pragma omp parallel for if (npts > __MIN_SIZE_MEAN__)
    for (E_Int i = 0; i < npts; i++) f2p[i] = val;
  }

  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

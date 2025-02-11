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

// randomize a variable

# include "converter.h"

using namespace K_FLD;
using namespace std;

// =============================================================================
// Perturbe aleatoirement un champ autour de sa valeur a dans l'intervalle 
// [a*(1-alpha),a*(1+alpha)]
//==============================================================================
PyObject* K_CONVERTER::randomizeVar(PyObject* self, PyObject* args)
{
  PyObject* array; char* varname;
  E_Float deltamin, deltamax;
  if (!PYPARSETUPLE_(args, O_ S_ RR_, &array, &varname, &deltamin, &deltamax)) return NULL;

  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; 
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, eltType, true);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "randomizeVar: array is invalid.");
    return NULL;
  }
  // Construit l'array resultat et l'initialise par copie
  PyObject* tpl;
  E_Int npts = f->getSize();
  E_Int nfld = f->getNfld();
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nfld, varString, 
                              ni, nj, nk);
  } 
  else //unstructured 
  {
    E_Int csize = cn->getSize()*cn->getNfld(); 
    tpl = K_ARRAY::buildArray(nfld, varString,
                              npts, cn->getSize(),
                              -1, eltType, false, csize);
  }
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true); fn.setAllValuesAt(*f);

  if (res == 2)
  {
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    E_Int* cnp = cn->begin();
    E_Int size = cn->getSize()*cn->getNfld();
    E_Int i;
#pragma omp parallel for shared(size,cnnp,cnp) private(i)
    for (i = 0; i < size; i++) cnnp[i] = cnp[i];
  }

  E_Int posvar = K_ARRAY::isNamePresent(varname, varString)+1;
  if (posvar == 0)
  {
    printf("Warning: randomizeVar: variable %s not found. Skipped...\n",varname);
  }
  else 
  {
    E_Float* fnp = fn.begin(posvar);
    E_Int i;
    E_LONG idum = -1;
    E_Float rfactor;
    for (i = 0; i < npts; i++) 
    {
      rfactor = (deltamax+deltamin)*K_NOISE::stdRand(&idum) - deltamin;
      fnp[i] += rfactor;
    }
  }

  RELEASESHAREDB(res, array, f, cn);
  return tpl;
} 
    

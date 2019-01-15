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
#include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert HO mesh to LO mesh */
// ============================================================================
PyObject* K_CONVERTER::convertHO2LO(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int mode;
  if (!PYPARSETUPLEI(args, "Ol", "Od", &array, &mode)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray2(array, varString, 
                               f, ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2)
  {
     PyErr_SetString(PyExc_TypeError, 
                     "convertHO2LO: array is invalid.");
     return NULL;
  }
  if (res == 1)
  {   
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertHO2LO: array must be unstructured.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, 4, "NGON") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "convertHO2LO: array must not be NGON.");
    return NULL;
  }
  
  // Caracteristiques de l'array input
  E_Int nelts = cn->getSize();
  E_Int nfld = f->getNfld();
  E_Int nvertex = f->getSize();
  E_Int api = f->getApi();

  // Caracteristique de l'array de sortie
  char outEltType[128]; E_Int d;
  K_ARRAY::eltString2TypeId(eltType, outEltType, d, d, d);

  // directement buildArray2
  PyObject* o = K_ARRAY::buildArray2(nfld, varString, nvertex, nelts, -1, outEltType, false, 0, 0, 0, api);

  FldArrayF* fo; FldArrayI* co;
  K_ARRAY::getFromArray2(o, fo, co);

  // Ne pas utiliser (*fo) = (*f); peut reallouer 
  fo->copy(*f, 1, nfld); // copie les champs est ok

  K_CONNECT::connectHO2LO(eltType, *cn, *co, 0);
  
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDU(o, fo, co);
  
  return o;
}

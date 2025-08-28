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

# include "../converter.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
// Extract a sub field from array2 from a list of index, return array2
// IN: array [array2]
// IN: index [numpy]
// OUT: array2
//=============================================================================
PyObject* K_CONVERTER::extractFields(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* index;
  if (!PYPARSETUPLE_(args, OO_,
                     &array, &index)) return NULL;
  
  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, nil, njl, nkl, 
                                     cn, eltType);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "extractFields: invalid array.");
    return NULL;
  }

  // Check index
  FldArrayI* inds;
  E_Int res2 = K_NUMPY::getFromNumpyArray(index, inds);

  if (res2 == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "extractFields: index numpy is invalid.");
    return NULL;
  }

  E_Int api = f->getApi();
  E_Int nfld = f->getNfld();
  E_Int npts = inds->getSize();
  E_Int* indp = inds->begin();
  
  // Cree un array, avec les memes champs mais structure a plat
  PyObject* o = K_ARRAY::buildArray2(nfld, varString, npts,1,1, api);
  FldArrayF* fo; K_ARRAY::getFromArray2(o, fo);
  
 #pragma omp parallel default(shared)
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* fpo = fo->begin(n);
#pragma omp for 
      for (E_Int i = 0; i < npts; i++)
      {
        fpo[i] = fp[indp[i]];
      }
    }
  }

  RELEASESHAREDS(o, fo);
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDN(index, inds);
  return o;
}

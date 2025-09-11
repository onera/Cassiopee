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

// exteriorEltsStructured

# include <stdio.h>
# include <string.h>

# include "post.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Cree une seule grille structuree. */
// ============================================================================
PyObject* K_POST::exteriorEltsStructured(PyObject* self, PyObject* args)
{
  /* Check */ 
  PyObject* array; E_Int depth;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &depth)) return NULL;
  
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorEltsStructured: invalid array.");
    return NULL;
  }
  if (res != 1)
  {
    if (res == 2) { RELEASESHAREDU(array, f, cn); }
    PyErr_SetString(PyExc_TypeError,
                    "exteriorEltsStructured: array must be structured.");
    return NULL;
  }

  /* Exterior grid */
  E_Int nfld = f->getNfld();
  
  E_Int api = f->getApi();
  E_Int nis = std::min(2*depth+2, ni);
  E_Int njs = std::min(2*depth+2, nj);
  E_Int nks = std::min(2*depth+2, nk);
  PyObject* o = K_ARRAY::buildArray3(nfld, varString, nis, njs, nks, api);
  FldArrayF* fo; K_ARRAY::getFromArray3(o, fo);
  
  E_Int nij = ni*nj;
  E_Int nijs = nis*njs;
  E_Int ip, jp, kp, ind, indp;
  for (E_Int k = 0; k < nks; k++)
  for (E_Int j = 0; j < njs; j++)
  for (E_Int i = 0; i < nis; i++)
  {
    ind = i+nis*j+nijs*k;
    ip = i;
    if (i >= nis-depth-1) ip = ni-nis+i-1;
    ip = std::max(ip,E_Int(0));
    ip = std::min(ip,ni-1);
    jp = j;
    if (j >= njs-depth-1) jp = nj-njs+j-1;
    jp = std::max(jp,E_Int(0));
    jp = std::min(jp,nj-1);
    kp = k;
    if (k >= nks-depth-1) kp = nk-nks+k-1;
    kp = std::max(kp,E_Int(0));
    kp = std::min(kp,nk-1);
    indp = ip+ni*jp+nij*kp;
    //printf("%d %d %d\n", ip,jp,kp);
    for (E_Int n = 1; n <= nfld; n++) (*fo)(ind,n) = (*f)(indp,n);
  }
  RELEASESHAREDS(array, f);
  RELEASESHAREDS(o, fo);
  return o;
}

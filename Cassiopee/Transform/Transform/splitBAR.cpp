/*
    Copyright 2013-2026 ONERA.

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

# include "transform.h"
# include "Nuga/include/BARSplitter.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/*
   splitBAR:
   Split a bar into 2 bits given 2 nodes of the bar.
*/
//=============================================================================
PyObject* K_TRANSFORM::splitBAR(PyObject* self, PyObject* args)
{
  PyObject* array; E_Int N; E_Int N2;
  if (!PYPARSETUPLE_(args, O_ II_, &array, &N, &N2))
  {
    return NULL;
  }
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitBAR: can not be used on a structured array.");
    return NULL;
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitBAR: invalid array.");
    return NULL;
  }

  if (strcmp(eltType, "BAR") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitBAR: must be used on a BAR-array.");
    return NULL;
  }

  // Duplique le pt et c'est tout. Le reste sera fait par un splitConnexity.
  E_Int npts = f->getSize();
  E_Int N1 = N;

  if (N1 < 0 || N1 > npts-1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "splitBAR: N is incorrect.");
    return NULL;
  }

  E_Int nptsf = npts+1;
  if (N2 > 0) nptsf = npts+2;
  E_Int nfld = f->getNfld();
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString,
                                       nptsf, cn->getSize(),
                                       eltType, false, f->getApi());

  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  for (E_Int v = 1; v <= nfld; v++)
  {
    E_Float* fp1 = f->begin(v);
    E_Float* fp2 = f2->begin(v);

    for (E_Int i = 0; i < npts; i++) fp2[i] = fp1[i];
    fp2[npts] = fp1[N1];
    if (N2 > 0) fp2[npts+1] = fp1[N2];
  }

  // copie connectivite
  for (E_Int n = 1; n <= cn2->getNfld(); n++)
    for (E_Int i = 0; i < cn2->getSize(); i++)
      (*cn2)(i,n) = (*cn)(i,n);

  // change le pts
  for (E_Int i = 0; i < cn->getSize(); i++)
  {
    if ((*cn2)(i,1) == N1+1) { (*cn2)(i,1) = npts+1; }
    if ((*cn2)(i,1) == N2+1) { (*cn2)(i,1) = npts+2; }
  }

  RELEASESHAREDU(tpl, f2, cn2);
  RELEASESHAREDU(array, f, cn);
  return tpl;
}

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
# include "converter.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
PyObject* K_CONVERTER::convertBAR2Struct(PyObject* self, PyObject* args)
{
  E_Float eps = 1.e-10;
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, nil, njl, nkl, 
                                    cn, eltType);
  if (res != 2)
  {
    if (res == 1) delete f;
    return array;
  }
  if (K_STRING::cmp(eltType, "BAR") != 0)
  {
    PyErr_SetString(PyExc_TypeError, "convertBAR2Struct: array must be of BAR type.");
    delete f; delete cn; return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
  if (posx == 0 || posy == 0 || posz == 0)
  {
    PyErr_SetString(PyExc_TypeError, "convertBAR2Struct: array must contain coordinates.");
    delete f; delete cn; return NULL;
  }
  K_CONNECT::cleanConnectivity(posx, posy, posz, eps, "BAR", *f, *cn);
  E_Int npts = f->getSize(); E_Int nelts = cn->getSize();
  E_Int nfld = f->getNfld();
  E_Int ni = npts; E_Int nj = 1; E_Int nk = 1;
  if (nelts != npts-1) ni = ni+1; // BAR = loop
  FldArrayF* fout = new FldArrayF(ni, nfld);
  K_CONNECT::orderBAR2Struct(posx, posy, posz, *f, *cn, *fout);
  delete f; delete cn;
  ni = fout->getSize();
  PyObject* tpl = K_ARRAY::buildArray(*fout, varString, ni,nj,nk);
  delete fout; 
  return tpl;
}

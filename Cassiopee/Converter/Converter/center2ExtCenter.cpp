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
#include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Conversion d un array en centres en centres �tendus 
   Extrapolation des centres sur les faces externes
   Attention : si les coordonnees sont fournies, sont aussi extrapol�es ! */
//=============================================================================
PyObject* K_CONVERTER::center2ExtCenter(PyObject* self, PyObject* args )
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* FCenter; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, FCenter, ni, nj, nk, 
                                     cn, eltType); 
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError, "center2ExtCenter: array must be structured."); 
    if (res == 2) RELEASESHAREDU(array, FCenter, cn);
    return NULL; 
  }
  E_Int nfld = FCenter->getNfld();
  E_Int nie = ni+2; E_Int nje = nj+2; E_Int nke = nk+2;
  PyObject* tpl = K_ARRAY::buildArray(nfld, varString, nie, nje, nke);
  E_Float* fep = K_ARRAY::getFieldPtr(tpl);
  FldArrayF FExtCenter(nie*nje*nke, nfld, fep, true);
  E_Int ret = K_LOC::center2ExtCenterStruct(ni, nj, nk, *FCenter, nie, nje, nke, FExtCenter);
  RELEASESHAREDS(array, FCenter);
  if (ret == 0) return NULL;
  return tpl;
}

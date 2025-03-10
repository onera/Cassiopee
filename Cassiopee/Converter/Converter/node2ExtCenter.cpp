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
/* Conversion d un array en noeuds en centres �tendus*/
//=============================================================================
PyObject* K_CONVERTER::node2ExtCenter(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* FNode; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, FNode, ni, nj, nk, 
                                    cn, eltType, true); 
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError, "node2ExtCenter: array must be structured."); 
    if (res == 2) RELEASESHAREDU(array, FNode, cn);
    return NULL; 
  }
  E_Int nfld = FNode->getNfld();
  E_Int ni2 = ni; E_Int nj2 = nj; E_Int nk2 = nk;
  if (ni == 1) {ni2 = nj; nj2 = nk; nk2 = 1;}
  else if (nj == 1) {ni2 = ni; nj2 = nk; nk2 = 1;}  
  E_Int nie = ni2+1, nje = nj2+1, nke = nk2+1;
  if (ni2 == 1) nie = 1;
  if (nj2 == 1) nje = 1;
  if (nk2 == 1) nke = 1;

  PyObject* tpl = K_ARRAY::buildArray(nfld, varString, nie, nje, nke);
  E_Float* fep = K_ARRAY::getFieldPtr(tpl);
  FldArrayF FExtCenter(nie*nje*nke, nfld, fep, true);
  E_Int ret = K_LOC::node2ExtCenterStruct(ni2, nj2, nk2, *FNode, nie, nje, nke, FExtCenter);
  RELEASESHAREDS(array, FNode);
  if (ret == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "node2ExtCenter: Fail to compute extended center mesh."); 
    return NULL;
  }
  return tpl;
}

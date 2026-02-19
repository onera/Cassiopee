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
#include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Conversion d un array en centres �tendus en noeuds */
//=============================================================================
PyObject* K_CONVERTER::extCenter2Node(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;
  
  // Check array
  E_Int nie, nje, nke;
  FldArrayF* fc; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, fc, 
                                     nie, nje, nke, cn, eltType); 
  if (res != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "node2ExtCenter: array must be structured."); 
    if (res == 2) RELEASESHAREDU(array, fc, cn);
    return NULL; 
  }

  E_Int nfld = fc->getNfld();
  E_Int ni = nie-1; E_Int nj = nje-1; E_Int nk = nke-1;
  if (ni == 0) ni = 1;
  if (nj == 0) nj = 1;
  if (nk == 0) nk = 1;

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, ni, nj, nk);
  FldArrayF* fn2;
  K_ARRAY::getFromArray3(tpl, fn2);

  E_Int ret = K_LOC::extCenters2NodeStruct(nie, nje, nke, *fc, 
                                           ni, nj, nk, *fn2);
  RELEASESHAREDS(array, fc);
  RELEASESHAREDS(tpl, fn2);

  if (ret == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "extCenter2Node: Fail to compute node mesh."); 
    return NULL;
  }
  return tpl;
}

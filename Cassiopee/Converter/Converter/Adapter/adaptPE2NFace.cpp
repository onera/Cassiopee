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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert a NGon parentElement numpy to a NGon NFace numpy */
//=============================================================================
PyObject* K_CONVERTER::adaptPE2NFace(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int ngonType;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &ngonType)) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  E_Int res = K_NUMPY::getFromNumpyArray(array, cFE);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "adaptPE2NFace: numpy is invalid.");
    return NULL;
  }
  
  FldArrayI cNFace, off; E_Int nelts;
  if (ngonType <= 3) K_CONNECT::connectFE2NFace3(*cFE, cNFace, off, nelts);
  else K_CONNECT::connectFE2NFace4(*cFE, cNFace, off, nelts);
  
  RELEASESHAREDN(array, cFE);

  PyObject* tpl = K_NUMPY::buildNumpyArray(cNFace, false);
  PyObject* tplO = K_NUMPY::buildNumpyArray(off, false);
  PyObject* tpl2 = Py_BuildValue("OOl", tpl, tplO, nelts);
  return tpl2;
}


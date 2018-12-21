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

using namespace K_FLD;

//=============================================================================
/* Convert a NGon parentElement numpy to a NGon NFace numpy */
//=============================================================================
PyObject* K_CONVERTER::adaptPE2NFace(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check numpy (parentElement)
  FldArrayI* cFE;
  E_Int res = K_NUMPY::getFromNumpyArray(array, cFE, true);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptPE2NFace: numpy is invalid.");
    return NULL;
  }
  
  FldArrayI cNFace; E_Int nelts;
  K_CONNECT::connectFE2NFace(*cFE, cNFace, nelts);
  
  RELEASESHAREDN(array, cFE);

  // Construction du numpy de sortie
  PyObject* tpl = K_NUMPY::buildNumpyArray(cNFace, false);
  PyObject* tpl2 = Py_BuildValue("Ol", tpl, nelts);
  return tpl2;
}

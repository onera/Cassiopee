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
/* Convert a NGon numpy to FaceIndex numpy */
//=============================================================================
PyObject* K_CONVERTER::adaptNGon2Index(PyObject* self, PyObject* args)
{
  PyObject* arrayNG; E_Int nfaces;
  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arrayNG, &nfaces)) return NULL;

  // Check numpy (NGon)
  FldArrayI* cNGon;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayNG, cNGon, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNGon2Index: cNGon numpy is invalid.");
    return NULL;
  }

  PyObject* tpl = K_NUMPY::buildNumpyArray(nfaces, 1, 1, 1);
  E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
  K_CONNECT::getIndex(cNGon->begin(), 0, nfaces, pos);

  RELEASESHAREDN(arrayNG, cNGon);
  return tpl;
}

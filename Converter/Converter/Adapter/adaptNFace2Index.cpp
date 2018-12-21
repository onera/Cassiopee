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
/* Convert a NFace numpy to FaceIndex numpy */
//=============================================================================
PyObject* K_CONVERTER::adaptNFace2Index(PyObject* self, PyObject* args)
{
  PyObject* arrayNF; E_Int nelts;
  if (!PYPARSETUPLEI(args, "Ol", "Oi", &arrayNF, &nelts)) return NULL;

  // Check numpy (NFace)
  FldArrayI* cNFace;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayNF, cNFace, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNFace2Index: cNFace numpy is invalid.");
    return NULL;
  }

  PyObject* tpl = K_NUMPY::buildNumpyArray(nelts, 1, 1, 1);
  E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
  K_CONNECT::getIndex(cNFace->begin(), 0, nelts, pos);

  RELEASESHAREDN(arrayNF, cNFace);
  return tpl;
}

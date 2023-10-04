/*    
    Copyright 2013-2023 Onera.

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
/* Sign faces in an NGon */
/* IN: NGonv3 or NGonv4 */
//=============================================================================
PyObject* K_CONVERTER::signNGonFaces(PyObject* self, PyObject* args)
{
  PyObject* o; 
  if (!PYPARSETUPLE_(args, O_, &o)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray3(o, varString, f, ni, nj, nk, c, eltType);
  
  if (ret <= 0)
  { PyErr_SetString(PyExc_TypeError, "signNGonFaces: only for NGons."); return NULL; }

  if (ret == 1)
  { 
    PyErr_SetString(PyExc_TypeError, "signNGonFaces: only for NGons."); 
    RELEASESHAREDS(o, f);
    return NULL; 
  }

  RELEASESHAREDU(o, f, c);

  return NULL;
}

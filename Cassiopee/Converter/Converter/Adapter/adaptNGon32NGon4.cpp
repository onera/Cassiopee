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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert a NGon1 (NGONv3) to NGon2 connectivity+offset (NGONv4) */
//=============================================================================
PyObject* K_CONVERTER::adaptNGon32NGon4(PyObject* self, PyObject* args)
{
  PyObject* arrayConnect;
  if (!PYPARSETUPLE_(args, O_, &arrayConnect)) return NULL;

  // Check numpy (connect)
  FldArrayI* connect;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayConnect, connect, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNGon32NGon4: connect numpy is invalid.");
    return NULL;
  }
  
  // Find the size and number of elements
  E_Int* c = connect->begin();
  E_Int np;
  E_Int ne = 0; E_Int size = 0;
  E_Int sizeTot = connect->getSize();
  E_Int p = 0; E_Int p2 = 0;
  while (p < sizeTot)
  {
    np = c[0];
    size += np;
    ne++;
    p += np+1;
    c += np+1;
  }
  
  // creation sortie
  PyObject* tpl = K_NUMPY::buildNumpyArray(size, 1, 1, 1);
  E_Int* co = K_NUMPY::getNumpyPtrI(tpl);
  PyObject* tpl2 = K_NUMPY::buildNumpyArray(ne+1, 1, 1, 1);
  E_Int* off = K_NUMPY::getNumpyPtrI(tpl2);
  
  // remplissage
  p = 0; p2 = 0; c = connect->begin(); off[0] = 0;
  for (E_Int i = 0; i < ne; i++)
  {
      np = c[0];
      for (E_Int j = 0; j < np; j++) co[j] = c[j+1];
      p2 += np;
      off[i+1] = p2;
      c += np+1;
      co += np;
  }
 
  RELEASESHAREDN(arrayConnect, connect);

  // Retour du numpy de sortie
  return Py_BuildValue("(OO)", tpl,tpl2);
}

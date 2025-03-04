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
/* Convert a NGon2+ElementStartOffSet (NGONv4) to NGON1 (NGONv3)  */
//=============================================================================
PyObject* K_CONVERTER::adaptNGon42NGon3(PyObject* self, PyObject* args)
{
  PyObject* arrayConnect; PyObject* arrayOffset;
  if (!PYPARSETUPLE_(args, OO_, 
                     &arrayConnect, &arrayOffset)) return NULL;

  // Check numpy (connect)
  FldArrayI* connect;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayConnect, connect, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNGon42NGon3: connect numpy is invalid.");
    return NULL;
  }
  // Check numpy (Offset)
  FldArrayI* offset;
  res = K_NUMPY::getFromNumpyArray(arrayOffset, offset, true);
  if (res == 0)
  {
    RELEASESHAREDN(arrayConnect, connect);
    PyErr_SetString(PyExc_TypeError, 
                    "adaptNGon42NGon3: offset numpy is invalid.");
    return NULL;
  }
  
  // nbre d'element dans la connectivite
  E_Int ne = offset->getSize()-1;
  // taille de la connectivite de sortie contenant le no + les indices pour chaque element
  E_Int size = connect->getSize()+ne;

  // creation sortie
  PyObject* tpl = K_NUMPY::buildNumpyArray(size, 1, 1, 1);
  E_Int* co = K_NUMPY::getNumpyPtrI(tpl);

  // remplissage
  E_Int* c = connect->begin();
  E_Int* o = offset->begin();
#pragma omp parallel
  {
    E_Int d, d1, d2;
    #pragma omp for
    for (E_Int i = 0; i < ne; i++)
    {
      d = o[i+1]-o[i];
      d1 = o[i];
      d2 = o[i]+i;
      co[d2] = d;
      for (E_Int j = 0; j < d; j++) co[d2+j+1] = c[d1+j];
    }
  }

  // Create index array from offset
  for (E_Int i = 0; i < ne; i++) o[i] = o[i]+i;
 
  RELEASESHAREDN(arrayConnect, connect);
  RELEASESHAREDN(arrayOffset, offset);

  // Retour du numpy de sortie
  return tpl;
}

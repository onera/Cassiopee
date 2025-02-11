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
#include "Connect/connect.h"

using namespace K_FLD;

//=============================================================================
/* Sign faces in an NGon */
/* IN: NGonv3 or NGonv4 */
//=============================================================================
PyObject* K_CONVERTER::makeParentElements(PyObject* self, PyObject* args)
{
  PyObject* o;
  if (!PYPARSETUPLE_(args, O_, &o)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray3(o, varString, f, ni, nj, nk, cn, eltType);
  
  if (ret <= 0) {
    PyErr_SetString(PyExc_TypeError, "makeParentElements: bad mesh.");
    return NULL;
  }

  if (ret == 1) { 
    PyErr_SetString(PyExc_TypeError, "makeParentElements: only for NGons."); 
    RELEASESHAREDS(o, f);
    return NULL; 
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(ret, o, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "makeParentElements: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // Coordonnees
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  E_Int *nface = cn->getNFace();
  E_Int *indPH = cn->getIndPH();
  E_Int ncells = cn->getNElts();

  // Reset sign of all faces
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = cn->getElt(i, stride, nface, indPH);
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j];
      if (face < 0) pf[j] = -face;
    }
  }

  // Check mesh
  ret = K_CONNECT::check_overlapping_cells(*cn);
  if (ret == 1) {
    PyErr_SetString(PyExc_ValueError,
                    "makeParentElements: non-valid mesh.");
    return NULL;
  }

  // Orient external faces outwards
  K_CONNECT::orient_boundary_ngon(x, y, z, *cn);
  
  // Deduce parent elements
  E_Int nfaces = cn->getNFaces();
  PyObject *pe = K_NUMPY::buildNumpyArray(nfaces, 2, 1, 1);
  E_Int *ptr = K_NUMPY::getNumpyPtrI(pe);
  E_Int *owner = ptr;
  E_Int *neigh = ptr + nfaces;
  K_CONNECT::build_parent_elements_ngon(*cn, owner, neigh);
 
  for (E_Int i = 0; i < 2*nfaces; i++) ptr[i] += 1;

  RELEASESHAREDU(o, f, cn);

  return (PyObject *)pe;
}

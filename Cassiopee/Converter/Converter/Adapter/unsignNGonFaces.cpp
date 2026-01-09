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
# include "kcore.h"

using namespace K_FLD;

//=============================================================================
/* Unsign faces in an NGon v4 */
//=============================================================================
PyObject* K_CONVERTER::unsignNGonFaces(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f,
                                     ni, nj, nk, cn, eltType);
  
  // Test non structure
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, "unsignNGonFaces: invalid array.");
    return NULL;
  }

  if (res == 1)
  { 
    PyErr_SetString(PyExc_TypeError,
                    "unsignNGonFaces: input array must be unstructured."); 
    RELEASESHAREDS(array, f);
    return NULL;
  }
  
  // Test type
  if (strcmp(eltType, "NGON") != 0 and strcmp(eltType, "NGON*") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "unsignNGonFaces: unstructured array must be NGon.");
    return NULL;
  }
  
  E_Int nf;
  E_Int ncells = cn->getNElts();
  E_Int sizeEF = cn->getSizeNFace();
  E_Int* nface = cn->getNFace();
  E_Int* indPH = cn->getIndPH();
  
  E_Int isSigned = 0;
  ncells = K_FUNC::E_min(ncells, 100); // search for a negative sign in a few elts
  
  for (E_Int i = 0; i < ncells; i++)
  {
    E_Int* elt = cn->getElt(i, nf, nface, indPH);
    for (E_Int j = 0; j < nf; j++)
    {
      if (elt[j] < 0) { isSigned = 1; break; }
    }
    if (isSigned == 1) break;
  }
  
  if (isSigned == 1)
  {
    // Unsign all faces
    #pragma omp parallel for
    for (E_Int i = 0; i < sizeEF; i++) nface[i] = K_FUNC::E_abs(nface[i]);
  }
  
  RELEASESHAREDU(array, f, cn);
  return Py_BuildValue("i", isSigned);
}

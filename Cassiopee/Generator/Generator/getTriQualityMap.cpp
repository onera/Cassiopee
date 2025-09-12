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

// getTriQualityMap

# include "generator.h"
#include "Nuga/include/Triangle.h"

using namespace K_FLD;

// ============================================================================
/* Return TRI quality map */
// Definition of the returned value: Q = 4 * sqrt(3) * S / (L * p), S suface, L longest edge, p perimeter.
// Q=0 for a flat triangle, Q=1 for a equilateral triangle.
// ============================================================================
PyObject* K_GENERATOR::getTriQualityMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, 
                               eltType);

  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getTriQualityMap: unknown type of array.");
    return NULL;
  }

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "getTriQualityMap: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // pointeurs associes aux coordonnees
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  
  E_Int nelts = cn->getSize();
  E_Int nnodes = cn->getNfld(); // nb de noeuds ds 1 element
  
  PyObject* tpl = K_ARRAY::buildArray(1, "quality", nelts, 
					nelts, -1, eltType, true);
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  K_KCORE::memcpy__(cnnp, cn->begin(), nelts*nnodes);
    
  E_Float* qualitiesp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF qualities(nelts,1, qualitiesp, true);
  
  if (strcmp(eltType, "TRI") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                      "getTriQualityMap: must be TRI.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }

  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3);
  E_Int ind1, ind2, ind3;
  E_Float P0[3], P1[3], P2[3];
  //
  for (E_Int et = 0; et < nelts; et++)
  {
    ind1 = cn1[et]-1;
    ind2 = cn2[et]-1;
    ind3 = cn3[et]-1;
    
    P0[0] = x[ind1]; P0[1] = y[ind1]; P0[2] = z[ind1];
    P1[0] = x[ind2]; P1[1] = y[ind2]; P1[2] = z[ind2];
    P2[0] = x[ind3]; P2[1] = y[ind3]; P2[2] = z[ind3];
    
    qualities[et]=K_MESH::Triangle::qualityG<3>(&P0[0], &P1[0], &P2[0]);
  }
  
  RELEASESHAREDU(array, f, cn); 
  return tpl;
}

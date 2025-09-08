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

# include "generator.h"
using namespace K_CONST;
using namespace K_FLD;

//============================================================================
/* Retourne le rayon du cercle circonscrit de tous les triangles d un array */
//============================================================================
PyObject* K_GENERATOR::getCircumCircleMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  if ( !PYPARSETUPLE_(args, O_, &array) )
  {
    return NULL;
  }
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 2 || strcmp(eltType, "TRI") != 0) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "getCircumCircleMap: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "getCircumCircleMap: cannot find coordinates in array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  
  E_Int ncells = cn->getSize();
  E_Int nnodes = cn->getNfld(); // nb de noeuds ds 1 element

  PyObject* tpl = K_ARRAY::buildArray(1, "ccradius", ncells, ncells, 
				      -1, eltType, true);
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  K_KCORE::memcpy__(cnnp, cn->begin(), ncells*nnodes);
  E_Float* ccradp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF ccrad(ncells, 1, ccradp, true);

  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  
  E_Int ind1, ind2, ind3;
  E_Float p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z;
  for (E_Int et = 0; et < ncells; et++)
  {
    ind1 = cn1[et]-1; ind2 = cn2[et]-1; ind3 = cn3[et]-1;
    p1x = xt[ind1]; p1y = yt[ind1]; p1z = zt[ind1];  
    p2x = xt[ind2]; p2y = yt[ind2]; p2z = zt[ind2];
    p3x = xt[ind3]; p3y = yt[ind3]; p3z = zt[ind3];
    ccrad[et] = K_COMPGEOM::circumCircleRadius(p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z);
  }
  
  RELEASESHAREDU(array, f, cn); 
  return tpl;
}

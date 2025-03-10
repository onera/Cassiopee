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
/* Retourne le rayon du cercle inscrit de tous les triangles d'un array */
//============================================================================
PyObject* K_GENERATOR::getInCircleMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true);

  if (res != 2 || strcmp(eltType, "TRI") != 0) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "getInCircleMap: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "getInCircleMap: cannot find coordinates in array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  
  E_Int ncells = cn->getSize();
  E_Int nnodes = cn->getNfld(); // nb de noeuds ds 1 element
  PyObject* tpl = K_ARRAY::buildArray(1, "icradius", ncells, ncells, -1, eltType, true);
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  K_KCORE::memcpy__(cnnp, cn->begin(), ncells*nnodes);
  E_Float* ccradp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF ccrad(ncells,1, ccradp, true);

  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  E_Int* cn1 = cn->begin(1);
  E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3);
  
  E_Float p1[3]; E_Float p2[3]; E_Float p3[3];
  for (E_Int et = 0; et < ncells; et++)
  {
    
    E_Int ind1 = cn1[et]-1; E_Int ind2 = cn2[et]-1; E_Int ind3 = cn3[et]-1;
    p1[0] = xt[ind1]; p1[1] = yt[ind1]; p1[2] = zt[ind1];  
    p2[0] = xt[ind2]; p2[1] = yt[ind2]; p2[2] = zt[ind2];  
    p3[0] = xt[ind3]; p3[1] = yt[ind3]; p3[2] = zt[ind3];
    ccrad[et] = K_COMPGEOM::inscribedCircleRadius(p1, p2, p3);
  }
  RELEASESHAREDU(array, f, cn); 
  return tpl;
}

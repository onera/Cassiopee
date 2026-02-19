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

// getEdgeRatio

# include "generator.h"

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Return edge ratio */
// ============================================================================
PyObject* K_GENERATOR::getEdgeRatio(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int dimPb;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &dimPb)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, 
                                     eltType);
  PyObject* tpl = NULL;
  E_Int dim = E_Int(dimPb);
  if (dim != 2 && dim != 3)
  {
    PyErr_SetString(PyExc_ValueError,
                    "getEdgeRatio: dim must be equal to 2 or 3.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  if (res != 1 && res != 2)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "getEdgeRatio: unknown type of array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "getEdgeRatio: cannot find coordinates in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  // Coordonnees
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  E_Int api = f->getApi();

  if (res == 1) // cas structure
  {
    E_Int im1 = im-1;
    E_Int jm1 = jm-1;
    E_Int km1 = km-1;
    if (im == 1) im1 = 1;
    if (jm == 1) jm1 = 1;
    if (km == 1) km1 = 1;
    tpl = K_ARRAY::buildArray3(1, "EdgeRatio", im1, jm1, km1, api);
    FldArrayF* field; 
    K_ARRAY::getFromArray3(tpl, field);

    E_Int ret = K_COMPGEOM::getEdgeLength(xt, yt, zt,
                                          im, jm, km, 
                                          NULL, NULL,
                                          dim, 2,
                                          field->begin());
    RELEASESHAREDS(tpl, field);
    RELEASESHAREDS(array, f);
    if (ret == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getEdgeRatio: failed.");
        return NULL;
    }
    return tpl;
  }
  else
  {
    tpl = K_ARRAY::buildArray3(1, "EdgeRatio", f->getSize(), *cn, eltType, 1, api);
    FldArrayF* field; 
    K_ARRAY::getFromArray3(tpl, field);
    E_Int ret = K_COMPGEOM::getEdgeLength(
      xt, yt, zt, -1, -1, -1, 
      cn, eltType, dim, 2, field->begin()
    );

    RELEASESHAREDS(tpl, field);
    RELEASESHAREDU(array, f, cn);
    if (ret == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getEdgeRatio: failed.");
      return NULL;          
    }
    return tpl; 
  }
}

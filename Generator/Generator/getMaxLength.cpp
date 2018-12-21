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

// getMaxLength

# include "generator.h"

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Return max length of a cell */
// ============================================================================
PyObject* K_GENERATOR::getMaxLength(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int dimPb;
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Ol", &array, &dimPb)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oi", &array, &dimPb)) return NULL;
#endif
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, 
                                    eltType, true);
  PyObject* tpl = NULL;
  E_Int dim = E_Int(dimPb);
  if (dim != 2 && dim != 3) 
  {
    PyErr_SetString(PyExc_ValueError,
                    "getMaxLength: dim must be equal to 2 or 3.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  if (res != 1 && res != 2)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "getMaxLength: unknown type of array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "getMaxLength: cannot find coordinates in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  // Coordonnees
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
 
  if (res == 1) // cas structure
  {
    E_Int im1 = im-1;
    E_Int jm1 = jm-1;
    E_Int km1 = km-1;
    if (im == 1) im1 = 1;
    if (jm == 1) jm1 = 1;
    if (km == 1) km1 = 1;
      
    tpl = K_ARRAY::buildArray(1, "MaxLength", im1, jm1, km1);
    E_Float* fieldp = K_ARRAY::getFieldPtr(tpl);

    E_Int ret = K_COMPGEOM::getEdgeLength(xt, yt, zt,
                                          im, jm, km, 
                                          NULL, NULL,
                                          dim, 0,
                                          fieldp);
    RELEASESHAREDS(array, f);
    if (ret == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                        "getMaxLength: failed.");
      return NULL;
    }
    return tpl;
  }
  else
  {
    E_Int nelts;
    if (strcmp(eltType,"NGON") == 0) 
    { E_Int* cnp = cn->begin(); E_Int sizeFN = cnp[1]; nelts = cnp[sizeFN+2]; }
    else nelts = cn->getSize();
    tpl = K_ARRAY::buildArray(1, "MaxLength",
                              f->getSize(), nelts, -1, eltType, true,
                              cn->getSize());
    E_Float* fieldp = K_ARRAY::getFieldPtr(tpl);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI cnn(nelts, cn->getNfld(), cnnp, true); cnn = *cn;
    
    E_Int ret = K_COMPGEOM::getEdgeLength(xt, yt, zt,
                                          -1, -1, -1, 
                                          cn, eltType,
                                          dim, 0, fieldp);

    RELEASESHAREDU(array, f, cn);
    if (ret == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                        "getMaxLength: failed.");
      return NULL;          
    }
    return tpl;
  }
}

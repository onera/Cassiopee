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

// Information on geometries

# include <string.h>
# include "geom.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Get the curvature height of all points */
// ============================================================================
PyObject* K_GEOM::getCurvatureHeight(PyObject* self,
                                     PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f;
  FldArrayI* cn;
  char* varString;
  char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, 
                                    im, jm, km, cn, eltType, true);
  if ( res != 1 && res != 2 ) 
  {
    PyErr_SetString(PyExc_TypeError, "getCurvatureHeight: invalid array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if ( posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getCurvatureHeight: array must contain coordinates.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  if ( res == 1 ) //1D ou 2D seulement
  {
    if ( im > 1 && jm > 1 && km > 1 ) 
    {
      PyErr_SetString(PyExc_TypeError,
                    "getCurvatureHeight: array must be 1D or 2D.");
      RELEASESHAREDS(array, f); return NULL;
    }
  }
  else if ( res == 2 ) 
  {
    if ( strcmp(eltType,"BAR") != 0 && strcmp(eltType,"TRI") != 0 && strcmp(eltType,"QUAD") != 0 )
    {
      PyErr_SetString(PyExc_TypeError,
                    "getCurvatureHeight: unstructured array must be BAR, TRI or QUAD.");
      RELEASESHAREDU(array, f, cn); return NULL;
    }
  }

  // build array
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);
  PyObject* tpl = NULL;

  if (res == 1)
  {
    tpl = K_ARRAY::buildArray(1,"hmax", im, jm, km);     
    E_Float* hmaxt = K_ARRAY::getFieldPtr(tpl);                                  
    if (im>1 && jm==1 && km== 1) K_COMPGEOM::compStructCurvatureHeight1D(im, xt, yt, zt, hmaxt);
    else if (im>1 && jm>1  && km==1 ) K_COMPGEOM::compStructCurvatureHeight2D(im, jm, xt, yt, zt, hmaxt);
    else 
    {
      PyErr_SetString(PyExc_TypeError,
                    "getCurvatureHeight: array must be (ni,1,1) in 1D or (ni,nj,1) in 2D.");
      RELEASESHAREDS(array, f); return NULL;
    }
    RELEASESHAREDS(array, f);
  }
  else
  {
    tpl = K_ARRAY::buildArray(1, "hmax", f->getSize(), cn->getSize(), -1, eltType);
    E_Float* hmaxt = K_ARRAY::getFieldPtr(tpl);                                  
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI cnn(cn->getSize(), cn->getNfld(), cnnp, true); cnn = *cn;
    
    if (strcmp(eltType,"BAR") == 0) K_COMPGEOM::compCurvatureHeightForBAR(f->getSize(), xt, yt, zt, cnn, hmaxt);
    else if (strcmp(eltType,"TRI") == 0 || strcmp(eltType,"QUAD") == 0) K_COMPGEOM::compCurvatureHeightForTRIQUAD(f->getSize(), xt, yt, zt, cnn, hmaxt);
    RELEASESHAREDU(array, f, cn);
  }
  return tpl; 
}     

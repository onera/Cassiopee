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

using namespace K_FLD;
 
//=============================================================================
/* Returns the bounding box of each cell of an array */
//=============================================================================
PyObject* K_GENERATOR::getBBOfCells(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, 
                                    ni, nj, nk, cn, eltType);
  
  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "getBBoxOfCells: invalid array.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getBBoxOfCells: cannot find coordinates in array.");
    return NULL;        
  }
  posx++; posy++; posz++;
  
  if (res == 1) 
  {
    E_Int ni1 = K_FUNC::E_max(1,ni-1);
    E_Int nj1 = K_FUNC::E_max(1,nj-1);
    E_Int nk1 = K_FUNC::E_max(1,nk-1);
    PyObject* tpl = K_ARRAY::buildArray3(6, "xmin,ymin,zmin,xmax,ymax,zmax", 
      ni1, nj1, nk1);

    char* varString2; FldArrayF* f2;
    K_ARRAY::getFromArray3(tpl, varString2, f2);

    K_COMPGEOM::boundingBoxOfStructCells(ni, nj, nk, *f, *f2);
    RELEASESHAREDS(array, f);
    RELEASESHAREDS(tpl, f2);
    return tpl;
  }
  else 
  {
    if (strcmp(eltType, "NGON") == 0)
    {
      E_Int api = f->getApi();
      E_Int nvertex = f->getSize();

      PyObject* tpl = K_ARRAY::buildArray3(6, "xmin,ymin,zmin,xmax,ymax,zmax", 
        nvertex, *cn, "NGON", 1, 3, true); // center=1, api=3, copyConnect=true

      FldArrayF* f2;  FldArrayI* cn2;
      K_ARRAY::getFromArray3(tpl, f2, cn2);

      K_COMPGEOM::boundingBoxOfNGonCells(*cn2, f->begin(posx), f->begin(posy),
                                          f->begin(posz), *f2);
      RELEASESHAREDU(array, f, cn);
      RELEASESHAREDU(tpl, f2, cn2);
      return tpl;
    }
    else
    {
      E_Int nvertex = f->getSize();
      
      PyObject* tpl = K_ARRAY::buildArray3(6, "xmin,ymin,zmin,xmax,ymax,zmax", 
        nvertex, *cn, eltType, 1, 3, true); // center=1, api=3, copyConnect=true

      FldArrayF* f2;  FldArrayI* cn2;
      K_ARRAY::getFromArray3(tpl, f2, cn2);

      K_COMPGEOM::boundingBoxOfUnstrCells(*cn2, f->begin(posx), f->begin(posy),
                                          f->begin(posz), *f2);
      RELEASESHAREDU(array, f, cn);
      RELEASESHAREDU(tpl, f2, cn2);
      return tpl;
    }
  }
}

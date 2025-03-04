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
  E_Int res = K_ARRAY::getFromArray(array, varString, f, 
                                    ni, nj, nk, cn, eltType, true);
  
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
    PyObject* tpl = K_ARRAY::buildArray(6, "xmin,ymin,zmin,xmax,ymax,zmax", 
					ni1, nj1, nk1);
    E_Float* bbp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF bboxOfCells(ni1*nj1*nk1, 6, bbp, true);
    K_COMPGEOM::boundingBoxOfStructCells(ni, nj, nk, *f, bboxOfCells);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else 
  {
    E_Int net = cn->getSize(); E_Int nvert = cn->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(6, "xmin,ymin,zmin,xmax,ymax,zmax", 
					net, net, -1, eltType, true);
    E_Float* bbp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF bboxOfCells(net, 6, bbp, true);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
    FldArrayI cnn(net, nvert, cnnp, true);
    K_COMPGEOM::boundingBoxOfUnstrCells(cnn, f->begin(posx), f->begin(posy),
                                        f->begin(posz), bboxOfCells);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
}

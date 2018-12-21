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

// Information on grids: check if point is in CEBB

# include "generator.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Check if a given point is in the Cartesian Elements Bounding Box of a 
   structured mesh */
// ============================================================================
PyObject* K_GENERATOR::checkPointInCEBBOfMesh(PyObject* self,
                                              PyObject* args)
{
  PyObject* array;
  double x, y, z;
  if (!PyArg_ParseTuple(args, "O(ddd)", &array, &x, &y, &z)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  FldArrayF coord;
  FldArrayF cartEltArray;// alloue ds compCartEltsArray
  E_Float *xtm, *ytm, *ztm, *xtp, *ytp, *ztp;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true);

  E_Int npts;
  short isok;
  E_Float x1, y1, z1, x2, y2, z2;
  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "checkPointInCEBB: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    // Calcul de la bounding box de l'array
    E_Int found = 0;
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    K_COMPGEOM::boundingBox(im, jm, km, posx, posy, posz, *f,
                            xmin, ymin,  zmin, xmax, ymax,zmax);

    if (x > xmax || y > ymax || z > zmax ||
        x < xmin || y < ymin || z < zmin) goto end;
    npts = f->getSize();
    coord.malloc(npts, 3);
    coord.setOneField(*f, posx, 1);
    coord.setOneField(*f, posy, 2);
    coord.setOneField(*f, posz, 3);

    // Construction de la CEBB
    for (E_Int dir = 1; dir <= 3; dir++) 
    {
      isok = K_COMPGEOM::compCartEltsArray(dir, im, jm, km,
                                           xmin, ymin, zmin, 
                                           xmax, ymax, zmax,
                                           coord, cartEltArray);
      if (isok == 1) break;
    }
    
    // Check if given point is in CEBB
    xtm = cartEltArray.begin(1);
    ytm = cartEltArray.begin(2);
    ztm = cartEltArray.begin(3);
    xtp = cartEltArray.begin(4);
    ytp = cartEltArray.begin(5);
    ztp = cartEltArray.begin(6);   
    for (E_Int e1 = 0; e1 < cartEltArray.getSize(); e1++)
    {
      x1 = xtm[e1]; y1 = ytm[e1]; z1 = ztm[e1];
      x2 = xtp[e1]; y2 = ytp[e1]; z2 = ztp[e1];
      if (x <= x2 && y <= y2 && z <= z2 &&
          x >= x1 && y >= y1 && z >= z1) 
      {
        found = 1;
        goto end;
      }
    }
    
    end:;
    RELEASESHAREDS(array, f);
#ifdef E_DOUBLEINT
    return Py_BuildValue("l", long(found));
#else
    return Py_BuildValue("i", found);
#endif
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "checkPointInCEBB: not for unstructured arrays.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "checkPointInCEBB: unrecognised type of array.");
    return NULL;
  }
}

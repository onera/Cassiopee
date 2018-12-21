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
# include <string.h>
# include "post.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Calcul un perlin noise en fonction de x,y,z */
//=============================================================================
PyObject* K_POST::perlinNoise(PyObject* self,PyObject* args)
{
  PyObject* array;
  E_Float alpha, beta;
  E_Int n;
  if (!PYPARSETUPLE(args,
                    "Oddl", "Oddi",
                    "Offl", "Offi",
                    &array, &alpha, &beta, &n))
  {
      return NULL;
  }
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array

  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, cn, 
                                    eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "perlin: invalid array.");
    return NULL;
  }

  //check validity alpha, beta, n
  
  // Check coordinates
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "perlin: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  char varStringOut[256];
  strcpy(varStringOut, varString);
  strcat(varStringOut, ",perlin");

  posx++; posy++; posz++;
  E_Int nfld = f->getNfld()+1;
  E_Int npts = f->getSize();
  PyObject* tpl;
  if (res == 1) //structured
  {
    tpl = K_ARRAY::buildArray(nfld, varStringOut, 
                              ni, nj, nk);
  } 
  else //unstructured 
  {
    E_Int csize = cn->getSize()*cn->getNfld(); 
    tpl = K_ARRAY::buildArray(nfld, varStringOut,
                              npts, cn->getSize(),
                              -1, eltType, false, csize);
  }

  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fn(npts, nfld, fnp, true);
  // Recopie
  E_Float* fp = f->begin();
  for (E_Int v = 0; v < nfld-1; v++)
  {
    for (E_Int i = 0; i < npts; i++) fnp[i] = fp[i];
    fnp += npts; fp += npts;
  }
  // Perlin
  K_NOISE::PDS data;
  K_NOISE::initPerlinNoise(100, data);

  fnp = fn.begin(nfld);
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);
  E_Float xmin=1.e6, xmax=-1.e6, ymin=1.e6, ymax=-1.e6, zmin=1.e6, zmax=-1.e6;
  for (E_Int i = 0; i < npts; i++)
  {
    xmin = K_FUNC::E_min(xmin, x[i]);
    xmax = K_FUNC::E_max(xmax, x[i]);
    ymin = K_FUNC::E_min(ymin, y[i]);
    ymax = K_FUNC::E_max(ymax, y[i]);
    zmin = K_FUNC::E_min(zmin, z[i]);
    zmax = K_FUNC::E_max(zmax, z[i]);
  }
  E_Float dx = xmax-xmin; if (dx == 0.) { dx = 1.; xmin = 0; } else dx = 1./dx;
  E_Float dy = ymax-ymin; if (dy == 0.) { dy = 1.; ymin = 0; } else dy = 1./dy;
  E_Float dz = zmax-zmin; if (dz == 0.) { dz = 1.; zmin = 0; } else dz = 1./dz;
  
  for (E_Int i = 0; i < npts; i++)
  {
    fnp[i] = K_NOISE::perlinNoise3D( (x[i]-xmin)*dx, (y[i]-ymin)*dy, 
                                     (z[i]-zmin)*dz, alpha, beta, n, data);
  }

  if (res == 2)
  {
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());
  }

  RELEASESHAREDB(res, array,f, cn);
  return tpl;
}

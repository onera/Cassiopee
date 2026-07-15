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
#include "geom.h"

using namespace K_FLD;

#include "kadolc.h"

//===========================================================================
/* dLength: differentiated getLength */
//===========================================================================
PyObject* K_GEOM::dLength(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km,
                                     cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "dLength: invalid array.");
    return NULL;
  }

  if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "dLength: only for structured arrays.");
    return NULL;
  }

  E_Int dim = 0;
  if (im > 1) dim += 1;
  if (jm > 1) dim += 1;
  if (km > 1) dim += 1;
  if (dim != 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "dLength: only for 1D structured arrays.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "getLength: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // coordinates
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  // tangent vector
  E_Float* dxt = f->begin(4);
  E_Float* dyt = f->begin(5);
  E_Float* dzt = f->begin(6);

  E_Float length; // return length
  //typedef adouble E_Float;

  E_Int npts = im*jm*km;
  // point
  auto x = new E_Float [3*npts];
  for (E_Int i = 0; i < npts; i++)
  {
    x[3*i] = xt[i];
    x[3*i+1] = yt[i];
    x[3*i+2] = zt[i];
  }

  // tangent
  auto dx = new E_Float [3*npts];
  for (E_Int i = 0; i < npts; i++)
  {
    dx[3*i] = dxt[i];
    dx[3*i+1] = dyt[i];
    dx[3*i+2] = dzt[i];
  }

  // gradient
  auto df = new E_Float [3*npts];

  // actives
  auto ax = new adouble[3*npts];

  trace_on(0);

  for (E_Int i = 0; i < 3*npts; i++) 
  {
    ax[i] <<= x[i];
  }

  adouble adx, ady, adz, alength=0.;
  E_Int i1;
  for (E_Int i = 1; i < npts; i++)
  {
    i1 = i-1;
    adx = ax[3*i] - ax[3*i1];
    ady = ax[3*i+1] - ax[3*i1+1];
    adz = ax[3*i+2] - ax[3*i1+2];
    alength += sqrt(adx*adx+ady*ady+adz*adz);
  }

  alength >>= length;

  trace_off();

  //typedef double E_Float;

  //fos_forward(0, 1, 3*npts, 0, x, dx, &length, df);
  gradient(0, 3*npts, x, df);

  //printf("length=%g\n", length);

  // Export derivatives
  PyObject* tpl;
  tpl = K_ARRAY::buildArray3(3, "dFx,dFy,dFz", im, jm, km, f->getApi());
  K_FLD::FldArrayF* fo;  
  K_ARRAY::getFromArray3(tpl, fo);
  E_Float* dFxt = fo->begin(1);
  E_Float* dFyt = fo->begin(2);
  E_Float* dFzt = fo->begin(3);

  for (E_Int i = 0; i < npts; i++)
  {
    dFxt[i] = df[3*i]; // dL/dxi
    dFyt[i] = df[3*i+1]; // dL/dyi
    dFzt[i] = df[3*i+2]; // dL/dzi
  }

  delete [] ax; delete [] df;

  RELEASESHAREDS(array, f);
  RELEASESHAREDS(tpl, fo);
  
  return tpl;
}

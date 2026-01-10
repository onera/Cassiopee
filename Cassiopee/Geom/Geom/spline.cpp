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
# include "geom.h"
# include <vector>
using namespace K_FLD;
using namespace std;

//===========================================================================
/* Array spline */
//===========================================================================
PyObject* K_GEOM::spline(PyObject* self, PyObject* args)
{
  PyObject* Array;
  E_Int N,M;
  E_Int ordern, orderm;
  E_Float density;
  if (!PYPARSETUPLE_(args, O_ IIII_ R_,
                    &Array, &ordern, &N, &orderm, &M, &density))
  {
    return NULL;
  }

  if ((ordern < 1) || (orderm < 1))
  {
    PyErr_SetString(PyExc_TypeError,
                    "spline: order must be >= 1.");
    return NULL;
  }

  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(Array, varString, f, im, jm, km,
                                     cn, eltType);
  if (res != 1)
  {
    if (res == 2) RELEASESHAREDU(Array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "spline: input array not valid.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDS(Array, f);
    PyErr_SetString(PyExc_TypeError,
                    "spline: coordinates not found in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if (im < 2)
  {
    RELEASESHAREDS(Array, f);
    PyErr_SetString(PyExc_TypeError,
                    "spline: minimum 2 control points required.");
    return NULL;
  }

  if (im < ordern)
  {
    RELEASESHAREDS(Array, f);
    PyErr_SetString(PyExc_TypeError,
                    "spline: number of control points must be greater than order.");
    return NULL;
  }

  if (jm > 1 && jm < orderm)// ij-array
  {
    RELEASESHAREDS(Array, f);
    PyErr_SetString(PyExc_TypeError,
                    "spline: number of control points must be greater than order.");
    return NULL;
  }

  if (km > 1)
  {
    RELEASESHAREDS(Array, f);
    PyErr_SetString(PyExc_TypeError,
                    "spline: array must have one or two dimensions.");
    return NULL;
  }

  if (N < im || M < jm)
  {
    RELEASESHAREDS(Array, f);
    PyErr_SetString(PyExc_TypeError,
                    "spline: N and M must be greater than the number of control points.");
    return NULL;
  }
  /* Fin des tests */
  E_Int api = f->getApi();
  E_Float* xt = f->begin(posx);
  E_Float* yt = f->begin(posy);
  E_Float* zt = f->begin(posz);

  if (im != 1)
  {
    if (jm == 1)
    {
      K_FLD::FldArrayF PF;
      K_COMPGEOM::regularSpline(im, ordern, N, density, xt, yt, zt, PF);
      RELEASESHAREDS(Array, f);
      PyObject* tpl = K_ARRAY::buildArray3(PF, "x,y,z", PF.getSize(), 1, 1, api);
      return tpl;
    }
    else
    {
      K_FLD::FldArrayF PF;
      E_Int niout, njout;
      K_COMPGEOM::regularSpline2D(im, jm, ordern, N, orderm, M,
                                  density, xt, yt, zt, PF, niout, njout);
      RELEASESHAREDS(Array, f);
      PyObject* tpl = K_ARRAY::buildArray3(PF, "x,y,z", niout, njout, 1, api);
      return tpl;
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "spline: control points array must be 1D, 2D.");
    return  NULL;
  }
}


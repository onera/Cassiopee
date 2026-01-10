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

# include "transform.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// Binomial coefficients Ck^n
// may be better if return long
size_t binomialCoeff(E_Int k, E_Int n)
{
  if (k > n - k) k = n - k;  // Use symmetry
  size_t res = 1;
  for (E_Int i = 0; i < k; i++)
  {
    res *= (n - i);
    res /= (i + 1);
  }
  return res;
}
// Bernstein bases of degree i
E_Float B(E_Int p, E_Int n, E_Float u)
{
  size_t C = binomialCoeff(p, n);
  E_Float res1 = std::pow(u, p);
  E_Float res2 = std::pow(1.-u, n-p);
  return res1*res2*C;
}

// ============================================================================
/*
  set free form
  IN: array: deja avec dx,dy,dz cree
  IN: control: contient dx,dy,dz
  Rempli dx,dy,dz rempli in place.
*/
// ============================================================================
PyObject* K_TRANSFORM::_freeForm(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* control;
  if (!PYPARSETUPLE_(args, OO_, &array, &control))
    return NULL;

  // Check array
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray3(array, varString1, f1,
                                      im1, jm1, km1, cn1, eltType1);

  if (res1 != 1 && res1 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "freeForm: 1st argument is invalid.");
    return NULL;
  }

  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray3(control, varString2, f2,
                                      im2, jm2, km2, cn2, eltType2);

  if (res2 != 1 && res2 != 2)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "freeForm: 2nd argument is invalid.");
    return NULL;
  }

  if (res2 == 2)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    RELEASESHAREDB(res2, control, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "freeForm: 2nd argument must be structured.");
    return NULL;
  }

  // Presence des coordonnees ?
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    RELEASESHAREDB(res2, control, f2, cn2);
    PyErr_SetString(PyExc_ValueError,
                    "freeForm: coordinates not found in 1st argument.");
    return NULL;
  }
  posx1++; posy1++; posz1++;

  E_Int posdx1 = K_ARRAY::isNamePresent("dx", varString1);
  E_Int posdy1 = K_ARRAY::isNamePresent("dy", varString1);
  E_Int posdz1 = K_ARRAY::isNamePresent("dz", varString1);
  if (posdx1 == -1 || posdy1 == -1 || posdz1 == -1)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    RELEASESHAREDB(res2, control, f2, cn2);
    PyErr_SetString(PyExc_ValueError,
                    "freeForm: dx,dy,dz not found in 1st argument.");
    return NULL;
  }
  posdx1++; posdy1++; posdz1++;

  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    RELEASESHAREDB(res2, control, f2, cn2);
    PyErr_SetString(PyExc_ValueError,
                    "freeForm: coordinates not found in 2nd argument.");
    return NULL;
  }
  posx2++; posy2++; posz2++;

  E_Int posdx2 = K_ARRAY::isNamePresent("dx", varString2);
  E_Int posdy2 = K_ARRAY::isNamePresent("dy", varString2);
  E_Int posdz2 = K_ARRAY::isNamePresent("dz", varString2);
  if (posdx2 == -1 || posdy2 == -1 || posdz2 == -1)
  {
    RELEASESHAREDB(res1, array, f1, cn1);
    RELEASESHAREDB(res2, control, f2, cn2);
    PyErr_SetString(PyExc_ValueError,
                    "freeForm: dx,dy,dz not found in 2nd argument.");
    return NULL;
  }
  posdx2++; posdy2++; posdz2++;

  // input array
  E_Int npts = f1->getSize();

  E_Float* f1x = f1->begin(posx1);
  E_Float* f1y = f1->begin(posy1);
  E_Float* f1z = f1->begin(posz1);

  E_Float* d1x = f1->begin(posdx1);
  E_Float* d1y = f1->begin(posdy1);
  E_Float* d1z = f1->begin(posdz1);

  // control point array im2,jm2,km2,f2 (structured)
  E_Float* f2x = f2->begin(posx2);
  E_Float* f2y = f2->begin(posy2);
  E_Float* f2z = f2->begin(posz2);

  E_Float* d2x = f2->begin(posdx2);
  E_Float* d2y = f2->begin(posdy2);
  E_Float* d2z = f2->begin(posdz2);

  // Control points corners
  E_Float eps = 1.e-10;
  E_Float xo = f2x[0]-eps;
  E_Float yo = f2y[0]-eps;
  E_Float zo = f2z[0]-eps;
  E_Float xm = f2x[im2-1]+eps;
  E_Float ym = f2y[(jm2-1)*im2]+eps;
  E_Float zm = f2z[(km2-1)*im2*jm2]+eps;
  E_Float dxm = xm-xo;
  E_Float dym = ym-yo;
  E_Float dzm = zm-zo;

  //printf("%g %g\n", xo, xm);
  //printf("%g %g\n", yo, ym);
  //printf("%g %g\n", zo, zm);

  E_Int ind;
  E_Float x, y, z;
  E_Float r1x, r1y, r1z;
  E_Float r2x, r2y, r2z;
  E_Float r3x, r3y, r3z;
  E_Float u, v, w;
  E_Float bx, by, bz;
  //E_Float dx, dy, dz;

  // doit on interpoler les points ou la deformation?
  for (E_Int n = 0; n < npts; n++)
  {
    x = f1x[n]; y = f1y[n]; z = f1z[n];
    //ii = E_Int( (x - xo) / im2 );
    //jj = E_Int( (y - yo) / jm2 );
    //kk = E_Int( (z - zo) / km2 );
    //ind = ii + jj*im2 + kk*im2*jm2;

    u = (x - xo) / dxm;
    v = (y - yo) / dym;
    w = (z - zo) / dzm;

    r3x = 0.; r3y = 0.; r3z = 0.;
    for (E_Int k = 0; k < km2; k++)
    {
      bz = B(k, km2-1, w);
      r2x = 0.; r2y = 0.; r2z = 0.;

      for (E_Int j = 0; j < jm2; j++)
      {
        by = B(j, jm2-1, v);
        r1x = 0.; r1y = 0.; r1z = 0.;
        for (E_Int i = 0; i < im2; i++)
        {
          bx = B(i, im2-1, u);
          ind = i+j*im2+k*im2*jm2;
          //r1x += bx * (f2x[ind]+d2x[ind]);
          //r1y += bx * (f2y[ind]+d2y[ind]);
          //r1z += bx * (f2z[ind]+d2z[ind]);
          r1x += bx * (d2x[ind]);
          r1y += bx * (d2y[ind]);
          r1z += bx * (d2z[ind]);
        }
        r2x += by * r1x;
        r2y += by * r1y;
        r2z += by * r1z;
      }
      r3x += bz * r2x;
      r3y += bz * r2y;
      r3z += bz * r2z;
    }
    // sortie
    //dx = r3x-x; dy = r3y-y; dz = r3z-z;
    d1x[n] = r3x;
    d1y[n] = r3y;
    d1z[n] = r3z;
  }
  RELEASESHAREDB(res1, array, f1, cn1);
  RELEASESHAREDB(res2, control, f2, cn2);

  Py_INCREF(Py_None);
  return Py_None;
}

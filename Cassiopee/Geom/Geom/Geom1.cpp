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

// Analytical geometries creation

#include "geom.h"
#include <math.h>
using namespace K_FLD;
using namespace K_CONST;

// ============================================================================
/* Create a line of N points passing by P1 and P2 */
// ============================================================================
PyObject* K_GEOM::lineMesh(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float x1, y1, z1;
  E_Float x2, y2, z2;

  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ I_,
                    &x1, &y1, &z1, &x2, &y2, &z2, &N))
  {
      return NULL;
  }

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError, "line: insufficient number of point.");
    return NULL;
  }

  // Create a line
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", N, 1, 1);
  E_Float* coordx = K_ARRAY::getFieldPtr(tpl);
  E_Float* coordy = coordx + N;
  E_Float* coordz = coordy + N;

  E_Float delta = 1./(N-1.);
  E_Float dx12 = delta * (x2-x1);
  E_Float dy12 = delta * (y2-y1);
  E_Float dz12 = delta * (z2-z1);

  for (E_Int i = 0; i < N; i++)
  {
    coordx[i] = x1 + i*dx12;
    coordy[i] = y1 + i*dy12;
    coordz[i] = z1 + i*dz12;
  }

  return tpl;
}

// ============================================================================
/* Create a circle of center C and radius R,
   between tetas and tetae angles */
// ============================================================================
PyObject* K_GEOM::circleMesh(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float xc, yc, zc;
  E_Float R, tetas, tetae;

  if (!PYPARSETUPLE_(args, TRRR_ RRR_ I_,
                    &xc, &yc, &zc, &R, &tetas, &tetae, &N))
  {
      return NULL;
  }
  E_Float pi = 4*atan(1.);
  E_Float t1 = tetas*pi/180.;
  E_Float t2 = tetae*pi/180.;

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "circle: insufficient number of points.");
    return NULL;
  }
  
  // Create a portion of circle
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", N, 1, 1);
  E_Float* coordx = K_ARRAY::getFieldPtr(tpl);
  E_Float* coordy = coordx + N;
  E_Float* coordz = coordy + N;

  E_Float delta = 1./(N-1.);
  E_Float alpha;
  for (E_Int i = 0; i < N; i++)
  {
    alpha = t1+i*delta*(t2-t1);
    coordx[i] = xc + R*cos(alpha);
    coordy[i] = yc + R*sin(alpha);
    coordz[i] = zc;
  }

  return tpl;  
}

// ============================================================================
/* Create a sphere of center C and radius R */
// ============================================================================
PyObject* K_GEOM::sphereMesh(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float xc, yc, zc;
  E_Float R;
  if (!PYPARSETUPLE_(args, TRRR_ R_ I_,
                    &xc, &yc, &zc, &R, &N))
  {
      return NULL;
  }

  E_Float pi = 4*atan(1.);

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "sphere: insufficient number of points.");
    return NULL;
  }
  
  // Create a sphere
  E_Int P = 2*N;
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", N, P, 1);
  E_Float* xt = K_ARRAY::getFieldPtr(tpl);
  E_Float* yt = xt + N*P;
  E_Float* zt = yt + N*P;

  E_Float alpha, beta, cbeta, sbeta, x1, y1, z1;
  E_Int ind;

  E_Float delta = pi/(N-1.);
  E_Float deltap = (2*pi)/(P-1.);
  
  for (E_Int j = 0; j < P; j++)
    for (E_Int i = 0; i < N; i++)
    {
      alpha = i*delta;
      beta = j*deltap;
      x1 = R*cos(alpha);
      y1 = R*sin(alpha);
      z1 = 0.;
      ind = i + j*N;
      cbeta = cos(beta);
      sbeta = sin(beta);
      xt[ind] = xc + x1;
      yt[ind] = yc + cbeta*y1 - sbeta*z1;
      zt[ind] = zc + sbeta*y1 + cbeta*z1;
    }
  
  return tpl;
}

// ============================================================================
/* Create a cone of center C, basis radius Rb, vertex radius Rv, 
   and height H */
// ============================================================================
PyObject* K_GEOM::coneMesh(PyObject* self, PyObject* args)
{
  E_Int N;
  E_Float xc, yc, zc;
  E_Float Rb, Rv, H;
  if (!PYPARSETUPLE_(args, TRRR_ RRR_ I_,
                    &xc, &yc, &zc, &Rb, &Rv, &H, &N))
  {
      return NULL;
  }
  E_Float pi = 4*atan(1.);

  // Data check
  if (N < 2)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cone: insufficient number of points.");
    return NULL;
  }
  if (K_FUNC::fEqualZero(H) == true)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "cone: H must be non null.");
    return NULL;
  }
  
  // Create a cone
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", N, N, 1);

  E_Float alpha;
  E_Float x1, y1, z1;
  E_Float* xt = K_ARRAY::getFieldPtr(tpl);
  E_Float* yt = xt + N*N;
  E_Float* zt = yt + N*N;
  
  E_Float delta = 2.*pi/(N-1.);
  E_Float rapport = (Rb-Rv)/H;
  E_Float hk = H/(N-1.);
  E_Float Rk;
  z1 = 0.;
  
  for (E_Int k = 0; k < N; k++)
  {
    Rk = Rb-z1*rapport;
    for (E_Int i = 0; i < N; i++)
    {
      alpha = i*delta;
      x1 = Rk*cos(alpha);
      y1 = Rk*sin(alpha);
      xt[i+k*N] = xc + x1;
      yt[i+k*N] = yc + y1;
      zt[i+k*N] = zc + z1;
    }
    z1 = z1+hk; 
  }
  
  return tpl;  
}

// ============================================================================
/* Create a triangle */
// ============================================================================
PyObject* K_GEOM::triangleMesh(PyObject* self, PyObject* args)
{
  E_Float x1, y1, z1;
  E_Float x2, y2, z2;
  E_Float x3, y3, z3;
  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ TRRR_,
                    &x1, &y1, &z1, &x2, &y2, &z2, &x3, &y3, &z3))
  {
      return NULL;
  }
  
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z",
                                      3, 1, -1, "TRI");

  E_Float* xt = K_ARRAY::getFieldPtr(tpl);
  E_Float* yt = xt + 3;
  E_Float* zt = yt + 3;
  E_Int* cn1 = K_ARRAY::getConnectPtr(tpl);
  E_Int* cn2 = cn1 + 1;
  E_Int* cn3 = cn2 + 1;

  xt[0] = x1; yt[0] = y1; zt[0] = z1;
  xt[1] = x2; yt[1] = y2; zt[1] = z2;
  xt[2] = x3; yt[2] = y3; zt[2] = z3;
  cn1[0] = 1; cn2[0] = 2; cn3[0] = 3;
  return tpl;
}
// ============================================================================
/* Create a quadrangle */
// ============================================================================
PyObject* K_GEOM::quadrangleMesh(PyObject* self, PyObject* args)
{
  E_Float x1, y1, z1;
  E_Float x2, y2, z2;
  E_Float x3, y3, z3;
  E_Float x4, y4, z4;
  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ TRRR_ TRRR_,
                    &x1, &y1, &z1, &x2, &y2, &z2, &x3, &y3, &z3, &x4, &y4, &z4))
  {
      return NULL;
  }

  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z",
                                      4, 1, -1, "QUAD");

  E_Float* xt = K_ARRAY::getFieldPtr(tpl);
  E_Float* yt = xt + 4;
  E_Float* zt = yt + 4;
  E_Int* cn1 = K_ARRAY::getConnectPtr(tpl);
  E_Int* cn2 = cn1 + 1;
  E_Int* cn3 = cn2 + 1;
  E_Int* cn4 = cn3 + 1;

  xt[0] = x1; yt[0] = y1; zt[0] = z1;
  xt[1] = x2; yt[1] = y2; zt[1] = z2;
  xt[2] = x3; yt[2] = y3; zt[2] = z3;
  xt[3] = x4; yt[3] = y4; zt[3] = z4;
  cn1[0] = 1; cn2[0] = 2; cn3[0] = 3; cn4[0] = 4;

  return tpl;
}

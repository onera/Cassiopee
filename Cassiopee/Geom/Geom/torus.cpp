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

// Toroid surfaces creation

#include "geom.h"
#include <math.h>
using namespace K_FLD;
using namespace K_CONST;

// ============================================================================
/* Create a toroid surface of center C, radii R and r.
   - R is the distance from the center of the tube to the center of the torus,
   - r is the radius of the tube. */
// ============================================================================
PyObject* K_GEOM::torus( PyObject* self, PyObject* args )
{
  E_Int NR, Nr;
  E_Float xc, yc, zc;
  E_Float R, r;
  E_Float alphas, alphae, betas, betae;

  if (!PYPARSETUPLE_(args, TRRR_ RRRR_ RR_ II_, 
                    &xc, &yc, &zc, &R, &r,
                    &alphas, &alphae, &betas, &betae, &NR, &Nr))
  {
      return NULL;
  }
  E_Float pi = 4*atan(1.);

  alphas *= pi/180.; // deg. to rad. conversion.
  alphae *= pi/180.;
  betas *= pi/180.;
  betae *= pi/180.;

  // Data check
  if (NR < 2 || Nr < 2)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "torus: insufficient number of points.");
    return NULL;
  }
  if (r > R)
  {
    PyErr_SetString(PyExc_ValueError, 
                    "torus: degenerated surface (r must be lower than R).");
    return NULL;
  }
  
  E_Bool valid_angles = (alphas >= 0.) && (alphas <= 2*pi);
  valid_angles &= (alphae >= 0.) && (alphae <= 2*pi);
  valid_angles &= (betas >= 0.) && (betas <= 2*pi);
  valid_angles &= (betae >= 0.) && (betae <= 2*pi);
  valid_angles &= (alphas < alphae) && (betas < betae);

  if (!valid_angles)
  {
    PyErr_SetString(
      PyExc_ValueError, 
      "torus: invalid angles. They must be specified in the range [0,360] and the starting angles must be lower than the ending ones.");
    return NULL;
  }

  // Create a toroid surface
  E_Float alpha, beta, ca, sa, cb, sb;
  E_Float dalpha = (alphae - alphas)/(NR-1.);
  E_Float dbeta  = (betae - betas)/(Nr-1.);
  E_Int ind;

  FldArrayF coord(NR*Nr, 3);
  E_Float* xt = coord.begin(1);
  E_Float* yt = coord.begin(2);
  E_Float* zt = coord.begin(3);
  
  for (E_Int j = 0; j < NR; ++j)
  {
    alpha = alphas + j * dalpha;
    ca = ::cos(alpha);
    sa = ::sin(alpha);
    
    for (E_Int i = 0; i < Nr; ++i)
    {
      beta  = betae - i * dbeta; // to have the normals toward the exterior.
      cb = ::cos(beta);
      sb = ::sin(beta);
      ind = i + j * Nr;

      xt[ind] = xc + (R + r * cb) * ca;
      yt[ind] = yc + (R + r * cb) * sa;
      zt[ind] = zc + r * sb;
    }
  }
  
  // Build array
  PyObject* tpl = K_ARRAY::buildArray(coord, "x,y,z", Nr, NR, 1);
  return tpl;

  return NULL;
}


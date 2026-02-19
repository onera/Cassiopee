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
# include "metric.h"
#include <math.h>


//=============================================================================
// Surface vector and normal to the interfaces for structured grids.
// IN: im, im, km: Number of mesh vertices along %i, %j, %k
// IN: inti: Total number of interfaces along the %i direction
// IN: intij: Total number of interfaces along the %i direction and %j direction
// IN: x, y, z: Vertex coordinates
// OUT: surfx, surfy, surfz: Surface vector
// OUT: snorm: Norm of the surface vector
//=============================================================================
void K_METRIC::compIntSurf(
  const E_Int im, const E_Int jm, const E_Int km,
  const E_Int inti, const E_Int intij,
  const E_Float* x, const E_Float* y, const E_Float* z,
  E_Float* surfx, E_Float* surfy, E_Float* surfz, E_Float* snorm
)
{
  E_Int imjm = im * jm;
  E_Int im1 = im - 1;
  E_Int jm1 = jm - 1;
  E_Int km1 = km - 1;
  E_Int im1jm1 = im1 * jm1;

  // PARAMETRES D'INCREMENTATION POUR LES INTERFACES "I"
  E_Int incnodeij  = im;
  E_Int incnodeik  = imjm;
  E_Int incnodeijk = im + imjm;

  // PARAMETRES D'INCREMENTATION POUR LES INTERFACES "J"
  E_Int incnodejk  = imjm;
  E_Int incnodeji  = 1;
  E_Int incnodejik = 1 + imjm;

  // PARAMETRES D'INCREMENTATION POUR LES INTERFACES "K"
  E_Int incnodeki  = 1;
  E_Int incnodekj  = im;
  E_Int incnodekij = 1 + im;

  E_Float sgn = K_CONST::ONE;

  #pragma omp parallel
  {
    E_Int i, j, k, rem;
    E_Int li, lj, lk, l1, l2, l3, l4;
    E_Float d13x, d13y, d13z, d24x, d24y, d24z;
    E_Float nx, ny, nz;

    // I-direction faces
    #pragma omp for nowait
    for (E_Int idx = 0; idx < im * jm1 * km1; idx++)
    {
      k = idx / (jm1 * im);
      rem = idx % (jm1 * im);
      j = rem / im;
      i = rem % im;

      li = i + j * im + k * im * jm1;
      l1 = i + j * im + k * imjm;
      l2 = l1 + incnodeij;
      l3 = l1 + incnodeijk;
      l4 = l1 + incnodeik;

      d13x = x[l3] - x[l1];
      d13y = y[l3] - y[l1];
      d13z = z[l3] - z[l1];

      d24x = x[l4] - x[l2];
      d24y = y[l4] - y[l2];
      d24z = z[l4] - z[l2];

      nx = K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
      ny = K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
      nz = K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

      surfx[li] = sgn * nx;
      surfy[li] = sgn * ny;
      surfz[li] = sgn * nz;
      snorm[li] = sqrt(nx * nx + ny * ny + nz * nz);
    }

    // J-direction faces
    #pragma omp for nowait
    for (E_Int idx = 0; idx < km1 * jm * im1; idx++)
    {
      k = idx / (jm * im1);
      rem = idx % (jm * im1);
      j = rem / im1;
      i = rem % im1;

      lj = i + j * im1 + k * im1 * jm + inti;
      l1 = i + j * im + k * imjm;
      l2 = l1 + incnodejk;
      l3 = l1 + incnodejik;
      l4 = l1 + incnodeji;

      d13x = x[l3] - x[l1];
      d13y = y[l3] - y[l1];
      d13z = z[l3] - z[l1];

      d24x = x[l4] - x[l2];
      d24y = y[l4] - y[l2];
      d24z = z[l4] - z[l2];

      nx = K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
      ny = K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
      nz = K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

      surfx[lj] = sgn * nx;
      surfy[lj] = sgn * ny;
      surfz[lj] = sgn * nz;
      snorm[lj] = sqrt(nx * nx + ny * ny + nz * nz);
    }

    // K-direction faces
    #pragma omp for
    for (E_Int idx = 0; idx < km * jm1 * im1; idx++)
    {
      k = idx / (jm1 * im1);
      rem = idx % (jm1 * im1);
      j = rem / im1;
      i = rem % im1;

      lk = i + j * im1 + k * im1jm1 + intij;
      l1 = i + j * im + k * imjm;
      l2 = l1 + incnodeki;
      l3 = l1 + incnodekij;
      l4 = l1 + incnodekj;

      d13x = x[l3] - x[l1];
      d13y = y[l3] - y[l1];
      d13z = z[l3] - z[l1];

      d24x = x[l4] - x[l2];
      d24y = y[l4] - y[l2];
      d24z = z[l4] - z[l2];

      nx = K_CONST::ONE_HALF * (d13y * d24z - d13z * d24y);
      ny = K_CONST::ONE_HALF * (d13z * d24x - d13x * d24z);
      nz = K_CONST::ONE_HALF * (d13x * d24y - d13y * d24x);

      surfx[lk] = sgn * nx;
      surfy[lk] = sgn * ny;
      surfz[lk] = sgn * nz;
      snorm[lk] = sqrt(nx * nx + ny * ny + nz * nz);
    }
  }
}

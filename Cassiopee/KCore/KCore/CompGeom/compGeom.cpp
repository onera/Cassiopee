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

#include "Array/Array.h"
#include "CompGeom/compGeom.h"
#include <vector>
#include <math.h>

using namespace K_FLD;
using namespace K_FUNC;
using namespace std;

// ============================================================================
// Evaluates the orientation of two structured blocks. 
// Returns -1 if the normals are reversed, 1 otherwise.
// ============================================================================
void K_COMPGEOM::rectifyNormals(const E_Int ni1, const E_Int nj1, const E_Int ind1,
  const E_Float* x1, const E_Float* y1, const E_Float* z1,
  const E_Int ni2, const E_Int nj2, const E_Int ind2,
  const E_Float* x2, const E_Float* y2, const E_Float* z2,
  const E_Float distmin,
  E_Int& isopp)
{
  E_Int i1, j1, i2, j2;
  E_Int im1, im2, jm1, jm2, ip1, ip2, jp1, jp2;
  E_Int ind1p, ind2p;
  E_Int ind11, ind12, ind13, ind14;
  E_Int ind21, ind22, ind23, ind24;
  E_Float eps, ps1, ps2, ps3, psmin, matchTol;
  E_Float x13, y13, z13, x24, y24, z24;
  E_Float x11p, y11p, z11p, x22p, y22p, z22p;
  E_Float n1x, n1y, n1z, n2x, n2y, n2z;
  E_Float v1x, v1y, v1z, v2x, v2y, v2z;
  E_Float norm1, norm2, norm, norm1p, norm2p, normp;
  E_Int isclassical;

  isopp = 1;
  isclassical = 0;

  eps = 1.e-12;
  matchTol = 1.e-12;
  psmin = 0.1; // minimum value to avoid considering blocks orthogonal

  j1 = ind1 / ni1 + 1;
  i1 = ind1 - (j1 - 1) * ni1 + 1;
  j2 = ind2 / ni2 + 1;
  i2 = ind2 - (j2 - 1) * ni2 + 1;

  // A1B1C1D1 indices
  im1 = i1 - 1;
  ip1 = i1;
  jm1 = j1 - 1;
  jp1 = j1;
  if (i1 == ni1) {im1 = ni1 - 2; ip1 = ni1 - 1;}
  if (j1 == nj1) {jm1 = nj1 - 2; jp1 = nj1 - 1;}

  ind11 = im1 + jm1 * ni1; // A1 (i,j)
  ind12 = ip1 + jm1 * ni1; // B1 (i+1,j)
  ind13 = ip1 + jp1 * ni1; // C1 (i+1,j+1)
  ind14 = im1 + jp1 * ni1; // D1 (i,j+1)

  // A1C1 
  x13 = x1[ind13] - x1[ind11];
  y13 = y1[ind13] - y1[ind11];
  z13 = z1[ind13] - z1[ind11];

  // B1D1 
  x24 = x1[ind14] - x1[ind12];
  y24 = y1[ind14] - y1[ind12];
  z24 = z1[ind14] - z1[ind12];

  // A1C1 x B1D1 -> n1
  n1x = y13 * z24 - y24 * z13;
  n1y = x24 * z13 - x13 * z24;
  n1z = x13 * y24 - x24 * y13;
  norm1 = n1x * n1x + n1y * n1y + n1z * n1z;

  // A2B2C2D2 indices
  im2 = i2 - 1;
  ip2 = i2;
  jm2 = j2 - 1;
  jp2 = j2;
  if (i2 == ni2) {im2 = ni2 - 2; ip2 = ni2 - 1;}
  if (j2 == nj2) {jm2 = nj2 - 2; jp2 = nj2 - 1;}

  ind21 = im2 + jm2 * ni2; // A2 (i,j)
  ind22 = ip2 + jm2 * ni2; // B2 (i+1,j)
  ind23 = ip2 + jp2 * ni2; // C2 (i+1,j+1)
  ind24 = im2 + jp2 * ni2; // D2 (i,j+1)

  // A2C2 
  x13 = x2[ind23] - x2[ind21];
  y13 = y2[ind23] - y2[ind21];
  z13 = z2[ind23] - z2[ind21];

  // B2D2 
  x24 = x2[ind24] - x2[ind22];
  y24 = y2[ind24] - y2[ind22];
  z24 = z2[ind24] - z2[ind22];

  // A2C2 x B2D2 -> n2  
  n2x = y13 * z24 - y24 * z13;
  n2y = x24 * z13 - x13 * z24;
  n2z = x13 * y24 - x24 * y13;
  norm2 = n2x * n2x + n2y * n2y + n2z * n2z;

  // dot product normalized (ps1): compare cell normal orientation
  norm = 1. / sqrt(norm1 * norm2);
  ps1 = (n1x * n2x + n1y * n2y + n1z * n2z) * norm;

  if (distmin > matchTol) // if far enough, classical test is sufficient
  {
    if (ps1 < -eps) isopp = -1;
    return;
  }

  ind1p = -1;
  ind2p = -1;

  // classify ind1: edge vs interior/corner
  if (i1 == 1 && j1 != 1 && j1 != nj1) ind1p = ind1 + 1;
  else if (i1 == ni1 && j1 != 1 && j1 != nj1) ind1p = ind1 - 1;
  else if (j1 == 1 && i1 != 1 && i1 != ni1) ind1p = ind1 + ni1;
  else if (j1 == nj1 && i1 != 1 && i1 != ni1) ind1p = ind1 - ni1;
  else isclassical = 1; // interior or corner point

  // classify ind2: edge vs interior/corner
  if (i2 == 1 && j2 != 1 && j2 != nj2) ind2p = ind2 + 1;
  else if (i2 == ni2 && j2 != 1 && j2 != nj2) ind2p = ind2 - 1;
  else if (j2 == 1 && i2 != 1 && i2 != ni2) ind2p = ind2 + ni2;
  else if (j2 == nj2 && i2 != 1 && i2 != ni2) ind2p = ind2 - ni2;
  else isclassical = 1; // interior or corner point 

  if (isclassical == 1)
  {
    if (ps1 < -eps) isopp = -1;
    return;
  }

  // edge-specific refined test

  // A1A1p
  x11p = x1[ind1p] - x1[ind1];
  y11p = y1[ind1p] - y1[ind1];
  z11p = z1[ind1p] - z1[ind1];
  norm1p = x11p * x11p + y11p * y11p + z11p * z11p;

  // A1A1p x n1 -> v1
  v1x = y11p * n1z - z11p * n1y;
  v1y = z11p * n1x - x11p * n1z;
  v1z = x11p * n1y - y11p * n1x;

  // A2A2p
  x22p = x2[ind2p] - x2[ind2];
  y22p = y2[ind2p] - y2[ind2];
  z22p = z2[ind2p] - z2[ind2];
  norm2p = x22p * x22p + y22p * y22p + z22p * z22p;

  // A2A2p x n2 -> v2
  v2x = y22p * n2z - z22p * n2y;
  v2y = z22p * n2x - x22p * n2z;
  v2z = x22p * n2y - y22p * n2x;

  // dot product normalized (ps2): compare edge tangent alignment
  normp = 1. / sqrt(norm1p * norm2p);
  ps2 = (x11p * x22p + y11p * y22p + z11p * z22p) * normp;

  // dot product normalized (ps3): compare edge normal orientation
  ps3 = v1x * v2x + v1y * v2y + v1z * v2z;
  ps3 = ps3 * norm * normp;

  if ((ps1 > eps && ps2 > eps) || (ps1 < eps && ps2 < eps))
  {
    if (E_abs(ps3) > psmin)
    {
      if (ps3 < -eps) isopp = -1;
    }
    else
    {
      if (ps1 < -eps) isopp = -1;
    }
  }
  else
  {
    if (ps1 < -eps) isopp = -1;
  }

  return;
}

// ============================================================================
// Compute the minimum squared distance between 2 blocks and return
// the corresponding cell indices
// ============================================================================
void K_COMPGEOM::compMeanDist(const E_Int ni1, const E_Int nj1,
  const E_Float* x1, const E_Float* y1, const E_Float* z1,
  const E_Int ni2, const E_Int nj2,
  const E_Float* x2, const E_Float* y2, const E_Float* z2,
  E_Int& ind1s, E_Int& ind2s, E_Float& dmin)
{
  E_Int ind1, ind2;
  E_Float dist;
  E_Float x10, y10, z10;
  E_Float x20, y20, z20;
  E_Int ind11, ind12, ind13, ind14, ind21, ind22, ind23, ind24;

  ind11 = 0; //imin,jmin
  ind12 = (ni1-1); //imax,jmin
  ind13 = (nj1-1)*ni1; //imin,jmax
  ind14 = (ni1-1) + (nj1-1)*ni1; //imax,jmax

  ind21 = 0; //imin,jmin
  ind22 = (ni2-1); //imax,jmin
  ind23 = (nj2-1)*ni2; //imin,jmax
  ind24 = (ni2-1) + (nj2-1)*ni2; //imax,jmax

  E_Int cornerIndices1[4] = { ind11, ind12, ind13, ind14 };
  E_Int cornerIndices2[4] = { ind21, ind22, ind23, ind24 };

  dmin = K_CONST::E_MAX_FLOAT;

  // distance between non-corner points
  for (E_Int j1 = 0; j1 < nj1; j1++)
  {
    for (E_Int i1 = 0; i1 < ni1; i1++)
    {
      ind1 = i1 + j1 * ni1;

      if (ind1 == ind11 || ind1 == ind12 || 
          ind1 == ind13 || ind1 == ind14) continue; // skip corners

      x10 = x1[ind1];
      y10 = y1[ind1];
      z10 = z1[ind1];

      for (E_Int j2 = 0; j2 < nj2; j2++)
      {
        for (E_Int i2 = 0; i2 < ni2; i2++)
        {
          ind2 = i2 + j2 * ni2;

          if (ind2 == ind21 || ind2 == ind22 ||
              ind2 == ind23 || ind2 == ind24) continue; // skip corners

          x20 = x2[ind2];
          y20 = y2[ind2];
          z20 = z2[ind2];

          dist = (x20-x10)*(x20-x10) + 
                 (y20-y10)*(y20-y10) + 
                 (z20-z10)*(z20-z10);

          if (dist < dmin)
          {
            dmin = dist;
            ind1s = ind1;
            ind2s = ind2;

            // zero distance
            if (dist < 1.e-12) return;
          }
        }
      }
    }
  }

  // distance between points of block 1 and corners of block 2
  for (E_Int j1 = 0; j1 < nj1; j1++)
  {
    for (E_Int i1 = 0; i1 < ni1; i1++)
    {
      ind1 = i1 + j1 * ni1;

      if (ind1 == ind11 || ind1 == ind12 || 
          ind1 == ind13 || ind1 == ind14) continue; // skip corners

      x10 = x1[ind1];
      y10 = y1[ind1];
      z10 = z1[ind1];

      for (E_Int k = 0; k < 4; k++)
      {
        ind2 = cornerIndices2[k];

        x20 = x2[ind2];
        y20 = y2[ind2];
        z20 = z2[ind2];

        dist = (x20-x10)*(x20-x10) + 
               (y20-y10)*(y20-y10) + 
               (z20-z10)*(z20-z10);

        if (dist < dmin)
        {
          dmin = dist;
          ind1s = ind1;
          ind2s = ind2;

          // zero distance
          if (dist < 1.e-12) return;
        }
      }
    }
  }

  // distance between points of block 2 and corners of block 1
  for (E_Int j2 = 0; j2 < nj2; j2++)
  {
    for (E_Int i2 = 0; i2 < ni2; i2++)
    {
      ind2 = i2 + j2 * ni2;

      if (ind2 == ind21 || ind2 == ind22 ||
        ind2 == ind23 || ind2 == ind24) continue; // skip corners

      x20 = x2[ind2];
      y20 = y2[ind2];
      z20 = z2[ind2];

      for (E_Int k = 0; k < 4; k++)
      {
        ind1 = cornerIndices1[k];

        x10 = x1[ind1];
        y10 = y1[ind1];
        z10 = z1[ind1];

        dist = (x20-x10)*(x20-x10) + 
               (y20-y10)*(y20-y10) + 
               (z20-z10)*(z20-z10);

        if (dist < dmin)
        {
          dmin = dist;
          ind1s = ind1;
          ind2s = ind2;

          // zero distance
          if (dist < 1.e-12) return;
        }
      }
    }
  }
  
  return;
}

//===========================================================================
// Calcul de l'aire d un triangle ABC a partir des longueurs de ses 3 cotes
// par la formule de Heron
//===========================================================================
E_Float K_COMPGEOM::compTriangleArea(
  const E_Float a, const E_Float b, const E_Float c
)
{
  E_Float ps2 = (a + b + c)*0.5;
  return sqrt(ps2*(ps2 - a)*(ps2 - b)*(ps2 - c));
}

//===========================================================================
// Calcul de la bounding box d'un array structure
//===========================================================================
void K_COMPGEOM::boundingBoxUnstruct(
  const E_Int npts, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& xmin, E_Float& ymin, E_Float& zmin,
  E_Float& xmax, E_Float& ymax, E_Float& zmax
)
{
  xmin =  K_CONST::E_MAX_FLOAT;
  ymin =  K_CONST::E_MAX_FLOAT;
  zmin =  K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int nthreads = __NUMTHREADS__;
  E_Float* xminl = new E_Float [nthreads];
  E_Float* yminl = new E_Float [nthreads];
  E_Float* zminl = new E_Float [nthreads];
  E_Float* xmaxl = new E_Float [nthreads];
  E_Float* ymaxl = new E_Float [nthreads];
  E_Float* zmaxl = new E_Float [nthreads];

  #pragma omp parallel
  {
    E_Int ithread = __CURRENT_THREAD__;

    xminl[ithread] =  K_CONST::E_MAX_FLOAT;
    yminl[ithread] =  K_CONST::E_MAX_FLOAT;
    zminl[ithread] =  K_CONST::E_MAX_FLOAT;
    xmaxl[ithread] = -K_CONST::E_MAX_FLOAT;
    ymaxl[ithread] = -K_CONST::E_MAX_FLOAT;
    zmaxl[ithread] = -K_CONST::E_MAX_FLOAT;

    #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      xminl[ithread] = K_FUNC::E_min(xminl[ithread], xt[ind]);
      yminl[ithread] = K_FUNC::E_min(yminl[ithread], yt[ind]);
      zminl[ithread] = K_FUNC::E_min(zminl[ithread], zt[ind]);
      xmaxl[ithread] = K_FUNC::E_max(xmaxl[ithread], xt[ind]);
      ymaxl[ithread] = K_FUNC::E_max(ymaxl[ithread], yt[ind]);
      zmaxl[ithread] = K_FUNC::E_max(zmaxl[ithread], zt[ind]);
    }
  }

  //final reduction
  for (E_Int ithread = 0; ithread < nthreads; ithread++) 
  {
    xmax = K_FUNC::E_max(xmaxl[ithread], xmax);
    ymax = K_FUNC::E_max(ymaxl[ithread], ymax);
    zmax = K_FUNC::E_max(zmaxl[ithread], zmax);
    xmin = K_FUNC::E_min(xminl[ithread], xmin);
    ymin = K_FUNC::E_min(yminl[ithread], ymin);
    zmin = K_FUNC::E_min(zminl[ithread], zmin);
  }

  //clean
  delete [] xminl;
  delete [] xmaxl;
  delete [] yminl;
  delete [] ymaxl;
  delete [] zminl;
  delete [] zmaxl;
}

//=============================================================================
// Calcul de la bounding box d'un array structure
//=============================================================================
void K_COMPGEOM::boundingBoxStruct(E_Int im, E_Int jm, E_Int km, 
                                   E_Float* x, E_Float* y, E_Float* z,
                                   E_Float& xmin, E_Float& ymin, E_Float& zmin,
                                   E_Float& xmax, E_Float& ymax, E_Float& zmax)
{

  E_Int im1 = K_FUNC::E_max(1, im-1);
  E_Int jm1 = K_FUNC::E_max(1, jm-1);
  E_Int km1 = K_FUNC::E_max(1, km-1);
  E_Int im1jm1 = im1*jm1;

  xmin =  K_CONST::E_MAX_FLOAT;
  ymin =  K_CONST::E_MAX_FLOAT;
  zmin =  K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int nthreads = __NUMTHREADS__;
  E_Float* xminl = new E_Float [nthreads];
  E_Float* yminl = new E_Float [nthreads];
  E_Float* zminl = new E_Float [nthreads];
  E_Float* xmaxl = new E_Float [nthreads];
  E_Float* ymaxl = new E_Float [nthreads];
  E_Float* zmaxl = new E_Float [nthreads];

  #pragma omp parallel
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Int i, j, k;
    E_Int ind;

    E_Float xmincell, ymincell, zmincell, xmaxcell, ymaxcell, zmaxcell;

    xminl[ithread] =  K_CONST::E_MAX_FLOAT;
    yminl[ithread] =  K_CONST::E_MAX_FLOAT;
    zminl[ithread] =  K_CONST::E_MAX_FLOAT;
    xmaxl[ithread] = -K_CONST::E_MAX_FLOAT;
    ymaxl[ithread] = -K_CONST::E_MAX_FLOAT;
    zmaxl[ithread] = -K_CONST::E_MAX_FLOAT;

    //faces imin & imax
    #pragma omp for collapse(2)
    for (E_Int k = 0; k < km1; k++)
    {
      for (E_Int j = 0; j < jm1; j++)
      {
        // imin --------------------------------------------------
        i = 0;
        ind = i+j*im1+k*im1jm1;

        boundingBoxOfStructCell(
          ind, im, jm, km,
          x, y, z,
          xmincell, ymincell, zmincell,
          xmaxcell, ymaxcell, zmaxcell, 1
        );

        xmaxl[ithread] = K_FUNC::E_max(xmaxl[ithread], xmaxcell);
        ymaxl[ithread] = K_FUNC::E_max(ymaxl[ithread], ymaxcell);
        zmaxl[ithread] = K_FUNC::E_max(zmaxl[ithread], zmaxcell);
        xminl[ithread] = K_FUNC::E_min(xminl[ithread], xmincell);
        yminl[ithread] = K_FUNC::E_min(yminl[ithread], ymincell);
        zminl[ithread] = K_FUNC::E_min(zminl[ithread], zmincell);

        // imax --------------------------------------------------
        i = im1-1;
        ind = i+j*im1+k*im1jm1;

        boundingBoxOfStructCell(
          ind, im, jm, km,
          x, y, z,
          xmincell, ymincell, zmincell,
          xmaxcell, ymaxcell, zmaxcell, 1
        );

        xmaxl[ithread] = K_FUNC::E_max(xmaxl[ithread], xmaxcell);
        ymaxl[ithread] = K_FUNC::E_max(ymaxl[ithread], ymaxcell);
        zmaxl[ithread] = K_FUNC::E_max(zmaxl[ithread], zmaxcell);
        xminl[ithread] = K_FUNC::E_min(xminl[ithread], xmincell);
        yminl[ithread] = K_FUNC::E_min(yminl[ithread], ymincell);
        zminl[ithread] = K_FUNC::E_min(zminl[ithread], zmincell);
      }
    }

    //faces jmin & jmax
    #pragma omp for collapse(2)
    for (E_Int k = 0; k < km1; k++)
    {
      for (E_Int i = 0; i < im1; i++)
      {
        // jmin --------------------------------------------------
        j = 0;
        ind = i+j*im1+k*im1jm1;

        boundingBoxOfStructCell(
          ind, im, jm, km,
          x, y, z,
          xmincell, ymincell, zmincell,
          xmaxcell, ymaxcell, zmaxcell, 1
        );

        xmaxl[ithread] = K_FUNC::E_max(xmaxl[ithread], xmaxcell);
        ymaxl[ithread] = K_FUNC::E_max(ymaxl[ithread], ymaxcell);
        zmaxl[ithread] = K_FUNC::E_max(zmaxl[ithread], zmaxcell);
        xminl[ithread] = K_FUNC::E_min(xminl[ithread], xmincell);
        yminl[ithread] = K_FUNC::E_min(yminl[ithread], ymincell);
        zminl[ithread] = K_FUNC::E_min(zminl[ithread], zmincell);

        // jmax --------------------------------------------------
        j = jm1-1;
        ind = i+j*im1+k*im1jm1;

        boundingBoxOfStructCell(
          ind, im, jm, km,
          x, y, z,
          xmincell, ymincell, zmincell,
          xmaxcell, ymaxcell, zmaxcell, 1
        );

        xmaxl[ithread] = K_FUNC::E_max(xmaxl[ithread], xmaxcell);
        ymaxl[ithread] = K_FUNC::E_max(ymaxl[ithread], ymaxcell);
        zmaxl[ithread] = K_FUNC::E_max(zmaxl[ithread], zmaxcell);
        xminl[ithread] = K_FUNC::E_min(xminl[ithread], xmincell);
        yminl[ithread] = K_FUNC::E_min(yminl[ithread], ymincell);
        zminl[ithread] = K_FUNC::E_min(zminl[ithread], zmincell);
      }
    }

    //faces kmin & kmax
    #pragma omp for collapse(2)
    for (E_Int j = 0; j < jm1; j++)
    {
      for (E_Int i = 0; i < im1; i++)
      {
        // kmin --------------------------------------------------
        k = 0;
        ind = i+j*im1+k*im1jm1;

        boundingBoxOfStructCell(ind, im, jm, km,
          x, y, z,
          xmincell, ymincell, zmincell,
          xmaxcell, ymaxcell, zmaxcell, 1
        );

        xmaxl[ithread] = K_FUNC::E_max(xmaxl[ithread], xmaxcell);
        ymaxl[ithread] = K_FUNC::E_max(ymaxl[ithread], ymaxcell);
        zmaxl[ithread] = K_FUNC::E_max(zmaxl[ithread], zmaxcell);
        xminl[ithread] = K_FUNC::E_min(xminl[ithread], xmincell);
        yminl[ithread] = K_FUNC::E_min(yminl[ithread], ymincell);
        zminl[ithread] = K_FUNC::E_min(zminl[ithread], zmincell);

        // kmax --------------------------------------------------
        k = km1-1;
        ind = i+j*im1+k*im1jm1;

        boundingBoxOfStructCell(
          ind, im, jm, km,
          x, y, z,
          xmincell, ymincell, zmincell,
          xmaxcell, ymaxcell, zmaxcell, 1
        );

        xmaxl[ithread] = K_FUNC::E_max(xmaxl[ithread], xmaxcell);
        ymaxl[ithread] = K_FUNC::E_max(ymaxl[ithread], ymaxcell);
        zmaxl[ithread] = K_FUNC::E_max(zmaxl[ithread], zmaxcell);
        xminl[ithread] = K_FUNC::E_min(xminl[ithread], xmincell);
        yminl[ithread] = K_FUNC::E_min(yminl[ithread], ymincell);
        zminl[ithread] = K_FUNC::E_min(zminl[ithread], zmincell);
      }
    }
  }

  //final reduction
  for (E_Int ithread = 0; ithread < nthreads; ithread++) 
  {
    xmax = K_FUNC::E_max(xmaxl[ithread], xmax);
    ymax = K_FUNC::E_max(ymaxl[ithread], ymax);
    zmax = K_FUNC::E_max(zmaxl[ithread], zmax);
    xmin = K_FUNC::E_min(xminl[ithread], xmin);
    ymin = K_FUNC::E_min(yminl[ithread], ymin);
    zmin = K_FUNC::E_min(zminl[ithread], zmin);
  }

  //clean
  delete [] xminl;
  delete [] xmaxl;
  delete [] yminl;
  delete [] ymaxl;
  delete [] zminl;
  delete [] zmaxl;
}

//===========================================================================
/* Calcul de la bounding box d'un ensemble de grilles non structurees
*/ 
//===========================================================================
void K_COMPGEOM::globalBoundingBox(
  vector<E_Int>& posxt, vector<E_Int>& posyt, vector<E_Int>& poszt,
  vector<FldArrayF*>& listOfFields,
  E_Float& xmin, E_Float& ymin, E_Float& zmin,
  E_Float& xmax, E_Float& ymax, E_Float& zmax)
{
  xmin = +K_CONST::E_MAX_FLOAT;
  ymin = +K_CONST::E_MAX_FLOAT;
  zmin = +K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int nd = listOfFields.size();
  E_Float xminl, yminl, zminl, xmaxl, ymaxl, zmaxl;

  for (E_Int i = 0; i < nd; i++)
  {
    FldArrayF* field = listOfFields[i];
    E_Int size = field->getSize();
    E_Int posx = posxt[i]; 
    E_Int posy = posyt[i];
    E_Int posz = poszt[i];
    
    boundingBoxUnstruct(
      size, field->begin(posx), field->begin(posy), field->begin(posz),
      xminl, yminl, zminl,
      xmaxl, ymaxl, zmaxl
    );

    xmin = K_FUNC::E_min(xmin, xminl);
    ymin = K_FUNC::E_min(ymin, yminl);
    zmin = K_FUNC::E_min(zmin, zminl);
    xmax = K_FUNC::E_max(xmax, xmaxl);
    ymax = K_FUNC::E_max(ymax, ymaxl);
    zmax = K_FUNC::E_max(zmax, zmaxl);
  }
}

//==========================================================================
/* Bounding box de toutes les cellules d'une grille structuree
   Stockage de l'information aux centres des cellules
   IN: im, jm, km: dimensions de l'array definissant la grille
   IN: coord: coordonnees de la grille
   OUT: bbox(ncells, 6): xmin, ymin, zmin, xmax, ymax, zmax
*/
//==========================================================================
void K_COMPGEOM::boundingBoxOfStructCells(E_Int im, E_Int jm, E_Int km, 
                                          E_Float* x, E_Float* y, E_Float* z,
                                          K_FLD::FldArrayF& bbox)
{

  E_Int im1 = K_FUNC::E_max(1, im-1);
  E_Int jm1 = K_FUNC::E_max(1, jm-1);
  E_Int km1 = K_FUNC::E_max(1, km-1);

  if (bbox.getSize() == 0) bbox.malloc(im1*jm1*km1, 6);

  #pragma omp parallel
  {
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    #pragma omp for
    for (E_Int ind = 0; ind < im1*jm1*km1; ind++)
    {
      boundingBoxOfStructCell(
        ind, im, jm, km,
        x, y, z,
        xmin, ymin, zmin,
        xmax, ymax, zmax, 1
      );

      bbox(ind,1) = xmin;
      bbox(ind,2) = ymin;
      bbox(ind,3) = zmin;
      bbox(ind,4) = xmax;
      bbox(ind,5) = ymax;
      bbox(ind,6) = zmax;
    }
  }
}
//======================================================================
/* Bounding box de toutes les cellules d'une grille non structuree
   Stockage de l'information aux centres des cellules
   IN: cn: connectivite de la grille
   IN: coord: coordonnees de la grille
   OUT: bbox(nelts, 6): xmin, ymin, zmin, xmax, ymax, zmax
   bbox est alloue ici. */
//======================================================================
void K_COMPGEOM::boundingBoxOfUnstrCells(K_FLD::FldArrayI& cn,
                                         E_Float* xt, E_Float* yt, E_Float* zt,
                                         K_FLD::FldArrayF& bbox)
{
  E_Int nc = cn.getNConnect();
  E_Int ntotElts = 0;
  std::vector<E_Int> nepc(nc+1);
  nepc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
    ntotElts += nelts;
  }

  if (bbox.getSize() == 0) bbox.malloc(ntotElts, 6);

  #pragma omp parallel
  {
    E_Int nelts, elOffset;
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;

    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      elOffset = nepc[ic];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        boundingBoxOfUnstrCell(i, cm, xt, yt, zt,
                              xmin, ymin, zmin, xmax, ymax, zmax);

        bbox(i+elOffset,1) = xmin;
        bbox(i+elOffset,2) = ymin;
        bbox(i+elOffset,3) = zmin;
        bbox(i+elOffset,4) = xmax;
        bbox(i+elOffset,5) = ymax;
        bbox(i+elOffset,6) = zmax;
      }
    }
  }
}

//======================================================================
/* Bounding box de toutes les cellules d'une grille NGon
   Stockage de l'information aux centres des cellules
   IN: connect: connectivite de la grille
   IN: coord: coordonnees de la grille
   OUT: bbox(nelts, 6): xmin, ymin, zmin, xmax, ymax, zmax
   bbox est alloue ici. */
//======================================================================
void K_COMPGEOM::boundingBoxOfNGonCells(K_FLD::FldArrayI& connect,
  E_Float* xt, E_Float* yt, E_Float* zt,
  K_FLD::FldArrayF& bbox)
{
  E_Int nelts = connect.getNElts();
  if (bbox.getSize() == 0) bbox.malloc(nelts, 6);

  #pragma omp parallel default(shared)
  {
    E_Float xmin, ymin, zmin, xmax, ymax, zmax;
    #pragma omp for
    for (E_Int et = 0; et < nelts; et++)
    {
      boundingBoxOfNGonCell(et, connect, xt, yt, zt,
        xmin, ymin, zmin, xmax, ymax, zmax);

      bbox(et,1) = xmin;
      bbox(et,2) = ymin;
      bbox(et,3) = zmin;
      bbox(et,4) = xmax;
      bbox(et,5) = ymax;
      bbox(et,6) = zmax;
    }
  }
}

//=============================================================================
/* Calcul de la bounding box d'une cellule structuree */
//=============================================================================
void K_COMPGEOM::boundingBoxOfStructCell(E_Int ind, E_Int im, E_Int jm, E_Int km,
  E_Float* x, E_Float* y, E_Float* z,
  E_Float& xmin, E_Float& ymin, E_Float& zmin,
  E_Float& xmax, E_Float& ymax, E_Float& zmax,
  E_Int loc) 
{

  xmin =  K_CONST::E_MAX_FLOAT;
  ymin =  K_CONST::E_MAX_FLOAT;
  zmin =  K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int im1 = K_FUNC::E_max(1, im-1);
  E_Int jm1 = K_FUNC::E_max(1, jm-1);
  E_Int km1 = K_FUNC::E_max(1, km-1);
  E_Int imjm = im*jm;
  E_Int im1jm1 = im1*jm1;

  E_Int stepi = 1;
  E_Int stepj = im;
  E_Int stepk = imjm;

  E_Int i,j,k;

  if (loc == 1) // ind is the indice of a cell center
  {  
    k = ind / im1jm1;
    j = (ind - k*im1jm1) / im1;
    i = ind - j*im1 - k*im1jm1;

    if (im == 1) stepi = 0;
    if (jm == 1) stepj = 0;
    if (km == 1) stepk = 0;  
  }
  else // ind is the indice of a cell vertex
  {
    k = ind / imjm;
    j = (ind - k*imjm) / im;
    i = ind - j*im - k*imjm;

    if (i == im1) stepi *= -1;
    if (j == jm1) stepj *= -1;
    if (k == km1) stepk *= -1;
  
    if (im == 1) stepi = 0;
    if (jm == 1) stepj = 0;
    if (km == 1) stepk = 0;
  }

  E_Int ind0 = i + j*im + k*imjm;
  E_Int ind1 = ind0 + stepi;
  E_Int ind2 = ind0 + stepj;
  E_Int ind3 = ind2 + stepi;
  E_Int ind4 = ind0 + stepk;
  E_Int ind5 = ind1 + stepk;
  E_Int ind6 = ind2 + stepk;
  E_Int ind7 = ind3 + stepk;

  // vertex (i,j,k)
  xmin = K_FUNC::E_min(xmin, x[ind0]);
  xmax = K_FUNC::E_max(xmax, x[ind0]);
  ymin = K_FUNC::E_min(ymin, y[ind0]);
  ymax = K_FUNC::E_max(ymax, y[ind0]);
  zmin = K_FUNC::E_min(zmin, z[ind0]);
  zmax = K_FUNC::E_max(zmax, z[ind0]);

  // vertex (i+1,j,k)
  xmin = K_FUNC::E_min(xmin, x[ind1]);
  xmax = K_FUNC::E_max(xmax, x[ind1]);
  ymin = K_FUNC::E_min(ymin, y[ind1]);
  ymax = K_FUNC::E_max(ymax, y[ind1]);
  zmin = K_FUNC::E_min(zmin, z[ind1]);
  zmax = K_FUNC::E_max(zmax, z[ind1]);

  // vertex (i+1,j+1,k)
  xmin = K_FUNC::E_min(xmin, x[ind2]);
  xmax = K_FUNC::E_max(xmax, x[ind2]);
  ymin = K_FUNC::E_min(ymin, y[ind2]);
  ymax = K_FUNC::E_max(ymax, y[ind2]);
  zmin = K_FUNC::E_min(zmin, z[ind2]);
  zmax = K_FUNC::E_max(zmax, z[ind2]);

  // vertex (i,j+1,k)
  xmin = K_FUNC::E_min(xmin, x[ind3]);
  xmax = K_FUNC::E_max(xmax, x[ind3]);
  ymin = K_FUNC::E_min(ymin, y[ind3]);
  ymax = K_FUNC::E_max(ymax, y[ind3]);
  zmin = K_FUNC::E_min(zmin, z[ind3]);
  zmax = K_FUNC::E_max(zmax, z[ind3]);

  // vertex (i,j,k+1)
  xmin = K_FUNC::E_min(xmin, x[ind4]);
  xmax = K_FUNC::E_max(xmax, x[ind4]);
  ymin = K_FUNC::E_min(ymin, y[ind4]);
  ymax = K_FUNC::E_max(ymax, y[ind4]);
  zmin = K_FUNC::E_min(zmin, z[ind4]);
  zmax = K_FUNC::E_max(zmax, z[ind4]);

  // vertex (i+1,j,k+1)
  xmin = K_FUNC::E_min(xmin, x[ind5]);
  xmax = K_FUNC::E_max(xmax, x[ind5]);
  ymin = K_FUNC::E_min(ymin, y[ind5]);
  ymax = K_FUNC::E_max(ymax, y[ind5]);
  zmin = K_FUNC::E_min(zmin, z[ind5]);
  zmax = K_FUNC::E_max(zmax, z[ind5]);

  // vertex (i+1,j+1,k+1)
  xmin = K_FUNC::E_min(xmin, x[ind6]);
  xmax = K_FUNC::E_max(xmax, x[ind6]);
  ymin = K_FUNC::E_min(ymin, y[ind6]);
  ymax = K_FUNC::E_max(ymax, y[ind6]);
  zmin = K_FUNC::E_min(zmin, z[ind6]);
  zmax = K_FUNC::E_max(zmax, z[ind6]);

  // vertex (i,j+1,k+1)
  xmin = K_FUNC::E_min(xmin, x[ind7]);
  xmax = K_FUNC::E_max(xmax, x[ind7]);
  ymin = K_FUNC::E_min(ymin, y[ind7]);
  ymax = K_FUNC::E_max(ymax, y[ind7]);
  zmin = K_FUNC::E_min(zmin, z[ind7]);
  zmax = K_FUNC::E_max(zmax, z[ind7]);

}

//=============================================================================
/* Calcul de la bounding box d'une cellule non structuree */
//=============================================================================
void K_COMPGEOM::boundingBoxOfUnstrCell(
  E_Int noet, FldArrayI& connect, 
  E_Float* xt, E_Float* yt, E_Float* zt,
  E_Float& xmin, E_Float& ymin, E_Float& zmin,
  E_Float& xmax, E_Float& ymax, E_Float& zmax) 
{

  xmin = K_CONST::E_MAX_FLOAT;
  ymin = K_CONST::E_MAX_FLOAT;
  zmin = K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int nvert = connect.getNfld();
  
  for (E_Int vert = 1; vert <= nvert; vert++)
  {
    E_Int ind = connect(noet,vert)-1;
    xmin = K_FUNC::E_min(xmin,xt[ind]);
    ymin = K_FUNC::E_min(ymin,yt[ind]);
    zmin = K_FUNC::E_min(zmin,zt[ind]);
    xmax = K_FUNC::E_max(xmax,xt[ind]);
    ymax = K_FUNC::E_max(ymax,yt[ind]);
    zmax = K_FUNC::E_max(zmax,zt[ind]);
  }
}

//=============================================================================
/* Calcul de la bounding box d'une cellule NGon */
//=============================================================================
void K_COMPGEOM::boundingBoxOfNGonCell(
  E_Int noet, FldArrayI& connect, 
  E_Float* xt, E_Float* yt, E_Float* zt,
  E_Float& xmin, E_Float& ymin, E_Float& zmin,
  E_Float& xmax, E_Float& ymax, E_Float& zmax) 
{

  xmin = K_CONST::E_MAX_FLOAT;
  ymin = K_CONST::E_MAX_FLOAT;
  zmin = K_CONST::E_MAX_FLOAT;
  xmax = -K_CONST::E_MAX_FLOAT;
  ymax = -K_CONST::E_MAX_FLOAT;
  zmax = -K_CONST::E_MAX_FLOAT;

  E_Int* ngon = connect.getNGon();
  E_Int* nface = connect.getNFace();
  E_Int* indPG = connect.getIndPG(); 
  E_Int* indPH = connect.getIndPH();

  E_Int sizeElt;
  E_Int sizeFace;
  E_Int ind;

  E_Int* elt = connect.getElt(noet, sizeElt, nface, indPH);
  for (E_Int i = 0; i < sizeElt; i++)
  {
    E_Int* face = connect.getFace(elt[i]-1, sizeFace, ngon, indPG);
    for (E_Int j = 0; j < sizeFace; j++)
    {
      ind = face[j]-1;
      xmin = K_FUNC::E_min(xmin,xt[ind]);
      ymin = K_FUNC::E_min(ymin,yt[ind]);
      zmin = K_FUNC::E_min(zmin,zt[ind]);
      xmax = K_FUNC::E_max(xmax,xt[ind]);
      ymax = K_FUNC::E_max(ymax,yt[ind]);
      zmax = K_FUNC::E_max(zmax,zt[ind]);
    }
  }
}

//=============================================================================
// Intersection de bbox de 2 grilles structurees
// retourne 1 si les bounding boxes s intersectent, 0 sinon
// similaire a BlkOverlapData::testBBIntersection
//=============================================================================
E_Int K_COMPGEOM::compBoundingBoxIntersection(E_Int ni1, E_Int nj1, E_Int nk1, 
                                        E_Int posx1, E_Int posy1, 
                                        E_Int posz1, FldArrayF& f1, 
                                        E_Int ni2, E_Int nj2, E_Int nk2, 
                                        E_Int posx2, E_Int posy2, 
                                        E_Int posz2, FldArrayF& f2, 
                                        E_Float& xmin1, E_Float& xmax1, 
                                        E_Float& ymin1, E_Float& ymax1, 
                                        E_Float& zmin1, E_Float& zmax1, 
                                        E_Float& xmin2, E_Float& xmax2, 
                                        E_Float& ymin2, E_Float& ymax2, 
                                        E_Float& zmin2, E_Float& zmax2,
                                        E_Float tol)
{
  // bbox1 ds repere absolu

  boundingBoxStruct(ni1, nj1, nk1,
    f1.begin(posx1), f1.begin(posy1), f1.begin(posz1),
    xmin1, ymin1, zmin1, xmax1, ymax1, zmax1);
  
  //bbox2 ds repere absolu
  boundingBoxStruct(ni2, nj2, nk2,
    f2.begin(posx2), f2.begin(posy2), f2.begin(posz2),
    xmin2, ymin2, zmin2, xmax2, ymax2, zmax2);

  if ( xmin1  <=  xmax2+tol && xmax1  >=  xmin2-tol &&
       ymin1  <=  ymax2+tol && ymax1  >=  ymin2-tol &&
       zmin1  <=  zmax2+tol && zmax1  >=  zmin2-tol )
    return 1;
  else return 0;
}
//=============================================================================
// Intersection de 2 bbox donnees
// retourne 1 si les bounding boxes s'intersectent, 0 sinon
//=============================================================================
E_Int K_COMPGEOM::compBoundingBoxIntersection(E_Float xmin1, E_Float xmax1, 
                                        E_Float ymin1, E_Float ymax1, 
                                        E_Float zmin1, E_Float zmax1, 
                                        E_Float xmin2, E_Float xmax2, 
                                        E_Float ymin2, E_Float ymax2, 
                                        E_Float zmin2, E_Float zmax2,
                                        E_Float tol)
{
  if (xmin1 > xmax2+tol) return 0;
  if (xmax1 < xmin2-tol) return 0;
  if (ymin1 > ymax2+tol) return 0;
  if (ymax1 < ymin2-tol) return 0;
  if (zmin1 > zmax2+tol) return 0;
  if (zmax1 < zmin2-tol) return 0;
  return 1;
}
//=============================================================================
/* Recherche si les CEBB  de 2 grilles s intersectent */
//=============================================================================
E_Int K_COMPGEOM::compCEBBIntersection(E_Int ni1, E_Int nj1, E_Int nk1, 
                                 E_Int posx1, E_Int posy1, E_Int posz1,
                                 FldArrayF& f1, 
                                 E_Int ni2, E_Int nj2, E_Int nk2,
                                 E_Int posx2, E_Int posy2, E_Int posz2,
                                 FldArrayF& f2, E_Float tol)
{
  // 1- test des bbox globales
  E_Float xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
  E_Float xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;

  E_Int isIntersect = compBoundingBoxIntersection(
    ni1, nj1, nk1, posx1, posy1, posz1, f1,
    ni2, nj2, nk2, posx2, posy2, posz2, f2, 
    xmin1, xmax1, ymin1, ymax1, zmin1, zmax1,
    xmin2, xmax2, ymin2, ymax2, zmin2, zmax2, tol);

  if ( isIntersect == 0 ) return 0;

  // 2 - test des CEBB 2 a 2
  FldArrayF cartMin1;
  FldArrayF cartMax1;
  FldArrayF cartMin2;
  FldArrayF cartMax2;

  // creation des tableaux de coordonnees
  E_Int npts1 = ni1*nj1*nk1;
  FldArrayF coord1(npts1, 3);
  E_Float* c1x = coord1.begin(1);
  E_Float* c1y = coord1.begin(2);
  E_Float* c1z = coord1.begin(3);
  E_Float* f1x = f1.begin(posx1);
  E_Float* f1y = f1.begin(posy1);
  E_Float* f1z = f1.begin(posz1);
  for (E_Int i1 = 0; i1 < npts1; i1++)
  {
    c1x[i1] = f1x[i1];
    c1y[i1] = f1y[i1];
    c1z[i1] = f1z[i1];
  }
  
  E_Int npts2 = ni2*nj2*nk2;
  FldArrayF coord2(npts2, 3);
  E_Float* c2x = coord2.begin(1);
  E_Float* c2y = coord2.begin(2);
  E_Float* c2z = coord2.begin(3);
  E_Float* f2x = f2.begin(posx2);
  E_Float* f2y = f2.begin(posy2);
  E_Float* f2z = f2.begin(posz2);
  for (E_Int i2 = 0; i2 < npts2; i2++)
  {
    c2x[i2] = f2x[i2];
    c2y[i2] = f2y[i2];
    c2z[i2] = f2z[i2];
  }
  
  FldArrayF cartEltArray1;
  FldArrayF cartEltArray2;
  E_Bool isDegenerated = true;

  for (E_Int dir1 = 1; dir1 <= 3; dir1++)
  {
    short isok1 = compCartEltsArray(
      dir1, ni1, nj1, nk1, xmin1, ymin1, zmin1, 
      xmax1, ymax1, zmax1, coord1, cartEltArray1);
    
    if ( isok1 == 1 ) // direction possible pour les elts cart 
    { 
      for (E_Int dir2 = 1; dir2 <= 3; dir2++)
      {
        short isok2 = compCartEltsArray(
          dir2, ni2, nj2, nk2, xmin2, ymin2, zmin2, 
          xmax2, ymax2, zmax2, coord2, cartEltArray2);

        if (isok2 == 1)
        {
          isDegenerated = false;
          // test des intersection des cebboxes
          for (E_Int elt1 = 0; elt1 < cartEltArray1.getSize(); elt1++)
          {
            E_Float xmin11 = cartEltArray1(elt1,1);
            E_Float ymin11 = cartEltArray1(elt1,2);
            E_Float zmin11 = cartEltArray1(elt1,3);
            E_Float xmax11 = cartEltArray1(elt1,4);
            E_Float ymax11 = cartEltArray1(elt1,5);
            E_Float zmax11 = cartEltArray1(elt1,6);
            
            for (E_Int elt2 = 0; elt2 < cartEltArray2.getSize(); elt2++)
            {
              E_Float xmin21 = cartEltArray2(elt2,1);
              E_Float ymin21 = cartEltArray2(elt2,2);
              E_Float zmin21 = cartEltArray2(elt2,3);
              E_Float xmax21 = cartEltArray2(elt2,4);
              E_Float ymax21 = cartEltArray2(elt2,5);
              E_Float zmax21 = cartEltArray2(elt2,6);
            
              // test intersection
              if ( xmin11 <=  xmax21+tol && xmax11 >=  xmin21-tol &&
                   ymin11 <=  ymax21+tol && ymax11 >=  ymin21-tol &&
                   zmin11 <=  zmax21+tol && zmax11 >=  zmin21-tol )
                return 1;
            }
          }
        }
      }
    } 
  }
  // Dans le cas ou on ne peut pas calculer la CEBB, on renvoit
  // le resultat de l'intersection des BB
  if (isDegenerated == true) return 1;
  return 0;
}
//============================ CompGeom/compGeom.cpp =======================

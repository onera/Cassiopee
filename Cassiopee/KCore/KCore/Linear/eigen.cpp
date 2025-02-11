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

#include "linear.h"
#include "Nuga/include/DelaunayMath.h"
# include <math.h>

#define SWAP(a, b) { temp = a; a = b; b = temp; }

//=============================================================================
// Symmetric Householder reduction to tridiagonal form.
//==============================================================================
void K_LINEAR::tred2(E_Float V[3][3], E_Float d[3], E_Float e[3]) 
{
//  This is derived from the Algol procedures tred2 by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  E_Float scale, h;

  for (E_Int j = 0; j < 3; j++) { d[j] = V[3-1][j]; }

  // Householder reduction to tridiagonal form.
  for (E_Int i = 3-1; i > 0; i--) 
  {
    // Scale to avoid under/overflow.
    scale = 0.0; h = 0.0;
    for (E_Int k = 0; k < i; k++) { scale += K_FUNC::E_abs(d[k]); }
    if (K_FUNC::fEqualZero(scale) == true)
    {
      e[i] = d[i-1];
      for (E_Int j = 0; j < i; j++) 
      {
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
        V[j][i] = 0.0;
      }
    } 
    else 
    {
      // Generate Householder vector.
      for (E_Int k = 0; k < i; k++) 
      {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      E_Float f = d[i-1];
      E_Float g = sqrt(h);
      if (f > 0) { g = -g; }
      e[i] = scale * g;
      h = h - f * g;
      d[i-1] = f - g;
      for (E_Int j = 0; j < i; j++) { e[j] = 0.0; }

      // Apply similarity transformation to remaining columns.
      for (E_Int j = 0; j < i; j++) 
      {
        f = d[j];
        V[j][i] = f;
        g = e[j] + V[j][j] * f;
        for (E_Int k = j+1; k <= i-1; k++) 
        {
          g += V[k][j] * d[k];
          e[k] += V[k][j] * f;
        }
        e[j] = g;
      }
      f = 0.0;
      for (E_Int j = 0; j < i; j++) 
      {
        e[j] /= h;
        f += e[j] * d[j];
      }
      E_Float hh = f / (h + h);
      for (E_Int j = 0; j < i; j++) { e[j] -= hh * d[j]; }
      for (E_Int j = 0; j < i; j++) 
      {
        f = d[j];
        g = e[j];
        for (E_Int k = j; k <= i-1; k++) { V[k][j] -= (f*e[k] + g*d[k]); }
        d[j] = V[i-1][j];
        V[i][j] = 0.0;
      }
    }
    d[i] = h;
  }

  // Accumulate transformations.
  for (E_Int i = 0; i < 3-1; i++) 
  {
    V[3-1][i] = V[i][i];
    V[i][i] = 1.0;
    h = d[i+1];
    if (K_FUNC::fEqualZero(h) == false) 
    {
      for (E_Int k = 0; k <= i; k++) { d[k] = V[k][i+1] / h; }
      for (E_Int j = 0; j <= i; j++) 
      {
        E_Float g = 0.0;
        for (E_Int k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
        for (E_Int k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
      }
    }
    for (E_Int k = 0; k <= i; k++) { V[k][i+1] = 0.0; }
  }
  for (E_Int j = 0; j < 3; j++) 
  {
    d[j] = V[3-1][j];
    V[3-1][j] = 0.0;
  }
  V[3-1][3-1] = 1.0;
  e[0] = 0.0;
} 

//=============================================================================
// Symmetric tridiagonal QL algorithm.
//==============================================================================
void K_LINEAR::tql2(E_Float V[3][3], E_Float d[3], E_Float e[3]) 
{
//  This is derived from the Algol procedures tql2, by
//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutine in EISPACK.

  for (E_Int i = 1; i < 3; i++) { e[i-1] = e[i]; }
  e[3-1] = 0.0;

  E_Float f = 0.0;
  E_Float tst1 = 0.0;
  E_Float eps = pow(2.0,-52.0);
  for (E_Int l = 0; l < 3; l++) 
  {
    // Find small subdiagonal element
    tst1 = K_FUNC::E_max(tst1, fabs(d[l]) + fabs(e[l]));
    E_Int m = l;
    while (m < 3) 
    {
      if (fabs(e[m]) <= eps*tst1) { break; }
      m++;
    }

    // If m == l, d[l] is an eigenvalue,
    // otherwise, iterate.

    if (m > l) 
    {
      E_Int iter = 0;
      do {
        iter += 1;  // (Could check iteration count here.)

        // Compute implicit shift
        E_Float g = d[l];
        E_Float p = (d[l+1] - g) / (2.0 * e[l]);
        E_Float r = sqrt(p*p+1.);
        if (p < 0) { r = -r; }
        d[l] = e[l] / (p + r);
        d[l+1] = e[l] * (p + r);
        E_Float dl1 = d[l+1];
        E_Float h = g - d[l];
        for (E_Int i = l+2; i < 3; i++) { d[i] -= h; }
        f = f + h;

        // Implicit QL transformation.
        p = d[m];
        E_Float c = 1.0;
        E_Float c2 = c;
        E_Float c3 = c;
        E_Float el1 = e[l+1];
        E_Float s = 0.0;
        E_Float s2 = 0.0;
        for (E_Int i = m-1; i >= l; i--) 
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * p;
          r = sqrt(p*p+e[i]*e[i]);
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * (c * g + s * d[i]);

          // Accumulate transformation.
          for (E_Int k = 0; k < 3; k++) 
          {
            h = V[k][i+1];
            V[k][i+1] = s * V[k][i] + c * h;
            V[k][i] = c * V[k][i] - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;

        // Check for convergence.
      } while (fabs(e[l]) > eps*tst1);
    }
    d[l] = d[l] + f;
    e[l] = 0.0;
  }
  
  // Sort eigenvalues and corresponding vectors.

  for (E_Int i = 0; i < 3-1; i++) 
  {
    E_Int k = i;
    E_Float p = d[i];
    for (E_Int j = i+1; j < 3; j++) 
    {
      if (d[j] < p) 
      {
        k = j;
        p = d[j];
      }
    }
    if (k != i) 
    {
      d[k] = d[i];
      d[i] = p;
      for (E_Int j = 0; j < 3; j++) 
      {
        p = V[j][i];
        V[j][i] = V[j][k];
        V[j][k] = p;
      }
    }
  }
}

//=============================================================================
/* Eigen value and vector of a 2x2 symetric matrix
*/
//=============================================================================
E_Int K_LINEAR::eigen2(E_Float a00, E_Float a11, E_Float a10, 
                       E_Float& lambda0, E_Float& lambda1, 
                       E_Float* v0, E_Float* v1)
{
  DelaunayMath::eigen_vectors(a00, a11, a10, lambda0, lambda1,
                              v0, v1);
  return 1;
}

//=============================================================================
/* Eigen value and vector of a 3x3 symetric matrix
*/
//=============================================================================
E_Int K_LINEAR::eigen3(E_Float a00, E_Float a01, E_Float a02,
                       E_Float a11, E_Float a12, E_Float a22,
                       E_Float& lambda0, E_Float& lambda1, E_Float& lambda2, 
                       E_Float* v0, E_Float* v1, E_Float* v2)
{
  root3(a00, a01, a02, a11, a12, a22, lambda0, lambda1, lambda2);
  E_Float b00, b01, b02;
  E_Float b10, b11, b12;
  E_Float b20, b21, b22;
  E_Float c00, c01, c02, c11, c12, c22;
  c00 = a00 - lambda0; c01 = a01; c02 = a02;
  c11 = a11 - lambda0; c12 = a12; c22 = a22 - lambda0;
  E_Int rank = rank3(c00, c01, c02, c11, c12, c22,
                     b00, b01, b02, b10, b11, b12, b20, b21, b22);
  if (rank == 0)
  {
    v0[0] = 1.; v0[1] = 0.; v0[2] = 0.;
    v1[0] = 0.; v1[1] = 1.; v1[2] = 0.;
    v2[0] = 0.; v2[1] = 0.; v2[2] = 1.;
    return rank;
  }
  if (rank == 1)
  {
    getComplement32(b00, b01, b02, v0[0], v0[1], v0[2], v1[0], v1[1], v1[2]);
    v2[0] = v0[1]*v1[2]-v0[2]*v1[1];
    v2[1] = v0[2]*v1[0]-v0[0]*v1[2];
    v2[2] = v0[0]*v1[1]-v0[1]*v1[0];
    return rank;
  }
  getComplement31(b00, b01, b02, b10, b11, b12, v0[0], v0[1], v0[2]);
  c00 = a00 - lambda1; c01 = a01; c02 = a02;
  c11 = a11 - lambda1; c12 = a12; c22 = a22 - lambda1;
  rank = rank3(c00, c01, c02, c11, c12, c22,
               b00, b01, b02, b10, b11, b12, b20, b21, b22);
  if (rank == 1)
  {
    getComplement32(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]);
    return rank;
  }

  // rank == 2
  getComplement31(b00, b01, b02, b10, b11, b12, v1[0], v1[1], v1[2]);
  v2[0] = v0[1]*v1[2]-v0[2]*v1[1];
  v2[1] = v0[2]*v1[0]-v0[0]*v1[2];
  v2[2] = v0[0]*v1[1]-v0[1]*v1[0];
  return rank;
}

//=============================================================================
/* Eigen value and vector of a 3x3 symetric matrix (other version)
*/
//=============================================================================
E_Int K_LINEAR::eigen3bis(E_Float A[3][3], E_Float V[3][3], E_Float d[3])
{
  E_Float e[3];
  for (E_Int i = 0; i < 3; i++) 
  {
    for (E_Int j = 0; j < 3; j++) { V[i][j] = A[i][j]; }
  }
  tred2(V, d, e);
  tql2(V, d, e);
  return 0;
}

//==============================================================================
E_Int K_LINEAR::root3(E_Float a00, E_Float a01, E_Float a02,
                      E_Float a11, E_Float a12, E_Float a22,
                      E_Float& lambda0, E_Float& lambda1, E_Float& lambda2)
{
  E_Float inv3 = 1./3.;
  E_Float root3 = sqrt(3.);
  E_Float c0 = a00*a11*a22 + 2.*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
  E_Float c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12;
  E_Float c2 = a00 + a11 + a22;

  E_Float c2Div3 = c2*inv3;
  E_Float aDiv3 = (c1 - c2*c2Div3)*inv3;
  if (aDiv3 > 0.) aDiv3 = 0.;
  E_Float mbDiv2 = 0.5*(c0 + c2Div3*(2.*c2Div3*c2Div3 - c1));
  E_Float q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3;
  if (q > 0.) q = 0.;
  E_Float magnitude = sqrt(-aDiv3);
  E_Float angle = atan2(sqrt(-q), mbDiv2)*inv3;
  E_Float cs = cos(angle);
  E_Float sn = sin(angle);
  lambda0 = c2Div3 + 2.*magnitude*cs;
  lambda1 = c2Div3 - magnitude*(cs + root3*sn);
  lambda2 = c2Div3 - magnitude*(cs - root3*sn);
  E_Float temp;
  if (lambda1 < lambda0) SWAP(lambda0, lambda1);
  if (lambda2 < lambda0) SWAP(lambda0, lambda2);
  if (lambda2 < lambda1) SWAP(lambda1, lambda2);
  return 1;
}

//==============================================================================
void K_LINEAR::getComplement31(
  E_Float u0, E_Float u1, E_Float u2,
  E_Float v0, E_Float v1, E_Float v2,
  E_Float& w0, E_Float& w1, E_Float& w2)
{
  E_Float n;
  w0 = u1*v2-u2*v1;
  w1 = u2*v0-u0*v2;
  w2 = u0*v1-u1*v0;
  n = w0*w0+w1*w1+w2*w2;
  n = sqrt(n);
  n = 1./K_FUNC::E_max(n, 1.e-10);
  w0 *= n; w1 *= n; w2 *= n;
}

//=============================================================================
void K_LINEAR::getComplement32(
  E_Float u0, E_Float u1, E_Float u2,
  E_Float& v0, E_Float& v1, E_Float& v2,
  E_Float& w0, E_Float& w1, E_Float& w2)
{
  E_Float n;
  n = u0*u0+u1*u1+u2*u2;
  n = sqrt(n);
  n = 1./K_FUNC::E_max(n, 1.e-10);
  u0 *= n; u1 *= n; u2 *= n;
  if (K_FUNC::E_abs(u0) >= K_FUNC::E_abs(u1))
  {
    n = sqrt(u0*u0+u2*u2);
    n = 1./K_FUNC::E_max(n, 1.e-10);
    v0 = -u2*n;
    v1 = 0.;
    v2 = u0*n;
    w0 = u1*v2;
    w1 = u2*v0-u0*v2;
    w2 = -u1*v0;
  }
  else
  {
    n = sqrt(u1*u1+u2*u2);
    n = 1./K_FUNC::E_max(n, 1.e-10);
    v0 = 0.;
    v1 = u2*n;
    v2 = -u1*n;
    w0 = u1*v2-u2*v1;
    w1 = -u0*v2;
    w2 = u0*v1;
  }
}

//=============================================================================
//
//=============================================================================
E_Int K_LINEAR::rank3(E_Float a00, E_Float a01, E_Float a02,
                      E_Float a11, E_Float a12, E_Float a22,
                      E_Float& b00, E_Float& b01, E_Float& b02,
                      E_Float& b10, E_Float& b11, E_Float& b12,
                      E_Float& b20, E_Float& b21, E_Float& b22)
{
  // Copie la matrice dans b
  b00 = a00; b01 = a01; b02 = a02;
  b10 = a01; b11 = a11; b12 = a12;
  b20 = a02; b21 = a12; b22 = a22;

  E_Float temp, n;
  // Calcul la plus grande valeur de la matrice en valeur absolue
  E_Float c00 = K_FUNC::E_abs(a00);
  E_Float c01 = K_FUNC::E_abs(a01);
  E_Float c02 = K_FUNC::E_abs(a02);
  E_Float c11 = K_FUNC::E_abs(a11);
  E_Float c12 = K_FUNC::E_abs(a12);
  E_Float c22 = K_FUNC::E_abs(a22);
  E_Float maxv = c00;
  E_Int maxRow = 0; E_Int maxCol = 0;
  if (c01 > maxv) { maxv = c01; maxRow = 0; maxCol = 1; }
  if (c02 > maxv) { maxv = c02; maxRow = 0; maxCol = 2; }
  if (c11 > maxv) { maxv = c11; maxRow = 1; maxCol = 1; }
  if (c12 > maxv) { maxv = c12; maxRow = 1; maxCol = 2; }
  if (c22 > maxv) { maxv = c22; maxRow = 2; maxCol = 2; }
  if (maxv < 1.e-5)
  {
    return 0; // rank=0, multiciplicity 3
  }
  if (maxRow == 1)
  {
    SWAP(b00, b10);
    SWAP(b01, b11);
    SWAP(b02, b12);
  }
  else if (maxRow == 2)
  {
    SWAP(b00, b20);
    SWAP(b01, b21);
    SWAP(b02, b22);
  }

  // Reduce
  if (maxCol == 0) n = 1./b00;
  else if (maxCol == 1) n = 1./b01;
  else n = 1./b02; // maxCol=2
  b00 *= n; b01 *= n; b02 *= n;
  if (maxCol == 0)
  {
    b11 -= b10*b01;
    b12 -= b10*b02;
    b21 -= b20*b01;
    b10 = 0.;
    b20 = 0.;
  }
  else if (maxCol == 1)
  {
    b10 -= b11*b00;
    b12 -= b11*b02;
    b20 -= b21*b00;
    b22 -= b21*b02;
    b11 = 0.;
    b21 = 0.;
  }
  else
  {
    b10 -= b12*b00;
    b11 -= b12*b01;
    b20 -= b22*b00;
    b21 -= b22*b01;
    b12 = 0.;
    b22 = 0.;
  }

  // max des 2 dernieres lignes
  E_Float c10 = K_FUNC::E_abs(b10);
  c11 = K_FUNC::E_abs(b11);
  c12 = K_FUNC::E_abs(b12);
  E_Float c20 = K_FUNC::E_abs(b20);
  E_Float c21 = K_FUNC::E_abs(b21);
  c22 = K_FUNC::E_abs(b22);
  maxv = c10;
  maxRow = 1; maxCol = 0;
  if (c11 > maxv) { maxv = c11; maxRow = 1; maxCol = 1; }
  if (c12 > maxv) { maxv = c12; maxRow = 1; maxCol = 2; }
  if (c20 > maxv) { maxv = c20; maxRow = 2; maxCol = 0; }
  if (c21 > maxv) { maxv = c21; maxRow = 2; maxCol = 1; }
  if (c22 > maxv) { maxv = c22; maxRow = 2; maxCol = 2; }
  if (maxv < 1.e-5)
  {
    return 1; 
  }

  if (maxRow == 2)
  {
    SWAP(b10, b20);
    SWAP(b11, b21);
    SWAP(b12, b22);
  }

  if (maxCol == 0) n = 1./b10;
  else if (maxCol == 1) n = 1./b11;
  else n = 1./b12;
  b10 *= n; b11 *= n; b12 *= n;
  return 2;
}

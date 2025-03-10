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
#include <cmath>
#include <cstring>

void K_LINEAR::sym3mat_eigen(const E_Float M[6], E_Float L[3],
  E_Float v1[3], E_Float v2[3], E_Float v3[3], const E_Float tol)
{
  // Init
  v1[0] = 1.0; v1[1] = 0.0; v1[2] = 0.0;
  v2[0] = 0.0; v2[1] = 1.0; v2[2] = 0.0;
  v3[0] = 0.0; v3[1] = 0.0; v3[2] = 1.0;
  L[0] = M[0];
  L[1] = M[3];
  L[2] = M[5];
  
  E_Float Ap[6], B[6], Ar[4];
  E_Float s1[3], s2[3], t2[3], t3[3], r1[3], r2[3], r3[3], tmp1[3], tmp2[3],
    u1[3], u2[3], w1[3];
  E_Float J2, J3, alpha, thirdTrA, norm1, norm2, norm3, coeff, dif, sum, sgn;
  
  // Check null matrix
  E_Float maxm = fabs(M[0]);
  E_Float valm;
  for (E_Int i = 1; i < 6; i++) {
    valm = fabs(M[i]);
    if (valm > maxm) maxm = valm;
  }
  if (maxm < 5e-6) return;
  //if (maxm < tol) return;
  
  // Normalize matrix
  E_Float dd = 1. / maxm;
  E_Float A[6];
  memcpy(A, M, 6*sizeof(E_Float));
  for (E_Int i = 0; i < 6; i++) A[i] *= dd;
  
  // Check for diagonal matrix
  E_Float maxd = fabs(A[1]);
  valm = fabs(A[2]);
  if (valm > maxd) maxd = valm;
  valm = fabs(A[4]);
  if (valm > maxd) maxd = valm;
  if (maxd < 1e-13) return;
  //if (maxd < tol) return; // Off-diagonal coeffs are smaller than tol

  thirdTrA = (A[0] + A[3] + A[5]) / 3.0;
  Ap[0] = A[0] - thirdTrA;
  Ap[1] = A[1];
  Ap[2] = A[2];
  Ap[3] = A[3] - thirdTrA;
  Ap[4] = A[4];
  Ap[5] = A[5] - thirdTrA;
  J2 = K_MATH::sym3mat_second_invariant(Ap);
  J3 = K_MATH::sym3mat_third_invariant(Ap);
  E_Float tmp = 0.5*J3 * pow(3.0/J2, 1.5);
  if (tmp > 1.0) tmp = 1.0;
  else if (tmp < -1.0) tmp = -1.0;
  alpha = acos(tmp) / 3.0;
  
  if (alpha < K_CONST::E_PI/6.0)
    L[0] = 2.0*sqrt(J2/3.0)*cos(alpha);
  else
    L[0] = 2.0*sqrt(J2/3.0)*cos(alpha + 4.0*K_CONST::E_PI/3.0);
  
  // Find eigenvector corresponding to L[0]
  B[0] = Ap[0] - L[0];
  B[1] = Ap[1];
  B[2] = Ap[2];
  B[3] = Ap[3] - L[0];
  B[4] = Ap[4];
  B[5] = Ap[5] - L[0];
  r1[0] = B[0]; r1[1] = B[1]; r1[2] = B[2];
  r2[0] = B[1]; r2[1] = B[3]; r2[2] = B[4];
  r3[0] = B[2]; r3[1] = B[4]; r3[2] = B[5];
  norm1 = K_MATH::norm(r1, 3);
  norm2 = K_MATH::norm(r2, 3);
  norm3 = K_MATH::norm(r3, 3);
  E_Float over_norm;
  if (norm1 >= norm2 && norm1 >= norm3) {
    over_norm = 1. / norm1;
    s1[0] = r1[0] * over_norm;
    s1[1] = r1[1] * over_norm;
    s1[2] = r1[2] * over_norm;

    coeff = K_MATH::dot(s1, r2, 3);
    t2[0] = r2[0] - coeff*s1[0];
    t2[1] = r2[1] - coeff*s1[1];
    t2[2] = r2[2] - coeff*s1[2];

    coeff = K_MATH::dot(s1, r3, 3);
    t3[0] = r3[0] - coeff*s1[0];
    t3[1] = r3[1] - coeff*s1[1];
    t3[2] = r3[2] - coeff*s1[2];
  } else if (norm2 >= norm1 && norm2 >= norm3) {
    over_norm = 1. / norm2;
    s1[0] = r2[0] * over_norm;
    s1[1] = r2[1] * over_norm;
    s1[2] = r2[2] * over_norm;

    coeff = K_MATH::dot(s1, r1, 3);
    t2[0] = r1[0] - coeff*s1[0];
    t2[1] = r1[1] - coeff*s1[1];
    t2[2] = r1[2] - coeff*s1[2];
    coeff = K_MATH::dot(s1, r3, 3);
    t3[0] = r3[0] - coeff*s1[0];
    t3[1] = r3[1] - coeff*s1[1];
    t3[2] = r3[2] - coeff*s1[2];
  } else {
    over_norm = 1. / norm3;
    s1[0] = r3[0] * over_norm;
    s1[1] = r3[1] * over_norm;
    s1[2] = r3[2] * over_norm;

    coeff = K_MATH::dot(s1, r2, 3);
    t2[0] = r2[0] - coeff*s1[0];
    t2[1] = r2[1] - coeff*s1[1];
    t2[2] = r2[2] - coeff*s1[2];
    coeff = K_MATH::dot(s1, r1, 3);
    t3[0] = r1[0] - coeff*s1[0];
    t3[1] = r1[1] - coeff*s1[1];
    t3[2] = r1[2] - coeff*s1[2];
  }
  
  norm2 = K_MATH::norm(t2, 3);
  norm3 = K_MATH::norm(t3, 3);
  if (norm2 >= norm3) {
    over_norm = 1. / norm2;
    s2[0] = t2[0] * over_norm;
    s2[1] = t2[1] * over_norm;
    s2[2] = t2[2] * over_norm;
  } else {
    over_norm = 1. / norm3;
    s2[0] = t3[0] * over_norm;
    s2[1] = t3[1] * over_norm;
    s2[2] = t3[2] * over_norm;   
  }

  K_MATH::cross(s1, s2, v1);
  norm1 = K_MATH::norm(v1, 3);
  v1[0] = v1[0] / norm1; v1[1] = v1[1] / norm1; v1[2] = v1[2] / norm1;
  
  // Reduced form of Ap
  K_MATH::sym3mat_dot_vec(Ap, s1, tmp1);
  K_MATH::sym3mat_dot_vec(Ap, s2, tmp2);
  Ar[0] = K_MATH::dot(s1, tmp1, 3);
  Ar[1] = K_MATH::dot(s1, tmp2, 3);
  Ar[2] = K_MATH::dot(s2, tmp1, 3);
  Ar[3] = K_MATH::dot(s2, tmp2, 3);
  
  // Wilkinson shift
  dif = Ar[0] - Ar[3];
  sum = Ar[0] + Ar[3];
  sgn = K_MATH::sign(dif);
  L[1] = 0.5*sum - 0.5*sgn*sqrt(dif*dif + 4.*Ar[1]*Ar[2]);
  L[2] = sum - L[1];
  
  // find eigenvector corresponding to L1
  B[0] = Ap[0] - L[1];
  B[1] = Ap[1];
  B[2] = Ap[2];
  B[3] = Ap[3] - L[1];
  B[4] = Ap[4];
  B[5] = Ap[5] - L[1];
  K_MATH::sym3mat_dot_vec(B, s1, u1);
  K_MATH::sym3mat_dot_vec(B, s2, u2);
  
  norm1 = K_MATH::norm(u1, 3);
  norm2 = K_MATH::norm(u2, 3);
  if (norm1 >= norm2) {
    over_norm = 1. / norm1;
    w1[0] = u1[0] * over_norm;
    w1[1] = u1[1] * over_norm;
    w1[2] = u1[2] * over_norm;
  } else {
    over_norm = 1. / norm2;
    w1[0] = u2[0] * over_norm;
    w1[1] = u2[1] * over_norm;
    w1[2] = u2[2] * over_norm;
  }

  K_MATH::cross(w1, v1, v2);
  norm2 = K_MATH::norm(v2, 3);
  v2[0] = v2[0] / norm2;
  v2[1] = v2[1] / norm2;
  v2[2] = v2[2] / norm2;

  K_MATH::cross(v1, v2, v3);
  norm3 = K_MATH::norm(v3, 3);
  v3[0] = v3[0] / norm3;
  v3[1] = v3[1] / norm3;
  v3[2] = v3[2] / norm3;

  for (E_Int i = 0; i < 3; i++)
    L[i] = (L[i] + thirdTrA) * maxm;
}

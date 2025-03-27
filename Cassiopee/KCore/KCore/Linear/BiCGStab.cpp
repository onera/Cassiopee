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
#include <cstring>
#include <cassert>

E_Int K_LINEAR::BiCGStab(const E_Float *A, const E_Float *b, const E_Int n,
  E_Float *x, const E_Float tol)
{
  for (E_Int i = 0; i < n; i++) x[i] = 0.0;

  E_Int converged = 1;

  E_Float bnorm2 = K_MATH::norm(b, n);
  if (K_MATH::feq(bnorm2, 0.0)) bnorm2 = 1.0;

  E_Float* r = (E_Float*)calloc(n, sizeof(E_Float));
  E_Float* Ax = (E_Float*)calloc(n, sizeof(E_Float));

  K_MATH::sqrmat_dot_vec(A, x, n, Ax);
  for (E_Int i = 0; i < n; i++) r[i] = b[i] - Ax[i];

  E_Float err = K_MATH::norm(r, n);
  if (err < tol) 
  {
    free(r);
    free(Ax);
    return converged;
  }

  E_Float omega, rho, beta, rho0, alpha;
  omega = rho = beta = rho0 = alpha = 1.;

  E_Float* rhat = (E_Float *)calloc(n, sizeof(E_Float));
  memcpy(rhat, r, n*sizeof(E_Float));

  E_Float* p = (E_Float*)calloc(n, sizeof(E_Float));
  E_Float* v = (E_Float*)calloc(n, sizeof(E_Float));
  E_Float* s = (E_Float*)calloc(n, sizeof(E_Float));
  E_Float* t = (E_Float*)calloc(n, sizeof(E_Float));

  // TODO(Imad): maxiter should be exactly n theoretically...
  E_Int maxiter = n*10;
  for (E_Int i = 1; i <= maxiter; i++) 
  {
    rho = K_MATH::dot(rhat, r, n);

    if (K_MATH::feq(rho, 0.0)) 
    {
      converged = 0;
      break;
    }

    if (i > 1) {
      beta = rho/rho0 * alpha/omega;
      for (E_Int j = 0; j < n; j++) p[j] = r[j] + beta*(p[j] - omega*v[j]);
    } else {
      memcpy(p, r, n*sizeof(E_Float));
    }

    K_MATH::sqrmat_dot_vec(A, p, n, v);

    alpha = rho / K_MATH::dot(rhat, v, n);

    for (E_Int j = 0; j < n; j++) s[j] = r[j] - alpha*v[j];

    err = K_MATH::norm(s, n);
    if (err < tol) {
      for (E_Int j = 0; j < n; j++) x[j] += alpha*p[j];
      err /= bnorm2;
      converged = 1;
      break;
    }

    K_MATH::sqrmat_dot_vec(A, s, n, t);

    omega = K_MATH::dot(t, s, n) / K_MATH::dot(t, t, n);

    for (E_Int j = 0; j < n; j++) x[j] += alpha*p[j] + omega*s[j];

    for (E_Int j = 0; j < n; j++) r[j] = s[j] - omega*t[j];

    err = K_MATH::norm(r, n);
    if (err < tol) {
      converged = 1;
      break;
    }

    rho0 = rho;
    assert(!K_MATH::feq(omega, 0.0));
  }

  free(r);
  free(Ax);
  free(s);
  free(t);
  free(p);
  free(v);
  free(rhat);

  return converged;
}

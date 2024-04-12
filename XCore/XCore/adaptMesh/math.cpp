/*    
    Copyright 2013-2024 Onera.

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
#include "proto.h"

E_Float dot(const E_Float *a, const E_Float *b, E_Int n)
{
  E_Float res = 0;
  for (E_Int i = 0; i < n; i++)
    res += a[i]*b[i];
  return res;
}

void cross(E_Float a[3], E_Float b[3], E_Float c[3])
{
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

E_Float norm(const E_Float *a, E_Int n)
{
  return sqrt(dot(a, a, n));
}

#define EPS 1e-37

E_Int feq(E_Float a, E_Float b)
{
  if (fabs(a-b) < EPS) return 1;
  return 0;
}

void sqrmat_times_vec(E_Float *A, E_Float *x, E_Float *Ax, E_Int n)
{
  E_Int i, j;
  for (i = 0; i < n; i++) {
    Ax[i] = 0.;
    for (j = 0; j < n; j++)
      Ax[i] += A[j*n + i]*x[j];
  }
}

E_Int BiCGStab(E_Float *A, E_Float *x, E_Float *b, E_Int n)
{
  // initial guess
  E_Int i;
  for (i = 0; i < n; i++) x[i] = 0.;
  const E_Float tol = 1.e-6;

  E_Int converged = 0;

  E_Float bnorm2 = norm(b, n);
  if (feq(bnorm2, 0.)) bnorm2 = 1.;

  E_Int j;


  E_Float *r = (E_Float *)XCALLOC(n, sizeof(E_Float));
  E_Float *Ax = (E_Float *)XCALLOC(n, sizeof(E_Float));

  sqrmat_times_vec(A, x, Ax, n);
  for (i = 0; i < n; i++) r[i] = b[i] - Ax[i];

  E_Float err = norm(r, n);
  if (err < tol) {
    XFREE(r);
    XFREE(Ax);
    return 1;
  }

  E_Float omega, rho, beta, rho0, alpha;
  omega = rho = beta = rho0 = alpha = 1.;

  E_Float *rhat = (E_Float *)XCALLOC(n, sizeof(E_Float));
  memcpy(rhat, r, n*sizeof(E_Float));

  E_Float *p = (E_Float *)XCALLOC(n, sizeof(E_Float));
  E_Float *v = (E_Float *)XCALLOC(n, sizeof(E_Float));
  E_Float *s = (E_Float *)XCALLOC(n, sizeof(E_Float));
  E_Float *t = (E_Float *)XCALLOC(n, sizeof(E_Float));

  E_Int maxiter = n*10;
  for (i = 1; i <= maxiter; i++) {
    rho = dot(rhat, r, n);

    if (feq(rho, 0.)) {
      converged = 0;
      break;
    }

    if (i > 1) {
      beta = rho/rho0 * alpha/omega;
      for (j = 0; j < n; j++) p[j] = r[j] + beta*(p[j] - omega*v[j]);
    } else {
      memcpy(p, r, n*sizeof(E_Float));
    }

    sqrmat_times_vec(A, p, v, n);

    alpha = rho / dot(rhat, v, n);

    for (j = 0; j < n; j++) s[j] = r[j] - alpha*v[j];

    err = norm(s, n);
    if (err < tol) {
      for (j = 0; j < n; j++) x[j] += alpha*p[j];
      err /= bnorm2;
      converged = 1;
      break;
    }

    sqrmat_times_vec(A, s, t, n);

    omega = dot(t, s, n) / dot(t, t, n);

    for (j = 0; j < n; j++) x[j] += alpha*p[j] + omega*s[j];

    for (j = 0; j < n; j++) r[j] = s[j] - omega*t[j];

    err = norm(r, n);
    if (err < tol) {
      converged = 1;
      break;
    }

    rho0 = rho;
    assert(!feq(omega, 0.));
  }

  XFREE(r);
  XFREE(Ax);
  XFREE(s);
  XFREE(t);
  XFREE(p);
  XFREE(v);
  XFREE(rhat);

  return converged;
}

#define MAXNEI 6

void compute_lsq_grad_matrices(mesh *M)
{
  if (M->lsqG) return;
  E_Int stride = 3*MAXNEI;
  M->lsqG = (E_Float *)XCALLOC(stride*M->ncells, sizeof(E_Float));
  M->lsqGG = (E_Float *)XCALLOC(9*M->ncells, sizeof(E_Float));
  E_Int *count_neis = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei < 0) continue;
    E_Int own = owner[i];

    E_Float *co = &M->cc[3*own];
    E_Float *cn = &M->cc[3*nei];

    E_Float d[3];
    for (E_Int j = 0; j < 3; j++)
      d[j] = cn[j] - co[j];

    E_Float *lo = &M->lsqG[stride*own + 3*count_neis[own]++];
    E_Float *ln = &M->lsqG[stride*nei + 3*count_neis[nei]++];

    for (E_Int j = 0; j < 3; j++) {
      lo[j] = d[j];
      ln[j] = -d[j];
    }
  }

  // TODO(Imad): account for boundary data

  if (!M->pnei_coords) {
    M->pnei_coords = (E_Float **)XCALLOC(M->nppatches, sizeof(E_Float *));
    for (E_Int i = 0; i < M->nppatches; i++) {
      M->pnei_coords[i] = (E_Float *)XCALLOC(
        3*M->ppatches[i].nfaces, sizeof(E_Float));
    }
    comm_interface_data_d(M, &M->cc[0], 3, M->pnei_coords);
  }

  // proc faces
  for (E_Int i = 0; i < M->nppatches; i++) {
    E_Int npfaces = M->ppatches[i].nfaces;
    E_Int *pfaces = M->ppatches[i].faces;
    E_Float *pnei_coords = M->pnei_coords[i];
    for (E_Int j = 0; j < npfaces; j++) {
      E_Int own = owner[pfaces[j]];

      E_Float *co = &M->cc[3*own];
      E_Float *cn = &pnei_coords[3*j];

      E_Float *lo = &M->lsqG[stride*own + 3*count_neis[own]++];
      assert(count_neis[own] <= MAXNEI);

      for (E_Int k = 0; k < 3; k++)
        lo[k] = cn[k] - co[k];
    }
  }

  // compute tA.A
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Float *pGG = &M->lsqGG[9*i];
    E_Float *pG = &M->lsqG[stride*i];
    E_Int idx = 0;
    for (E_Int j = 0; j < 3; j++) {
      for (E_Int k = 0; k < 3; k++) {
        pGG[idx] = 0;
        for (E_Int l = 0; l < count_neis[i]; l++)
          pGG[idx] += pG[l*3+j] * pG[l*3+k];
        idx++;
      }
    }
  }

  XFREE(count_neis);
}

E_Float *compute_grad(mesh *M, E_Float *fld)
{
  compute_cell_centers(M);
  compute_lsq_grad_matrices(M);

  E_Float *B = (E_Float *)XCALLOC(3, sizeof(E_Float));
  E_Float *G = (E_Float *)XCALLOC(3*M->ncells, sizeof(E_Float));
  E_Float *b = (E_Float *)XCALLOC(MAXNEI*M->ncells, sizeof(E_Float));

  E_Int *count_neis = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  // construct b vector
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei < 0) continue;
    E_Int own = owner[i];
    b[MAXNEI*own + count_neis[own]++] = fld[nei] - fld[own];
    b[MAXNEI*nei + count_neis[nei]++] = fld[own] - fld[nei];
  }

  if (!M->pnei_flds) {
    M->pnei_flds = (E_Float **)XCALLOC(M->nppatches, sizeof(E_Float *));
    for (E_Int i = 0; i < M->nppatches; i++) {
      M->pnei_flds[i] = (E_Float *)XCALLOC(
        M->ppatches[i].nfaces, sizeof(E_Float));
    }
  }
  comm_interface_data_d(M, fld, 1, M->pnei_flds);

  for (E_Int i = 0; i < M->nppatches; i++) {
    E_Float *pnei_fld = M->pnei_flds[i];
    for (E_Int j = 0; j < M->ppatches[i].nfaces; j++) {
      E_Int own = owner[M->ppatches[i].faces[j]];
      b[MAXNEI*own + count_neis[own]++] = pnei_fld[j] - fld[own];
    }
  }

  // construct B vector
  E_Int stride = 3*MAXNEI;
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Float *pG = &M->lsqG[stride*i];
    E_Float *pb = &b[MAXNEI*i];
    for (E_Int j = 0; j < 3; j++) {
      B[j] = 0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pG[3*k+j];
    }

    // solve
    E_Int converged = BiCGStab(&M->lsqGG[9*i], &G[3*i], B, 3);
    assert(converged);
  }

  XFREE(B);
  XFREE(b);
  XFREE(count_neis);

  return G;
}

void compute_lsq_hess_matrices(mesh *M)
{
  if (M->lsqH) return;
  E_Int stride = 6*MAXNEI;
  M->lsqH = (E_Float *)XCALLOC(stride*M->ncells, sizeof(E_Float));
  M->lsqHH = (E_Float *)XCALLOC(36*M->ncells, sizeof(E_Float));
  E_Int *count_neis = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  E_Float *lsqH = M->lsqH;
  E_Float *lsqHH = M->lsqHH;

  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  // E_Internal faces
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei < 0) continue;

    E_Int own = owner[i];

    E_Float *co = &M->cc[3*own];
    E_Float *cn = &M->cc[3*nei];

    E_Float dx = cn[0] - co[0];
    E_Float dy = cn[1] - co[1];
    E_Float dz = cn[2] - co[2];

    E_Float *lo = &lsqH[stride*own + 6*count_neis[own]++];
    E_Float *ln = &lsqH[stride*nei + 6*count_neis[nei]++];

    lo[0] = dx*dx;
    lo[1] = 2.*dx*dy;
    lo[2] = 2.*dx*dz;
    lo[3] = dy*dy;
    lo[4] = 2.*dy*dz;
    lo[5] = dz*dz;

    ln[0] = dx*dx;
    ln[1] = 2.*dx*dy;
    ln[2] = 2.*dx*dz;
    ln[3] = dy*dy;
    ln[4] = 2.*dy*dz;
    ln[5] = dz*dz;
  }

  // for now, only neighbouring cell contribution is accounted for
  // TODO: boundary data

  if (!M->pnei_coords) {
    M->pnei_coords = (E_Float **)XCALLOC(M->nppatches, sizeof(E_Float *));
    for (E_Int i = 0; i < M->nppatches; i++) {
      M->pnei_coords[i] = (E_Float *)XCALLOC(
        3*M->ppatches[i].nfaces, sizeof(E_Float));
    }
    comm_interface_data_d(M, &M->cc[0], 3, M->pnei_coords);
  }

  // proc faces
  E_Int npfaces, *pfaces;
  for (E_Int i = 0; i < M->nppatches; i++) {
    npfaces = M->ppatches[i].nfaces;
    pfaces = M->ppatches[i].faces;
    E_Float *pnei_coords = M->pnei_coords[i];
    for (E_Int j = 0; j < npfaces; j++) {
      E_Int own = owner[pfaces[j]];

      E_Float *co = &M->cc[3*own];
      E_Float *cn = &pnei_coords[3*j];

      E_Float dx = cn[0] - co[0];
      E_Float dy = cn[1] - co[1];
      E_Float dz = cn[2] - co[2];

      E_Float *lo = &lsqH[stride*own + 6*count_neis[own]++];
      assert(count_neis[own] <= MAXNEI);

      lo[0] = dx*dx;
      lo[1] = 2.*dx*dy;
      lo[2] = 2.*dx*dz;
      lo[3] = dy*dy;
      lo[4] = 2.*dy*dz;
      lo[5] = dz*dz;
    }
  }

  // compute tA.A
  E_Float *pH, *pHH;
  for (E_Int i = 0; i < M->ncells; i++) {
    pHH = &lsqHH[36*i];
    pH = &lsqH[stride*i];
    E_Int idx = 0;
    for (E_Int j = 0; j < 6; j++) {
      for (E_Int k = 0; k < 6; k++) {
        pHH[idx] = 0.;
        for (E_Int l = 0; l < count_neis[i]; l++)
          pHH[idx] += pH[l*6+j] * pH[l*6+k];
        idx++;
      }
    }
  }

  XFREE(count_neis);
}

E_Float *compute_hessian(mesh* M, E_Float *fld)
{
  compute_cell_centers(M);
  E_Float *G = compute_grad(M, fld);
  compute_lsq_hess_matrices(M);
  
  E_Float *B = (E_Float *)XCALLOC(6, sizeof(E_Float));
  E_Float *H = (E_Float *)XCALLOC(6*M->ncells, sizeof(E_Float));
  E_Float *b = (E_Float *)XCALLOC(MAXNEI*M->ncells, sizeof(E_Float));

  E_Int *count_neis = (E_Int *)XCALLOC(M->ncells, sizeof(E_Int));

  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  // construct b vector
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int nei = neigh[i];
    if (nei < 0) continue;

    E_Int own = owner[i];
  
    E_Float *co = &M->cc[3*own];
    E_Float *cn = &M->cc[3*nei];

    E_Float *go = &G[3*own];
    E_Float *gn = &G[3*nei];

    E_Float dfld = fld[nei] - fld[own];
    E_Float dx = cn[0] - co[0];
    E_Float dy = cn[1] - co[1];
    E_Float dz = cn[2] - co[2];

    b[MAXNEI*own + count_neis[own]++] =
      2.*( dfld - (go[0]*dx + go[1]*dy + go[2]*dz));
    b[MAXNEI*nei + count_neis[nei]++] =
      2.*(-dfld + (gn[0]*dx + gn[1]*dy + gn[2]*dz));
  }

  assert(M->pnei_flds);

  for (E_Int i = 0; i < M->nppatches; i++) {
    //E_Float *pnei_grad  = M->pnei_grads[i];
    E_Float *pnei_coord = M->pnei_coords[i];
    E_Float *pnei_fld   = M->pnei_flds[i];

    for (E_Int j = 0; j < M->ppatches[i].nfaces; j++) {
      E_Int own = owner[M->ppatches[i].faces[j]];

      E_Float *co = &M->cc[3*own];
      E_Float *cn = &pnei_coord[3*j];

      E_Float *go = &G[3*own];
      //E_Float *gn = &pnei_grad[3*j];

      E_Float dfld = pnei_fld[j] - fld[own];
      E_Float dx = cn[0] - co[0];
      E_Float dy = cn[1] - co[1];
      E_Float dz = cn[2] - co[2];

      b[MAXNEI*own + count_neis[own]++] =
        2.*(dfld - (go[0]*dx + go[1]*dy + go[2]*dz));
    }
  }

  // TODO(Imad): boundary data
  
  E_Int stride = 6*MAXNEI;
  for (E_Int i = 0; i < M->ncells; i++) {
    // construct B vector
    E_Float *pH = &M->lsqH[stride*i];
    E_Float *pb = &b[MAXNEI*i];
    for (E_Int j = 0; j < 6; j++) {
      B[j] = 0.;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pH[k*6+j];
    }

    // solve
    E_Int converged = BiCGStab(&M->lsqHH[36*i], &H[6*i], B, 6);
    if (!converged) {
      printf("BiCGStab not converged\n");
      printf("Cell: " SF_D_ "\n", i);
      E_Float *pt = &M->lsqHH[36*i];
      for (E_Int j = 0; j < 6; j++) {
        for (E_Int k = 0; k < 6; k++) {
          printf("%.3e ", pt[6*j+k]);
        }
        puts("");
      }
      for (E_Int j = 0; j < 6; j++) {
        printf("%.3e ", B[j]);
      }
      puts("");
      exit(1);
    }
    //assert(converged);
  }

  XFREE(B);
  XFREE(b);
  XFREE(count_neis);
  XFREE(G);

  return H;
}

void symmat_dot_vec(const E_Float *a, const E_Float *b, E_Float *c)
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[1]*b[0] + a[3]*b[1] + a[4]*b[2];
  c[2] = a[2]*b[0] + a[4]*b[1] + a[5]*b[2];
}

static
E_Float sign(E_Float a)
{
  if (a > 0.) return 1.;
  else if (a < 0.) return -1.;
  else return a;
}

static
void symmat_dot_symmat(E_Float *a, E_Float *b, E_Float *c)
{
  c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  c[1] = a[0]*b[1] + a[1]*b[3] + a[2]*b[4];
  c[2] = a[0]*b[2] + a[1]*b[4] + a[2]*b[5]; 
  c[3] = a[1]*b[0] + a[3]*b[1] + a[4]*b[2];
  c[4] = a[1]*b[1] + a[3]*b[3] + a[4]*b[4];
  c[5] = a[1]*b[2] + a[3]*b[4] + a[4]*b[5]; 
  c[6] = a[2]*b[0] + a[4]*b[1] + a[5]*b[2];
  c[7] = a[2]*b[1] + a[4]*b[3] + a[5]*b[4];
  c[8] = a[2]*b[2] + a[4]*b[4] + a[5]*b[5];
}

static
E_Float symmat_det(E_Float A[6])
{
  E_Float a = A[0];
  E_Float b = A[1];
  E_Float c = A[2];
  E_Float d = A[3];
  E_Float e = A[4];
  E_Float f = A[5]; 
  return a*d*f - (a*e*e + d*c*c + f*b*b) + 2.*b*c*e;
}

static inline
E_Float symmat_trace(E_Float *A)
{
  return A[0] + A[3] + A[5];
}

static
E_Float symmat_second_invariant(E_Float *A)
{
  E_Float AA[9];
  symmat_dot_symmat(A, A, AA);
  return 0.5*(AA[0] + AA[4] + AA[8]);
}

static
E_Float symmat_third_invariant(E_Float *A)
{
  return symmat_det(A);
}

void eigen(E_Float *M, E_Float *L, E_Float *v1, E_Float *v2, E_Float *v3)
{
  const E_Float TOL = 1e-12;
  
  // Init
  v1[0] = 1.; v1[1] = 0.; v1[2] = 0.;
  v2[0] = 0.; v2[1] = 1.; v2[2] = 0.;
  v3[0] = 0.; v3[1] = 0.; v3[2] = 1.;
  L[0] = M[0];
  L[1] = M[3];
  L[2] = M[5];
  
  E_Float Ap[6], B[6], Ar[4];
  E_Float s1[3], s2[3], t2[3], t3[3], r1[3], r2[3], r3[3], tmp1[3], tmp2[3],
    u1[3], u2[3], w1[3];
  E_Float J2, J3, alpha, thirdTrA, norm1, norm2, norm3, coeff, dif, sum, sgn;
  
  // maxm
  E_Float maxm = fabs(M[0]);
  E_Float valm;
  for (E_Int i = 1; i < 6; i++) {
    valm = fabs(M[i]);
    if (valm > maxm) maxm = valm;
  }
  if (maxm < TOL) return;
  
  // normalize matrix
  E_Float dd = 1. / maxm;
  E_Float A[6];
  memcpy(A, M, 6*sizeof(E_Float));
  for (E_Int i = 0; i < 6; i++) A[i] *= dd;
  
  // check for diagonal matrix
  E_Float maxd = fabs(A[1]);
  valm = fabs(A[2]);
  if (valm > maxd) maxd = valm;
  valm = fabs(A[4]);
  if (valm > maxd) maxd = valm;
  if (maxd < TOL) return; // off-diagonal coeffs are smaller than tol
  thirdTrA = (A[0] + A[3] + A[5]) / 3.;
  Ap[0] = A[0] - thirdTrA;
  Ap[1] = A[1];
  Ap[2] = A[2];
  Ap[3] = A[3] - thirdTrA;
  Ap[4] = A[4];
  Ap[5] = A[5] - thirdTrA;
  J2 = symmat_second_invariant(Ap);
  J3 = symmat_third_invariant(Ap);
  E_Float tmp = 0.5*J3 * pow(3./J2, 1.5);
  if (tmp > 1.) tmp = 1.;
  else if (tmp < -1) tmp = -1.;
  alpha = acos(tmp) / 3.;
  
  if (alpha < M_PI/6.) {
    // find L[0] first
    L[0] = 2.*sqrt(J2/3.)*cos(alpha);
  } else {
    // find L[2] first
    L[0] = 2.*sqrt(J2/3.)*cos(alpha + 4.*M_PI/3.);
  }
  
  // find eigenvector corresponding to L[0]
  B[0] = Ap[0] - L[0];
  B[1] = Ap[1];
  B[2] = Ap[2];
  B[3] = Ap[3] - L[0];
  B[4] = Ap[4];
  B[5] = Ap[5] - L[0];
  r1[0] = B[0]; r1[1] = B[1]; r1[2] = B[2];
  r2[0] = B[1]; r2[1] = B[3]; r2[2] = B[4];
  r3[0] = B[2]; r3[1] = B[4]; r3[2] = B[5];
  norm1 = norm(r1, 3);
  norm2 = norm(r2, 3);
  norm3 = norm(r3, 3);
  E_Float over_norm;
  if (norm1 >= norm2 && norm1 >= norm3) {
    over_norm = 1. / norm1;
    s1[0] = r1[0] * over_norm;
    s1[1] = r1[1] * over_norm;
    s1[2] = r1[2] * over_norm;

    coeff = dot(s1, r2, 3);
    t2[0] = r2[0] - coeff*s1[0];
    t2[1] = r2[1] - coeff*s1[1];
    t2[2] = r2[2] - coeff*s1[2];

    coeff = dot(s1, r3, 3);
    t3[0] = r3[0] - coeff*s1[0];
    t3[1] = r3[1] - coeff*s1[1];
    t3[2] = r3[2] - coeff*s1[2];
  } else if (norm2 >= norm1 && norm2 >= norm3) {
    over_norm = 1. / norm2;
    s1[0] = r2[0] * over_norm;
    s1[1] = r2[1] * over_norm;
    s1[2] = r2[2] * over_norm;

    coeff = dot(s1, r1, 3);
    t2[0] = r1[0] - coeff*s1[0];
    t2[1] = r1[1] - coeff*s1[1];
    t2[2] = r1[2] - coeff*s1[2];
    coeff = dot(s1, r3, 3);
    t3[0] = r3[0] - coeff*s1[0];
    t3[1] = r3[1] - coeff*s1[1];
    t3[2] = r3[2] - coeff*s1[2];
  } else {
    over_norm = 1. / norm3;
    s1[0] = r3[0] * over_norm;
    s1[1] = r3[1] * over_norm;
    s1[2] = r3[2] * over_norm;

    coeff = dot(s1, r2, 3);
    t2[0] = r2[0] - coeff*s1[0];
    t2[1] = r2[1] - coeff*s1[1];
    t2[2] = r2[2] - coeff*s1[2];
    coeff = dot(s1, r1, 3);
    t3[0] = r1[0] - coeff*s1[0];
    t3[1] = r1[1] - coeff*s1[1];
    t3[2] = r1[2] - coeff*s1[2];
  }
  
  norm2 = norm(t2, 3);
  norm3 = norm(t3, 3);
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

  cross(s1, s2, v1); // got it!
  
  // Reduced form of Ap
  symmat_dot_vec(Ap, s1, tmp1);
  symmat_dot_vec(Ap, s2, tmp2);
  Ar[0] = dot(s1, tmp1, 3);
  Ar[1] = dot(s1, tmp2, 3);
  Ar[2] = dot(s2, tmp1, 3);
  Ar[3] = dot(s2, tmp2, 3);
  
  // Wilkinson shift
  dif = Ar[0] - Ar[3];
  sum = Ar[0] + Ar[3];
  sgn = sign(dif);
  L[1] = 0.5*sum - 0.5*sgn*sqrt(dif*dif + 4.*Ar[1]*Ar[2]);
  L[2] = sum - L[1];
  
  // find eigenvector corresponding to L1
  B[0] = Ap[0] - L[1];
  B[1] = Ap[1];
  B[2] = Ap[2];
  B[3] = Ap[3] - L[1];
  B[4] = Ap[4];
  B[5] = Ap[5] - L[1];
  symmat_dot_vec(B, s1, u1);
  symmat_dot_vec(B, s2, u2);
  
  norm1 = norm(u1, 3);
  norm2 = norm(u2, 3);
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

  cross(w1, v1, v2);
  cross(v1, v2, v3);
  for (E_Int i = 0; i < 3; i++)
    L[i] = (L[i] + thirdTrA) * maxm;
}

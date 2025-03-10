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
#ifndef NUGA_DIAG_H
#define NUGA_DIAG_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

namespace NUGA
{

#define EPSD      1.e-30
#define EPS       1.e-06

/* seeking 1.e-05 accuracy */
#define  EIGENV_EPSD           1.e-13
#define  EIGENV_EPSD2          1.e-10
#define  EIGENV_EPS6           5.e-06
#define  EIGENV_EPS            1.e-06
#define  EIGENV_EPSX2          2.e-06
#define  MAXITER         50

/**
 * \def EGAL(x,y)
 * Check if numbers \a x and \a y are equal.
 */
#define EGAL(x,y)   (                                             \
    (  ((x) == 0.0f) ? (fabs(y) < EIGENV_EPS) :                      \
       ( ((y) == 0.0f) ? (fabs(x) < EIGENV_EPS) :                    \
         (fabs((x)-(y)) / (fabs(x) + fabs(y)) < EIGENV_EPSX2) )  ) )

/**
 * \brief Identity matrix.
 */
static double Id[3][3] = {
  {1.0, 0.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.0, 0.0, 1.0} };


/**
 * \fn static int newton3(double p[4],double x[3])
 * \brief Find root(s) of a polynomial of degree 3.
 * \param p polynomial coefficients (b=p[2], c=p[1], d=p[0]).
 * \param x root(s) of polynomial.
 * \return 0 if no roots.
 * \return 1 for 3 roots.
 * \return 2 for 2 roots.
 * \return 3 for 1 root.
 *
 * Find root(s) of a polynomial of degree 3: \f$P(x) = x^3+bx^2+cx+d\f$.
 *
 */
static int newton3(double p[4], double x[3]) 
{
  double      b,c,d,da,db,dc,epsd;
  double      delta,fx,dfx,dxx;
  double      fdx0,fdx1,dx0,dx1,x1,x2;
  int         it,n;
  static char mmgWarn=0;

  /* coeffs polynomial, a=1 */
  if ( p[3] != 1. ) {
    if ( !mmgWarn ) {
      fprintf(stderr,"\n  ## Warning: %s: bad use of newton3 function, polynomial"
              " must be of type P(x) = x^3+bx^2+cx+d.\n",
              __func__);
      mmgWarn = 1;
    }
    return 0;
  }

  b = p[2];
  c = p[1];
  d = p[0];

  /* 1st derivative of f */
  da = 3.0;
  db = 2.0*b;

  /* solve 2nd order eqn */
  delta = db*db - 4.0*da*c;
  epsd  = db*db*EIGENV_EPSD2;

  /* inflexion (f'(x)=0, x=-b/2a) */
  x1 = -db / 6.0f;

  n = 1;
  if ( delta > epsd ) {
    delta = sqrt(delta);
    dx0   = (-db + delta) / 6.0;
    dx1   = (-db - delta) / 6.0;
    /* Horner */
    fdx0 = d + dx0*(c+dx0*(b+dx0));
    fdx1 = d + dx1*(c+dx1*(b+dx1));

    if ( fabs(fdx0) < EIGENV_EPSD ) {
      /* dx0: double root, compute single root */
      n = 2;
      x[0] = dx0;
      x[1] = dx0;
      x[2] = -b - 2.0*dx0;
      /* check if P(x) = 0 */
      fx = d + x[2]*(c+x[2]*(b+x[2]));
      if ( fabs(fx) > EIGENV_EPSD2 ) {
#ifdef DEBUG
         fprintf(stderr,"\n  ## Error: %s: ERR 9100, newton3: fx= %E.\n",
                 __func__,fx);
#endif
        return 0;
      }
      return n;
    }
    else if ( fabs(fdx1) < EIGENV_EPSD ) {
      /* dx1: double root, compute single root */
      n = 2;
      x[0] = dx1;
      x[1] = dx1;
      x[2] = -b - 2.0*dx1;
      /* check if P(x) = 0 */
      fx = d + x[2]*(c+x[2]*(b+x[2]));
      if ( fabs(fx) > EIGENV_EPSD2 ) {
#ifdef DEBUG
        fprintf(stderr,"\n  ## Error: %s: ERR 9100, newton3: fx= %E.\n",
                __func__,fx);
#endif
        return 0;
      }
      return n;
    }
  }

  else if ( fabs(delta) < epsd ) {
    /* triple root */
    n = 3;
    x[0] = x1;
    x[1] = x1;
    x[2] = x1;
    /* check if P(x) = 0 */
    fx = d + x[0]*(c+x[0]*(b+x[0]));
    if ( fabs(fx) > EIGENV_EPSD2 ) {
#ifdef DEBUG
      fprintf(stderr,"\n  ## Error: %s: ERR 9100, newton3: fx= %E.\n",
              __func__,fx);
#endif
      return 0;
    }
    return n;
  }

  else {
#ifdef DEBUG
    fprintf(stderr,"\n  ## Error: %s: ERR 9101, newton3: no real roots.\n",
            __func__);
#endif
    return 0;
  }

  /* Newton method: find one root (middle)
     starting point: P"(x)=0 */
  x1  = -b / 3.0;
  dfx =  c + b*x1;
  fx  = d + x1*(c -2.0*x1*x1);
  it  = 0;
  do {
    x2 = x1 - fx / dfx;
    fx = d + x2*(c+x2*(b+x2));
    if ( fabs(fx) < EIGENV_EPSD ) {
      x[0] = x2;
      break;
    }
    dfx = c + x2*(db + da*x2);

    /* check for break-off condition */
    dxx = fabs((x2-x1) / x2);
    if ( dxx < 1.0e-10 ) {
      x[0] = x2;
      if ( fabs(fx) > EIGENV_EPSD2 ) {
        fprintf(stderr,"\n  ## Error: %s: ERR 9102, newton3, no root found"
                " (fx %E).\n",
                __func__,fx);
        return 0;
      }
      break;
    }
    else
      x1 = x2;
  }
  while ( ++it < MAXITER );

  if ( it == MAXITER ) {
    x[0] = x1;
    fx   = d + x1*(c+(x1*(b+x1)));
    if ( fabs(fx) > EIGENV_EPSD2 ) {
      fprintf(stderr,"\n  ## Error: %s: ERR 9102, newton3, no root found"
              " (fx %E).\n",
              __func__,fx);
      return 0;
    }
  }

  /* solve 2nd order equation
     P(x) = (x-sol(1))* (x^2+bb*x+cc)  */
  db    = b + x[0];
  dc    = c + x[0]*db;
  delta = db*db - 4.0*dc;

  if ( delta <= 0.0 ) {
    fprintf(stderr,"\n  ## Error: %s: ERR 9103, newton3, det = 0.\n",__func__);
    return 0;
  }

  delta = sqrt(delta);
  x[1] = 0.5 * (-db+delta);
  x[2] = 0.5 * (-db-delta);

#ifdef DEBUG
  /* check for root accuracy */
  fx = d + x[1]*(c+x[1]*(b+x[1]));
  if ( fabs(fx) > EIGENV_EPSD2 ) {
    fprintf(stderr,"\n  ## Error: %s: ERR 9104, newton3: fx= %E  x= %E.\n",
            __func__,fx,x[1]);
    return 0;
  }
  fx = d + x[2]*(c+x[2]*(b+x[2]));
  if ( fabs(fx) > EIGENV_EPSD2 ) {
    fprintf(stderr,"\n  ## Error: %s: ERR 9104, newton3: fx= %E  x= %E.\n",
            __func__,fx,x[2]);
    return 0;
  }
#endif

  return n;
}

/**
 * \param mat pointer toward a 3x3 matrix.
 * \param lambda eigenvalues.
 * \param v eigenvectors.
 * \param w1 temporary array to perform the matrix cross product.
 * \param w2 temporary array to perform the matrix cross product.
 * \param w3 temporary array to perform the matrix cross product.
 * \param maxm maximal value of the matrix used for normalization.
 * \param order order of eigenvalues (1,2,3) or 0 if failed.
 * \param symmat 0 if matrix is not symetric, 1 otherwise.
 *
 * \return 1 if success, 0 if fail.
 *
 * Check the accuracy of the eigenvalues and vectors computation of a 3x3 matrix
 * (symetric).
 *
 */
static
int check_accuracy(double mat[6],double lambda[3], double v[3][3],
                        double w1[3], double w2[3], double w3[3],
                        double maxm, int order, int symmat) {
  double  err,tmpx,tmpy,tmpz;
  float   m[6];
  int     i,j,k;

  if ( !symmat ) return 1;

  k = 0;
  for (i=0; i<3; i++) {
    for (j=i; j<3; j++) {
      m[k++] = lambda[0]*v[i][0]*v[j][0]
        + lambda[1]*v[i][1]*v[j][1]
        + lambda[2]*v[i][2]*v[j][2];
    }
  }
  err = fabs(mat[0]-m[0]);
  for (i=1; i<6; i++)
    if ( fabs(m[i]-mat[i]) > err )  err = fabs(m[i]-mat[i]);

  if ( err > 1.e03*maxm ) {
    fprintf(stderr,"\n  ## Error: %s:\nProbleme eigenv3: err= %f\n",__func__,err*maxm);
    fprintf(stderr,"\n  ## Error: %s:mat depart :\n",__func__);
    fprintf(stderr,"\n  ## Error: %s:%13.6f  %13.6f  %13.6f\n",__func__,mat[0],mat[1],mat[2]);
    fprintf(stderr,"\n  ## Error: %s:%13.6f  %13.6f  %13.6f\n",__func__,mat[1],mat[3],mat[4]);
    fprintf(stderr,"\n  ## Error: %s:%13.6f  %13.6f  %13.6f\n",__func__,mat[2],mat[4],mat[5]);
    fprintf(stderr,"\n  ## Error: %s:mat finale :\n",__func__);
    fprintf(stderr,"\n  ## Error: %s:%13.6f  %13.6f  %13.6f\n",__func__,m[0],m[1],m[2]);
    fprintf(stderr,"\n  ## Error: %s:%13.6f  %13.6f  %13.6f\n",__func__,m[1],m[3],m[4]);
    fprintf(stderr,"\n  ## Error: %s:%13.6f  %13.6f  %13.6f\n",__func__,m[2],m[4],m[5]);
    fprintf(stderr,"\n  ## Error: %s:lambda : %f %f %f\n",__func__,lambda[0],lambda[1],lambda[2]);
    fprintf(stderr,"\n  ## Error: %s: ordre %d\n",__func__,order);
    fprintf(stderr,"\n  ## Error: %s:\nOrtho:\n",__func__);
    fprintf(stderr,"\n  ## Error: %s:v1.v2 = %.14f\n",__func__,
            v[0][0]*v[1][0]+v[0][1]*v[1][1]+ v[0][2]*v[1][2]);
    fprintf(stderr,"\n  ## Error: %s:v1.v3 = %.14f\n",__func__,
            v[0][0]*v[2][0]+v[0][1]*v[2][1]+ v[0][2]*v[2][2]);
    fprintf(stderr,"\n  ## Error: %s:v2.v3 = %.14f\n",__func__,
            v[1][0]*v[2][0]+v[1][1]*v[2][1]+ v[1][2]*v[2][2]);

    fprintf(stderr,"\n  ## Error: %s:Consistency\n",__func__);
    for (i=0; i<3; i++) {
      tmpx = v[0][i]*m[0] + v[1][i]*m[1]
        + v[2][i]*m[2] - lambda[i]*v[0][i];
      tmpy = v[0][i]*m[1] + v[1][i]*m[3]
        + v[2][i]*m[4] - lambda[i]*v[1][i];
      tmpz = v[0][i]*m[2] + v[1][i]*m[4]
        + v[2][i]*m[5] - lambda[i]*v[2][i];
      fprintf(stderr,"\n  ## Error: %s: Av %d - lambda %d *v %d = %f %f %f\n",
              __func__,i,i,i,tmpx,tmpy,tmpz);

      fprintf(stderr,"\n  ## Error: %s:w1 %f %f %f\n",__func__,w1[0],w1[1],w1[2]);
      fprintf(stderr,"\n  ## Error: %s:w2 %f %f %f\n",__func__,w2[0],w2[1],w2[2]);
      fprintf(stderr,"\n  ## Error: %s:w3 %f %f %f\n",__func__,w3[0],w3[1],w3[2]);
    }
    return 0;
  }

  return 1;
}

/**
 * \brief Find eigenvalues and vectors of a 3x3 matrix.
 * \param symmat 0 if matrix is not symetric, 1 otherwise.
 * \param mat pointer toward the matrix.
 * \param lambda eigenvalues.
 * \param v eigenvectors.
 *
 * \return order of eigenvalues (1,2,3) or 0 if failed.
 *
 * \remark the i^{th} eigenvector is stored in v[i][.].
 *
 */
inline int eigenv(int symmat,double *mat,double lambda[3],double v[3][3]) 
{
  double    a11,a12,a13,a21,a22,a23,a31,a32,a33;
  double    aa,bb,cc,dd,ee,ii,vx1[3],vx2[3],vx3[3],dd1,dd2,dd3;
  double    maxd,maxm,valm,p[4],w1[3],w2[3],w3[3];
  int       k,n;

  /* default */
  memcpy(v,Id,9*sizeof(double));
  if ( symmat ) {
    lambda[0] = (double)mat[0];
    lambda[1] = (double)mat[3];
    lambda[2] = (double)mat[5];

    maxm = fabs(mat[0]);
    for (k=1; k<6; k++) {
      valm = fabs(mat[k]);
      if ( valm > maxm )  maxm = valm;
    }
    /* single float accuracy */
    if ( maxm < EIGENV_EPS6 )  return 1;

    /* normalize matrix */
    dd  = 1.0 / maxm;
    a11 = mat[0] * dd;
    a12 = mat[1] * dd;
    a13 = mat[2] * dd;
    a22 = mat[3] * dd;
    a23 = mat[4] * dd;
    a33 = mat[5] * dd;

    /* diagonal matrix */
    maxd = fabs(a12);
    valm = fabs(a13);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a23);
    if ( valm > maxd )  maxd = valm;
    if ( maxd < EIGENV_EPSD )  return 1;

    a21  = a12;
    a31  = a13;
    a32  = a23;

    /* build characteristic polynomial
       P(X) = X^3 - trace X^2 + (somme des mineurs)X - det = 0 */
    aa = a11*a22;
    bb = a23*a32;
    cc = a12*a21;
    dd = a13*a31;
    p[0] =  a11*bb + a33*(cc-aa) + a22*dd -2.0*a12*a13*a23;
    p[1] =  a11*(a22 + a33) + a22*a33 - bb - cc - dd;
    p[2] = -a11 - a22 - a33;
    p[3] =  1.0;
  }
  else {
    lambda[0] = (double)mat[0];
    lambda[1] = (double)mat[4];
    lambda[2] = (double)mat[8];

    maxm = fabs(mat[0]);
    for (k=1; k<9; k++) {
      valm = fabs(mat[k]);
      if ( valm > maxm )  maxm = valm;
    }
    if ( maxm < EIGENV_EPS6 )  return 1;

    /* normalize matrix */
    dd  = 1.0 / maxm;
    a11 = mat[0] * dd;
    a12 = mat[1] * dd;
    a13 = mat[2] * dd;
    a21 = mat[3] * dd;
    a22 = mat[4] * dd;
    a23 = mat[5] * dd;
    a31 = mat[6] * dd;
    a32 = mat[7] * dd;
    a33 = mat[8] * dd;

    /* diagonal matrix */
    maxd = fabs(a12);
    valm = fabs(a13);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a23);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a21);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a31);
    if ( valm > maxd )  maxd = valm;
    valm = fabs(a32);
    if ( valm > maxd )  maxd = valm;
    if ( maxd < EIGENV_EPSD )  return 1;

    /* build characteristic polynomial
       P(X) = X^3 - trace X^2 + (somme des mineurs)X - det = 0 */
    aa = a22*a33 - a23*a32;
    bb = a23*a31 - a21*a33;
    cc = a21*a32 - a31*a22;
    ee = a11*a33 - a13*a31;
    ii = a11*a22 - a12*a21;

    p[0] =  -a11*aa - a12*bb - a13*cc;
    p[1] =  aa + ee + ii;
    p[2] = -a11 - a22 - a33;
    p[3] =  1.0;
  }

  /* solve polynomial (find roots using newton) */
  n = newton3(p,lambda);
  if ( n <= 0 )  return 0;

  /* compute eigenvectors:
     an eigenvalue belong to orthogonal of Im(A-lambda*Id) */
  v[0][0] = 1.0; v[0][1] = v[0][2] = 0.0;
  v[1][1] = 1.0; v[1][0] = v[1][2] = 0.0;
  v[2][2] = 1.0; v[2][0] = v[2][1] = 0.0;

  w1[1] = a12;  w1[2] = a13;
  w2[0] = a21;  w2[2] = a23;
  w3[0] = a31;  w3[1] = a32;

  if ( n == 1 ) {
    /* vk = crsprd(wi,wj) */
    for (k=0; k<3; k++) {
      w1[0] = a11 - lambda[k];
      w2[1] = a22 - lambda[k];
      w3[2] = a33 - lambda[k];

      /* cross product vectors in (Im(A-lambda(i) Id) ortho */
      vx1[0] = w1[1]*w3[2] - w1[2]*w3[1];
      vx1[1] = w1[2]*w3[0] - w1[0]*w3[2];
      vx1[2] = w1[0]*w3[1] - w1[1]*w3[0];
      dd1    = vx1[0]*vx1[0] + vx1[1]*vx1[1] + vx1[2]*vx1[2];

      vx2[0] = w1[1]*w2[2] - w1[2]*w2[1];
      vx2[1] = w1[2]*w2[0] - w1[0]*w2[2];
      vx2[2] = w1[0]*w2[1] - w1[1]*w2[0];
      dd2    = vx2[0]*vx2[0] + vx2[1]*vx2[1] + vx2[2]*vx2[2];

      vx3[0] = w2[1]*w3[2] - w2[2]*w3[1];
      vx3[1] = w2[2]*w3[0] - w2[0]*w3[2];
      vx3[2] = w2[0]*w3[1] - w2[1]*w3[0];
      dd3    = vx3[0]*vx3[0] + vx3[1]*vx3[1] + vx3[2]*vx3[2];

      /* find vector of max norm */
      if ( dd1 > dd2 ) {
        if ( dd1 > dd3 ) {
          dd1 = 1.0 / sqrt(dd1);
          v[k][0] = vx1[0] * dd1;
          v[k][1] = vx1[1] * dd1;
          v[k][2] = vx1[2] * dd1;
        }
        else {
          dd3 = 1.0 / sqrt(dd3);
          v[k][0] = vx3[0] * dd3;
          v[k][1] = vx3[1] * dd3;
          v[k][2] = vx3[2] * dd3;
        }
      }
      else {
        if ( dd2 > dd3 ) {
          dd2 = 1.0 / sqrt(dd2);
          v[k][0] = vx2[0] * dd2;
          v[k][1] = vx2[1] * dd2;
          v[k][2] = vx2[2] * dd2;
        }
        else {
          dd3 = 1.0 / sqrt(dd3);
          v[k][0] = vx3[0] * dd3;
          v[k][1] = vx3[1] * dd3;
          v[k][2] = vx3[2] * dd3;
        }
      }
    }
  }

  /* (vp1,vp2) double,  vp3 simple root */
  else if ( n == 2 ) {
    w1[0] = a11 - lambda[2];
    w2[1] = a22 - lambda[2];
    w3[2] = a33 - lambda[2];

    /* cross product */
    vx1[0] = w1[1]*w3[2] - w1[2]*w3[1];
    vx1[1] = w1[2]*w3[0] - w1[0]*w3[2];
    vx1[2] = w1[0]*w3[1] - w1[1]*w3[0];
    dd1 = vx1[0]*vx1[0] + vx1[1]*vx1[1] + vx1[2]*vx1[2];

    vx2[0] = w1[1]*w2[2] - w1[2]*w2[1];
    vx2[1] = w1[2]*w2[0] - w1[0]*w2[2];
    vx2[2] = w1[0]*w2[1] - w1[1]*w2[0];
    dd2 = vx2[0]*vx2[0] + vx2[1]*vx2[1] + vx2[2]*vx2[2];

    vx3[0] = w2[1]*w3[2] - w2[2]*w3[1];
    vx3[1] = w2[2]*w3[0] - w2[0]*w3[2];
    vx3[2] = w2[0]*w3[1] - w2[1]*w3[0];
    dd3 = vx3[0]*vx3[0] + vx3[1]*vx3[1] + vx3[2]*vx3[2];

    /* find vector of max norm */
    if ( dd1 > dd2 ) {
      if ( dd1 > dd3 ) {
        dd1 = 1.0 / sqrt(dd1);
        v[2][0] = vx1[0] * dd1;
        v[2][1] = vx1[1] * dd1;
        v[2][2] = vx1[2] * dd1;
      }
      else {
        dd3 = 1.0 / sqrt(dd3);
        v[2][0] = vx3[0] * dd3;
        v[2][1] = vx3[1] * dd3;
        v[2][2] = vx3[2] * dd3;
      }
    }
    else {
      if ( dd2 > dd3 ) {
        dd2 = 1.0 / sqrt(dd2);
        v[2][0] = vx2[0] * dd2;
        v[2][1] = vx2[1] * dd2;
        v[2][2] = vx2[2] * dd2;
      }
      else {
        dd3 = 1.0 / sqrt(dd3);
        v[2][0] = vx3[0] * dd3;
        v[2][1] = vx3[1] * dd3;
        v[2][2] = vx3[2] * dd3;
      }
    }

    /* compute v1 and v2 in Im(A-vp3*Id) */
    dd1 = w1[0]*w1[0] + w1[1]*w1[1] + w1[2]*w1[2];
    dd2 = w2[0]*w2[0] + w2[1]*w2[1] + w2[2]*w2[2];
    if ( dd1 > dd2 ) {
      dd1 = 1.0 / sqrt(dd1);
      v[0][0] = w1[0]*dd1;
      v[0][1] = w1[1]*dd1;
      v[0][2] = w1[2]*dd1;
    }
    else {
      dd2 = 1.0 / sqrt(dd2);
      v[0][0] = w2[0]*dd2;
      v[0][1] = w2[1]*dd2;
      v[0][2] = w2[2]*dd2;
    }

    /* 3rd vector orthogonal */
    v[1][0] = v[2][1]*v[0][2] - v[2][2]*v[0][1];
    v[1][1] = v[2][2]*v[0][0] - v[2][0]*v[0][2];
    v[1][2] = v[2][0]*v[0][1] - v[2][1]*v[0][0];
    dd1 = v[1][0]*v[1][0] + v[1][1]*v[1][1] + v[1][2]*v[1][2];
    dd1 = 1.0 / sqrt(dd1);
    v[1][0] *= dd1;
    v[1][1] *= dd1;
    v[1][2] *= dd1;
  }

  lambda[0] *= maxm;
  lambda[1] *= maxm;
  lambda[2] *= maxm;

  /* check accuracy */
  if ( getenv("MEIGENV_DDEBUG") && symmat ) {
    if ( !check_accuracy ( mat, lambda, v, w1, w2, w3, maxm, n, symmat ) )
      return 0;
  }

  return n;
}

/**
 * \brief Find eigenvalues and vectors of a 2x2 matrix.
 * \param mm pointer toward the matrix.
 * \param lambda pointer toward the output eigenvalues.
 * \param vp eigenvectors.
 * \return 1.
 *
 * \warning not used for now
 */
inline int eigen2(double *mm,double *lambda,double vp[2][2]) {
  double   m[3],dd,a1,xn,ddeltb,rr1,rr2,ux,uy;

  /* normalize */
  memcpy(m,mm,3*sizeof(double));
  xn = fabs(m[0]);
  if ( fabs(m[1]) > xn )  xn = fabs(m[1]);
  if ( fabs(m[2]) > xn )  xn = fabs(m[2]);
  if ( xn < EIGENV_EPSD2 ) {
    lambda[0] = lambda[1] = 0.0;
    vp[0][0] = 1.0;
    vp[0][1] = 0.0;
    vp[1][0] = 0.0;
    vp[1][1] = 1.0;
    return 1;
  }
  xn = 1.0 / xn;
  m[0] *= xn;
  m[1] *= xn;
  m[2] *= xn;

  if ( EGAL(m[1],0.0) ) {
    rr1 = m[0];
    rr2 = m[2];
    goto vect;
  }

  /* eigenvalues of jacobian */
  a1     = -(m[0] + m[2]);
  ddeltb = a1*a1 - 4.0 * (m[0]*m[2] - m[1]*m[1]);

  if ( ddeltb < 0.0 ) {
    fprintf(stderr,"\n  ## Error: %s: Delta: %f\n",__func__,ddeltb);
    ddeltb = 0.0;
  }
  ddeltb = sqrt(ddeltb);

  if ( fabs(a1) < EIGENV_EPS ) {
    rr1 = 0.5 * sqrt(ddeltb);
    rr2 = -rr1;
  }
  else if ( a1 < 0.0 ) {
    rr1 = 0.5 * (-a1 + ddeltb);
    rr2 = (-m[1]*m[1] + m[0]*m[2]) / rr1;
  }
  else if ( a1 > 0.0 ) {
    rr1 = 0.5 * (-a1 - ddeltb);
    rr2 = (-m[1]*m[1] + m[0]*m[2]) / rr1;
  }
  else {
    rr1 = 0.5 * ddeltb;
    rr2 = -rr1;
  }

vect:
  xn = 1.0 / xn;
  lambda[0] = rr1 * xn;
  lambda[1] = rr2 * xn;

  /* eigenvectors */
  a1 = m[0] - rr1;
  if ( fabs(a1)+fabs(m[1]) < EIGENV_EPS ) {
    if (fabs(lambda[1]) < fabs(lambda[0]) ) {
      ux = 1.0;
      uy = 0.0;
    }
    else {
      ux = 0.0;
      uy = 1.0;
    }
  }
  else if ( fabs(a1) < fabs(m[1]) ) {
    ux = 1.0;
    uy = -a1 / m[1];
  }
  else if ( fabs(a1) > fabs(m[1]) ) {
    ux = -m[1] / a1;
    uy = 1.0;
  }
  else if ( fabs(lambda[1]) > fabs(lambda[0]) ) {
    ux = 0.0;
    uy = 1.0;
  }
  else {
    ux = 1.0;
    uy = 0.0;
  }

  dd = sqrt(ux*ux + uy*uy);
  dd = 1.0 / dd;
  if ( fabs(lambda[0]) > fabs(lambda[1]) ) {
    vp[0][0] =  ux * dd;
    vp[0][1] =  uy * dd;
  }
  else {
    vp[0][0] =  uy * dd;
    vp[0][1] = -ux * dd;
  }

  /* orthogonal vector */
  vp[1][0] = -vp[0][1];
  vp[1][1] =  vp[0][0];

  return 1;
}

/**
 * \param m terms of symetric matrix \f$2x2\f$.
 * \param lambda eigenvalues of \a m.
 * \param vp eigenvectors of \a m.
 * \return order of the eigenvalues.
 *
 * Compute eigenelements of a symetric matrix m. Eigenvectors are orthogonal.
 *
 */
inline int eigensym(double m[3],double lambda[2],double vp[2][2]) {
  double   sqDelta,dd,trm,vnorm;

  dd  = m[0]-m[2];
  trm = m[0]+m[2];
  sqDelta = sqrt(dd*dd + 4.0*m[1]*m[1]);
  lambda[0] = 0.5*(trm - sqDelta);

  /* Case when m = lambda[0]*I */
  if ( sqDelta < EPS ) {
    lambda[1] = lambda[0];
    vp[0][0] = 1.0;
    vp[0][1] = 0.0;

    vp[1][0] = 0.0;
    vp[1][1] = 1.0;
    return 2;
  }
  vp[0][0] = m[1];
  vp[0][1] = (lambda[0] - m[0]);
  vnorm = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);

  if ( vnorm < EPS ) {
    vp[0][0] = (lambda[0] - m[2]);
    vp[0][1] = m[1];
    vnorm = sqrt(vp[0][0]*vp[0][0] + vp[0][1]*vp[0][1]);
  }
  assert(vnorm > EPSD);

  vnorm = 1.0/vnorm;
  vp[0][0] *= vnorm;
  vp[0][1] *= vnorm;

  vp[1][0] = -vp[0][1];
  vp[1][1] = vp[0][0];

  lambda[1] = m[0]*vp[1][0]*vp[1][0] + 2.0*m[1]*vp[1][0]*vp[1][1]
    + m[2]*vp[1][1]*vp[1][1];

  return 1;
}

} // namespace

#endif

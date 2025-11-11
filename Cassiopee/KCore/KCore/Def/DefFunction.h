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
#ifndef _KCORE_DEF_DEFFUNCTION_H_
#define _KCORE_DEF_DEFFUNCTION_H_

#include "Def/DefTypes.h"
#include "Def/DefCplusPlusConst.h"
#include <math.h>

// ___________________________________________________
//
// Note: Use f suffix to specify a float constant
//       Use L suffix to specify a long constant
//
// Example: E_max(0.0f, a, b), E_abs(-50.2f)
// ---------------------------------------------------

namespace K_FUNC {

// Factorielle 
inline E_Int fact(E_Int n)
{
  switch (n)
  {
    case 1: return 1;
    case 2: return 2;
    case 3: return 6;
    case 4: return 24;
    case 5: return 120;
    case 6: return 720;
    case 7: return 5040;
    default: return (n)*fact(n-1);
  }
}

///+ Numerical Helper Functions
/** Minimum for 2 E_Float  */
inline E_Float E_min(E_Float a, E_Float b)
{
  return ( (b) < (a) ? (b) : (a) );
}

/** Minimum for 2 E_Int */
inline E_Int E_min(E_Int a, E_Int b)
{
  return ( (b) < (a) ? (b) : (a) );
}

/** Minimum for 2 short */
inline short E_min(short a, short b)
{
  return ( (b) < (a) ? (b) : (a) );
}

/** Maximum for 2 E_Float */
inline E_Float E_max(E_Float a, E_Float b)
{
  return ( (b) > (a) ? (b) : (a) );
}

/** Maximum for 2 E_Int */
inline E_Int E_max(E_Int a, E_Int b)
{
  return ( (b) > (a) ? (b) : (a) );
}

/** Maximum for 2 short */
inline short E_max(short a, short b)
{
  return ( (b) > (a) ? (b) : (a) );
}

/** Minimum for 3 E_Float */
inline E_Float E_min(E_Float a, E_Float b, E_Float c)
{
  return ( E_min( E_min(a,b), c ));
}

/** Minimum for 3 E_Int */
inline E_Int E_min(E_Int a, E_Int b, E_Int c)
{
  return ( E_min( E_min(a,b), c ));
}

/** Maximum for 3 E_Float */
inline E_Float E_max(E_Float a, E_Float b, E_Float c)
{
  return ( E_max( E_max(a,b), c ));
}

/** Maximum for 3 E_Int */
inline E_Int E_max(E_Int a, E_Int b, E_Int c)
{
  return ( E_max( E_max(a,b), c ));
}

/** Maximum for 4 E_Float */
inline E_Float E_max(E_Float a, E_Float b, E_Float c, E_Float d)
{
  return ( E_max( E_max(a,b), E_max(c,d) ));
}

/** Maximum for 4 E_Int */
inline E_Int E_max(E_Int a, E_Int b, E_Int c, E_Int d)
{
  return ( E_max( E_max(a,b), E_max(c,d) ));
}

/** Absolute value for E_Float */
inline E_Float E_abs(E_Float a)
{
  return (a < 0) ? -a : a;
}

#ifndef E_DOUBLEREAL
inline double E_abs(double a)
{
  return (a < 0) ? -a : a;
}
#endif

/** Absolute value for E_Int */ 
inline E_Int E_abs(E_Int a)
{
  return (a < 0) ? -a : a;
}

/** Sign for E_Float (return -1 or 1) */
inline E_Float E_sign(E_Float a)
{
  return (a < 0.0) ? -1.0 : 1.0;
}

/** Sign for E_Int (return -1 or 1) */
inline E_Int E_sign(E_Int a)
{
  return (a < 0) ? -1 : 1;
}

/** Sign for E_Int (-1,0,1) */
inline E_Int E_signum(E_Int a)
{
  return (a > 0)-(a < 0);
}

/** Swap for E_Float */
inline void E_swap(E_Float a, E_Float b)
{
  E_Float tmp = a;
  a = b;
  b = tmp;
}

/** Swap for E_Int */
inline void E_swap(E_Int a, E_Int b)
{
  E_Int tmp = a;
  a = b;
  b = tmp;
}

/** Equality test with given accuracy */
inline
E_Bool fEqual(E_Float lhs, E_Float rhs, 
              E_Float precision=K_CONST::E_CUTOFF)
{
  E_Float t = lhs - rhs;
  return (t >= -precision && t <= precision) ? true : false;
}   

/** Equality to 0 test with given accuracy */
inline
E_Bool fEqualZero(E_Float lhs, E_Float precision=K_CONST::E_CUTOFF)
{
  return (lhs >= -precision && lhs <= precision) ? true : false;
}

/** Square of the distance between to points with dim coordinates of type 
    E_Float stored in two sequence containers (array, vector...). */
template <typename InputIterator1, typename InputIterator2>
inline
E_Float sqrDistance (const InputIterator1 i1, const InputIterator2 i2, 
                     E_Int dim)
{
  InputIterator1 p1(i1);
  InputIterator2 p2(i2);
  E_Float result = 0.;

  for (E_Int i = 0; i < dim; ++i)
  {
    result += (*p2 - *p1)*(*p2 - *p1);
    ++p1;
    ++p2;
  }
  return result;
}

// z = x - y
template <E_Int dim, typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
diff (InputIterator1 x, InputIterator2 y, InputIterator3 z) 
{
  for (E_Int i = 0; i < dim; ++i) 
  {
    *(z+i) = *(x+i) - *(y+i); 
  }
  return z + dim;
}

template <E_Int dim, typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
diff (InputIterator1 x, InputIterator2 y, E_Int stride, InputIterator3 z) 
{
  for (E_Int i = 0; i < dim; ++i)
  {
    *(z+i) = *(x+i*stride) - *(y+i*stride);
  }
  return z + dim;
}

// z = x + y
template <E_Int dim, typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
sum (InputIterator1 x, InputIterator2 y, InputIterator3 z) 
{
  for (E_Int i = 0; i < dim; ++i)
    *(z+i) = *(x+i) + *(y+i);
  return z + dim;
}

// z = ax + by
template <E_Int dim, typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
sum (E_Float a, InputIterator1 x, E_Float b, InputIterator2 y,  InputIterator3 z) 
{
  for (E_Int i = 0; i < dim; ++i)
    *(z+i) = *(x+i)*a + *(y+i)*b;
  return z + dim;
}

// z = ax + by + c
template <E_Int dim, typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
sum (E_Float a, InputIterator1 x, E_Float b, InputIterator2 y,  InputIterator2 c, InputIterator3 z) 
{
  for (E_Int i = 0; i < dim; ++i)
    *(z+i) = *(x+i)*a + *(y+i)*b + *(c+i);
  return z + dim;
}

// z = ax + y
template <E_Int dim, typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
sum (E_Float a, InputIterator1 x, InputIterator2 y, InputIterator3 z) 
{
  for (E_Int i = 0; i < dim; ++i)
    *(z+i) = *(x+i)*a + *(y+i);
  return z + dim;
}

/** Square of the norm of a vector with dim coordinates of type E_Float
 * stored in a sequence container (array, vector...). */
template <E_Int dim, typename InputIterator>
inline
E_Float sqrNorm(InputIterator it)
{
  E_Float result = 0.;
  for (E_Int i = 0; i < dim; ++i)
  { result += (*it) * (*it); it++; }
  return result;
}

template <E_Int dim>
inline
void crossProduct(const E_Float* x, const E_Float* y, E_Float* z);

template <>
inline
void crossProduct<2>(const E_Float* x, const E_Float* y, E_Float* z)
{
  *z = *x * (*(y+1)) - *(x+1) * (*y);
}

template <>
inline
void crossProduct<3> (const E_Float* x, const E_Float* y, E_Float* z) 
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}

// | u1 v1 w1 |
// | u2 v2 w2 |
// | u3 v3 w3 |
//#define zzdet3(u1,u2,u3, v1,v2,v3, w1,w2,w3) (u1*(v2*w3 - v3*w2) + u2*(v3*w1 - v1*w3) + u3*(v1*w2 - v2*w1))
inline
E_Float zzdet3(E_Float u1, E_Float u2, E_Float u3, E_Float v1, E_Float v2, E_Float v3, E_Float w1, E_Float w2, E_Float w3)
{
  return (u1*(v2*w3 - v3*w2) + u2*(v3*w1 - v1*w3) + u3*(v1*w2 - v2*w1));
}
// Produit mixte u.(v ^ w) = det(u,v,w)
inline
E_Float tripleProduct(const E_Float* u, const E_Float* v, E_Float* w)
{
  return zzdet3(u[0],u[1],u[2], v[0],v[1],v[2], w[0],w[1],w[2]);
  //return (u[0]*(v[1]*w[2] - v[2]*w[1]) + u[1]*(v[2]*w[0] - v[0]*w[2]) + u[2]*(v[0]*w[1] - v[1]*w[0]));
}

// the sign of this determinant tells whether Q is over (+) or under (-) plan (P0,P1,P2).
inline
E_Float zzdet4(const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q)
{
  return zzdet3(Q[0]-P0[0], Q[1]-P0[1], Q[2]-P0[2],  P1[0]-P0[0], P1[1]-P0[1], P1[2]-P0[2], P2[0]-P0[0], P2[1]-P0[1], P2[2]-P0[2]);
}

template <E_Int dim>
inline
E_Float dot(const E_Float* x, const E_Float* y);

template <>
inline
E_Float dot<2> (const E_Float* x, const E_Float* y) 
{
  return (*x * (*y)) + (*(x+1) * (*(y+1)));
}

template <>
inline
E_Float dot<3> (const E_Float* x, const E_Float* y) 
{
  return (*x * (*y)) + (*(x+1) * (*(y+1))) + (*(x+2) * (*(y+2)));
}

template <E_Int dim, typename InputIterator>
inline
E_Float normalize (InputIterator it)
{
  E_Float L0 = ::sqrt(sqrNorm<dim>(it));
  if (L0 != 0.)
  {
    E_Float L1 = 1./L0;
    for (E_Int i = 0; i < dim; ++i) *(it+i) *= L1;
  }
  return L0;
}

// Returns (x^y).(x^y)
template <E_Int dim>
inline
E_Float sqrCross(const E_Float* x, const E_Float* y);

template <>
inline
E_Float sqrCross<2>(const E_Float* x, const E_Float* y)
{
  E_Float d;
  K_FUNC::crossProduct<2>(x, y, &d);
  return d*d;
}

template <>
inline
E_Float sqrCross<3>(const E_Float* x, const E_Float* y)
{
  E_Float d[3];
  K_FUNC::crossProduct<3>(x, y, d);
  return K_FUNC::sqrNorm<3>(d);
}

///-

} // K_FUNC
#endif
// ===== Def/DefFunction.h === Last line ===

/*
 
 
 
              NUGA 
 
 
Author : Sam Landier (sam.landier@onera.fr) 
 */

#ifndef NUGA_MATHS_HXX
#define NUGA_MATHS_HXX

#include "Nuga/include/defs.h"
#include "Nuga/include/macros.h"
#include <cmath>

namespace NUGA {

/** Square of the distance between to points with dim coordinates of type 
    E_Float stored in two sequence containers (array, vector...). */
template <typename InputIterator1, typename InputIterator2>
inline
E_Float sqrDistance (const InputIterator1 i1, const InputIterator2 i2,
                     E_Int dim)
{
  InputIterator1  p1(i1);
  InputIterator2  p2(i2);
  E_Float         result = 0.;

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
// Prouit mixte u.(v ^ w) = det(u,v,w)
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
  NUGA::crossProduct<2>(x, y, &d);
  return d*d;
}

template <>
inline
E_Float sqrCross<3>(const E_Float* x, const E_Float* y)
{
  E_Float d[3];
  NUGA::crossProduct<3>(x, y, d);
  return NUGA::sqrNorm<3>(d);
}

///
inline double project(double const * plane_pt, double const * plane_dir, double const * pt, double const * dir, double* proj_pt)
{
  double PPt[3];
  NUGA::diff<3>(pt, plane_pt, PPt);
  double k = -NUGA::dot<3>(PPt, plane_dir);
  double c = NUGA::dot<3>(plane_dir, dir);
  //assert(SIGN(c, EPSILON) != 0); //fixme
  k /= c;
  NUGA::sum<3>(1., pt, k, dir, proj_pt); //project
  return k;
}

///
inline double angle_measure
(const E_Float* ni, const E_Float* nj, const E_Float* E0, const E_Float* E1)
{
  //
  E_Float nk[3];
  NUGA::crossProduct<3>(ni, nj, nk);
  E_Float c = NUGA::dot<3>(ni, nj);

  E_Int s = zSIGN(::fabs(c) - 1., ZERO_M);

  if (s != 0) // non-degn case
  {
    E_Float E0E1[3];
    NUGA::diff<3>(E1, E0, E0E1);
    NUGA::normalize<3>(E0E1);
    E_Float K2 = -NUGA::dot<3>(E0E1, nk);
    E_Int signK2 = zSIGN(K2, ZERO_M);

    E_Float s2 = NUGA::sqrNorm<3>(nk);

#ifdef DEBUG_GEOM_ALGO
    if (signK2 == 0) std::cout << "ERROR : GeomAlgo::angle_measure : inconsistence between s and c" << std::endl;
    assert(signK2 != 0);
#endif

    E_Float alpha = ::atan2(::sqrt(s2), c);
    alpha = NUGA::PI - signK2 * alpha;

    return alpha;
  }
  else // (s == 0) : ni and nj are nearly colinear : 0, Pi or 2Pi
  {
    if (c > 0.) return NUGA::PI;

#ifdef DEBUG_GEOM_ALGO
    std::cout << "ERROR : GeomAlgo::angle_measure : 0 or 2Pi ?" << std::endl;
    assert(false);
#endif
    return 2.*NUGA::PI; //error : either 0 or 2Pi are wrong. return one of them as an arbitray choice.
  }

}

///
inline bool angular_weighted_normal(const double* Pim1, const double* Pi, const double* Pip1, double* n)
{
  double ray1[3], ray2[3];
  NUGA::diff<3>(Pip1, Pi, ray1);
  NUGA::normalize<3>(ray1);
  NUGA::diff<3>(Pim1, Pi, ray2);
  NUGA::normalize<3>(ray2);
    
  NUGA::crossProduct<3>(ray1, ray2, n);    // normal
  NUGA::normalize<3>(n);

  // E1 is such PiE1 is normal to Pim1PiPip1
  double ni[3], nj[3], E1[3];
  NUGA::sum<3>(Pi, n, E1);

  // ni = E0E1 ^ PiPip1
  NUGA::crossProduct<3>(n, ray1, ni);

  // nj = PiPim1 ^ E0E1
  NUGA::crossProduct<3>(ray2, n, nj);

  double alpha = angle_measure(ni, nj, Pi/*E0*/, E1);
  if (alpha == 2.*NUGA::PI) return true;

  n[0] *= alpha;
  n[1] *= alpha;
  n[2] *= alpha;

  return false;
}

} // NUGA
#endif


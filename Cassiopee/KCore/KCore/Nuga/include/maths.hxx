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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_MATHS_HXX
#define NUGA_MATHS_HXX

#include "Nuga/include/defs.h"
#include "Nuga/include/macros.h"
#include "Nuga/include/DynArray.h"

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

inline
void crossProduct2D(const long double* x, const long double* y, long double* z)
{
  *z = *x * (*(y + 1)) - *(x + 1) * (*y);
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

///
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

inline long double sqrCross_2D(const long double* x, const long double* y)
{
  long double d;
  NUGA::crossProduct2D(x, y, &d);
  return d * d;
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
inline double normals_angle (const E_Float* ni, const E_Float* nj)
{
  // Angle between 2 normals (conical tolerance)

  E_Float nk[3];
  NUGA::crossProduct<3>(ni, nj, nk);
  E_Float c = NUGA::dot<3>(ni, nj);

  E_Int s = zSIGN(::fabs(c) - 1., ZERO_M);

  if (s != 0) // non-degn case
  {
    E_Float s2 = NUGA::sqrNorm<3>(nk);

    E_Float alpha = ::atan2(::sqrt(s2), c);
    return alpha;
  }
  else // (s == 0) : ni and nj are nearly colinear : 0, Pi or 2Pi
  {
    if (c > 0.) return 0.;

#ifdef DEBUG_GEOM_ALGO
    std::cout << "ERROR : GeomAlgo::angle_measure : 0 or 2Pi ?" << std::endl;
    assert(false);
#endif
    return NUGA::PI; 
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

///
inline void __get_transform_matrix
(E_Float* U, E_Float*V, E_Float* W, K_FLD::FloatArray& P)
{
  P.resize(3, 3);

  for (E_Int i = 0; i < 3; ++i)
  {
    P(i, 0) = U[i];
    P(i, 1) = V[i];
    P(i, 2) = W[i];
  }
}

///
inline void __get_normal_to
(const E_Float* V, E_Float* N)
{
  N[0] = 1.0;
  for (E_Int k = 1; k < 3; ++k)
    N[k] = 0.;

  if ((V[0] != 0.) || (V[1] != 0.))
  {
    N[0] = -V[1];
    N[1] = V[0];
  }
}

/// Computes a transformation matrix to the coordinate system (having W as 3rd axis)
inline void computeAFrame(const E_Float* W, K_FLD::FloatArray& P)
{
  E_Float U[3], V[3], w[] = { W[0], W[1], W[2] };
  NUGA::normalize<3>(w);
  __get_normal_to(w, U);
  assert(NUGA::dot<3>(w, U) == 0.);
  NUGA::crossProduct<3>(w, U, V);
  NUGA::normalize<3>(U);
  NUGA::normalize<3>(V);

  __get_transform_matrix(U, V, w, P);
}

inline void transform(K_FLD::FloatArray& pos, const K_FLD::FloatArray& t)
{
  K_FLD::FloatArray::iterator pN;
  E_Float Q[3];
  E_Int i, j, n;
  for (i = 0; i < pos.cols(); ++i)
  {
    pN = pos.col(i);

    for (j = 0; j < 3; ++j)
    {
      Q[j] = 0.;
      for (n = 0; n < 3; ++n)
        Q[j] += t(j, n) * (*(pN + n));
    }

    for (j = 0; j < 3; ++j)
      pos(j, i) = Q[j];
  }
}

inline // assume z is (0,0,1)
void computeNodeRadiusAndAngles
(K_FLD::FloatArray& coord, E_Float x0, E_Float y0, std::vector<E_Float>& radius, std::vector<E_Float>& angles)
{
  radius.clear();
  radius.resize(coord.cols(), NUGA::FLOAT_MAX);
  angles.clear();
  angles.resize(coord.cols(), NUGA::FLOAT_MAX);

  for (E_Int i = 0; i < coord.cols(); ++i)
  {
    const E_Float* pt = coord.col(i);

    radius[i] = ::sqrt(((pt[0] - x0)*(pt[0] - x0)) + ((pt[1] - y0)*(pt[1] - y0)));

    E_Float c = (pt[0] - x0) / radius[i];
    E_Float s = (pt[1] - y0) / radius[i];

    angles[i] = ::atan2(s, c);
  }
}

inline void axial_rotate(K_FLD::FloatArray& crd, const double* axis_pt, const double* axis_dir, double angle)
{
  crd.pushBack(axis_pt, axis_pt + 3); // to embark it in the transfo

  K_FLD::FloatArray P(3, 3), iP(3, 3);
  computeAFrame(axis_dir, P);
  iP = P;
  K_FLD::FloatArray::inverse3(iP);
  transform(crd, iP);// Now we are in the reference cylindrical coordinate system.

  double * axi_pt = crd.col(crd.cols() - 1);

  for (E_Int i = 0; i < crd.cols(); ++i)
  {
    double* pt = crd.col(i);
    double X = pt[0] - axi_pt[0];
    double Y = pt[1] - axi_pt[1];
    pt[0] = ::cos(angle) * X - ::sin(angle) * Y + axi_pt[0];
    pt[1] = ::sin(angle) * X + ::cos(angle) * Y + axi_pt[1];
  }

  NUGA::transform(crd, P); // back to original ref frame  

}

template <typename T> inline
T abs(T a) { return ::abs(a); } //wrapper for l.64 V1_smoother.hxx
template <typename T> inline
T max(T a, T b) { return std::max(a, b); } //wrapper for l.361 hierarchical_mesh.hxx
template <typename T> inline
T min(T a, T b) { return std::min(a, b); } //wrapper for l.106 estimator.hxx



inline long szudzik_pairing(int x, int y)
{
  return ((x <y) ? (y * y) + x : (x * x) + x + y);
}

inline void szudzik_unpairing(E_Int szudzic_val, E_Int& x, E_Int& y)
{
  E_Int a = (E_Int)(::sqrt(szudzic_val));
  E_Int a2 = a * a;

  if ((szudzic_val - a2) < a)
  {
    x = szudzic_val - a2;
    y = a;
  }
  else
  {
    x = a;
    y = szudzic_val - a2 - a;
  }
}

template <typename T>
inline T signed_distance2D(const T* P, const T* Q0, const T* Q1)
{
  T normal_Q0Q1[] = { Q0[1] - Q1[1], Q1[0] - Q0[0] };
  normalize<2>(normal_Q0Q1);

  T V01[2];
  diff<2>(P, Q0, V01);

  T pow2D_P0_Q0Q1 = NUGA::dot<2>(V01, normal_Q0Q1);
  return pow2D_P0_Q0Q1;
}


} // NUGA
#endif


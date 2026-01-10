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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef __ANISO_METRIC_TYPE_H__
#define __ANISO_METRIC_TYPE_H__

#include "Nuga/include/defs.h"
#include <array>
#include <cmath>
#include <iostream>

namespace DELAUNAY{

template <short DIM>
class AnisoMetricType
{
public:
  typedef NUGA::size_type            size_type;
  typedef AnisoMetricType<DIM>       self_type;
  #define                            DIMANISO DIM*(DIM+1) / 2

public:

  explicit AnisoMetricType(void)               {for (E_Int k = 0; k < DIMANISO; ++k)_mij[k] = 0.;}
  explicit AnisoMetricType(E_Float h){ _mij[0] = _mij[2] = 1./(h*h);_mij[1] = 0.;} // fixme : 2D only
  explicit AnisoMetricType(const E_Float* m) { _mij[0] = m[0]; _mij[1] = m[1]; _mij[2] = m[2];} // fixme : 2D only

  ~AnisoMetricType(void){}

  inline E_Float operator[](size_type i) const {assert (i < DIMANISO); return _mij[i];}
  inline E_Float& operator[](size_type i) {assert (i < DIMANISO); return _mij[i];}
  inline self_type& operator=(const self_type& m){for (E_Int k = 0; k < DIMANISO; ++k)_mij[k] = m[k]; return *this;}
  inline self_type& operator=(const E_Float* mvals){for (E_Int k = 0; k < DIMANISO; ++k)_mij[k] = mvals[k]; return *this;}

  //fixme imad : not great to have a hard coded test value (should be passed as an argument so not appropriate for operator==)
  inline bool operator==(const self_type& other) const
  {
    for (E_Int k = 0; k < DIMANISO; k++) {
      if (::fabs(_mij[k] - other[k]) > 1e-6) return false;
    }
    return true;
  }

  inline bool operator!=(const self_type& other) const
  {
    return !(*this == other);
  }

  inline self_type operator*(const E_Float& a) const;
  //inline self_type operator*(E_Float* v) const;

  // product of 2 metrics is a 3x3 matrix
  inline std::array<E_Float, 9> operator*(const self_type& N); 


  inline self_type operator/(const E_Float& a) const;
  
  AnisoMetricType operator+(const AnisoMetricType&) const;

 void eigen_values(E_Float &lmax, E_Float & lmin) const;

  inline self_type inverse() const;

  inline E_Float det() const;

#ifndef DEBUG_METRIC
  private:
#else
  public:
#endif
	public:
  E_Float                 _mij[DIM*(DIM+1) / 2]; // 2D -> 3, 3D -> 6

};


using Aniso2D = AnisoMetricType<2>;
using Aniso3D = AnisoMetricType<3>;

///
template <short DIM>
AnisoMetricType<DIM>
AnisoMetricType<DIM>::operator*(const E_Float& a) const {

  self_type res(*this);
  for (E_Int i = 0; i < DIMANISO; ++i)
    res._mij[i] *= a;

  return res;
}

template <> inline
std::array<E_Float, 9>
AnisoMetricType<3>::operator*(const self_type& N)
{
  std::array<E_Float, 9> RES;

  const E_Float& a11 = _mij[0];
  const E_Float& a12 = _mij[1];
  const E_Float& a13 = _mij[2];
  const E_Float& a22 = _mij[3];
  const E_Float& a23 = _mij[4];
  const E_Float& a33 = _mij[5];

  const E_Float& b11 = N[0];
  const E_Float& b12 = N[1];
  const E_Float& b13 = N[2];
  const E_Float& b22 = N[3];
  const E_Float& b23 = N[4];
  const E_Float& b33 = N[5];
	
  RES[0] = a11*b11 + a12*b12 + a13*b13; RES[1] = a11*b12 + a12*b22 + a13*b23; RES[2] = a11*b13 + a12*b23 + a13*b33;
  RES[3] = a12*b11 + a22*b12 + a23*b13; RES[4] = a12*b12 + a22*b22 + a23*b23; RES[5] = a12*b13 + a22*b23 + a23*b33;
  RES[6] = a13*b11 + a23*b12 + a33*b13; RES[7] = a13*b12 + a23*b22 + a33*b23; RES[8] = a13*b13 + a23*b23 + a33*b33;

  return RES;
}



// doesn't check if a != 0
template <short DIM>
AnisoMetricType<DIM>
AnisoMetricType<DIM>::operator/(const E_Float& a) const
{
  self_type res(*this);
  for (E_Int i = 0; i < DIMANISO; ++i)
    res._mij[i] /= a;

  return res;
}

///
template <short DIM>
AnisoMetricType<DIM>
AnisoMetricType<DIM>::operator+(const AnisoMetricType& rhs) const
{
  self_type res(*this);
   for (E_Int i = 0; i < DIMANISO; ++i)
     res._mij[i] += rhs._mij[i];

   return res;
}

template <> inline
void AnisoMetricType<2>::eigen_values(E_Float &lmax, E_Float & lmin) const
{
  E_Float a = _mij[0] + _mij[2];                 //trace
  E_Float b = _mij[0]*_mij[2] - _mij[1]*_mij[1]; //det
  E_Float d = a*a - 4.*b;

  d = (d > 0.) ? sqrt(d) : 0.;
  lmin = 0.5*(a - d);
  lmax = lmin+d;
}

template <> inline
void AnisoMetricType<3>::eigen_values(E_Float &lmax, E_Float & lmin) const
{
  //todo Imad : fixme : asssume here a diagonal matrix!
  lmin = std::min(_mij[0], std::min(_mij[3], _mij[5]));
  lmax = std::max(_mij[0], std::max(_mij[3], _mij[5]));
}

#ifdef DEBUG_METRIC
inline std::ostream&
 operator<<(std::ostream& os, const DELAUNAY::AnisoMetricType<2>& m)  
{  
  os << m._mij[0] << '/' << m._mij[1] << '/' << m._mij[2];  
  return os;  
}  
#endif

inline std::ostream&
 operator<<(std::ostream& os, const DELAUNAY::AnisoMetricType<3>& m)  
{  
  os << m[0] << " " << m[1] << " " << m[2] << " " << 
  	m[3] << " " << m[4] << " " << m[5];
  return os;  
}  

template <> inline
E_Float AnisoMetricType<2>::det() const
{
  return _mij[0]*_mij[2] - _mij[1]*_mij[1];
}

template <> inline
AnisoMetricType<2> AnisoMetricType<2>::inverse() const
{
  AnisoMetricType<2> INV;

  INV[0] = _mij[2];
  INV[1] = _mij[0];  
  INV[2] = -_mij[1];
	
  return INV/det();  
}

template <> inline
E_Float AnisoMetricType<3>::det() const
{
  const E_Float& a = _mij[0];
  const E_Float& b = _mij[1];
  const E_Float& c = _mij[2];
  const E_Float& d = _mij[3];
  const E_Float& e = _mij[4];
  const E_Float& f = _mij[5];

  return a*d*f - (a*e*e + d*c*c + f*b*b) + 2.*b*c*e;
}

template <> inline
AnisoMetricType<3> AnisoMetricType<3>::inverse() const
{
  AnisoMetricType<3> INV;

  const E_Float& a11 = _mij[0];
  const E_Float& a12 = _mij[1];
  const E_Float& a13 = _mij[2];
  const E_Float& a22 = _mij[3];
  const E_Float& a23 = _mij[4];
  const E_Float& a33 = _mij[5];
	
  INV[0] = a33*a22 - a23*a23;
  INV[1] = a23*a13 - a33*a12;  
  INV[2] = a23*a12 - a22*a13;  
  INV[3] = a33*a11 - a13*a13;  
  INV[4] = a12*a13 - a23*a11;  
  INV[5] = a22*a11 - a12*a12;
	
  return INV/det();  
}


}

#endif

/*
 
 
 
              NUGA 
 
 
Author : Sam Landier (sam.landier@onera.fr) 
 */

#ifndef NUGA_VECTOR_HXX
#define NUGA_VECTOR_HXX

#include "Nuga/nuga/include_c/basic_types.h"
#include "Nuga/nuga/include_c/macros.h"
#include <cmath>

#define SQRT std::sqrt

namespace NUGA
{
  template <int_t DIM>
  class mathvec_t
  {
    public:
   /// 
   /// \param x
   /// \param y
   /// \param z
   /// \param stride
   /// \return 
    template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
    static inline InputIterator3
    diff (InputIterator1 x, InputIterator2 y, InputIterator3 z, int_t stride = 1);
    
    /// 
    /// \param x
    /// \param y
    /// \param z
    static inline
    void crossProduct(const real_t* x, const real_t* y, real_t* z);
    
    /// 
    /// \param x
    /// \return 
    template <typename InputIterator>
    static inline real_t 
    normalize (InputIterator x);
    
    /// 
    /// \param x
    /// \param y
    /// \return 
    static inline real_t
    dot(const real_t* x, const real_t* y);
    
    /// 
    /// \param a
    /// \param x
    /// \param y
    /// \param z
    /// \return 
    template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
    static inline InputIterator3
    sum (InputIterator1 x, InputIterator2 y, InputIterator3 z);
    
    /// 
    /// \param a
    /// \param x
    /// \param y
    /// \param z
    /// \return 
    template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
    static inline InputIterator3
    sum (real_t a, InputIterator1 x, InputIterator2 y, InputIterator3 z);
    
    /// 
    /// \param a
    /// \param x
    /// \param b
    /// \param y
    /// \param z
    /// \return 
    template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
    static inline InputIterator3
    sum (real_t a, InputIterator1 x, real_t b, InputIterator2 y,  InputIterator3 z);
    
    /// 
    /// \param x
    /// \param y
    /// \return 
    template <typename InputIterator1, typename InputIterator2>
    static inline real_t
    L22(InputIterator1 x, InputIterator2 y);
    
    /// 
    /// \param x
    /// \return 
    template <typename InputIterator>
    static inline real_t
    L22(InputIterator x);

  };

/// 
template<>
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
mathvec_t<3>::diff
(InputIterator1 x, InputIterator2 y, InputIterator3 z, int_t stride)
{
  *z = *x - *y;
  *(z+1) = *(x+stride) - *(y+stride);
  stride *= 2;
  *(z+2) = *(x+stride) - *(y+stride);
  
  return z +3;
}

/// 
/// \param x
/// \param y
/// \param z
template <>
inline
void mathvec_t<2>::crossProduct(const real_t* x, const real_t* y, real_t* z)
{
  *z = *x * (*(y+1)) - *(x+1) * (*y);
}
/// 
/// \param x
/// \param y
/// \param z
template <>
inline
void mathvec_t<3>::crossProduct(const real_t* x, const real_t* y, real_t* z) 
{
   z[0] = x[1]*y[2] - x[2]*y[1];
   z[1] = x[2]*y[0] - x[0]*y[2];
   z[2] = x[0]*y[1] - x[1]*y[0];
}

/// 
/// \param it
/// \return
template <>
template <typename InputIterator>
inline real_t 
mathvec_t<2>::normalize (InputIterator x)
{
  real_t L = (*x * *x) + (*(x+1) * *(x+1));
  
  if (L < ZERO_MAC*ZERO_MAC) return 0.;
  
  L = 1. / SQRT(L);
  
  *x     *= L;
  *(x+1) *= L;
  
  return L;
}

/// 
/// \param it
/// \return
template <>
template <typename InputIterator>
inline real_t 
mathvec_t<3>::normalize (InputIterator x)
{
  real_t L = (*x * *x) + (*(x+1) * *(x+1)) + (*(x+2) * *(x+2));
  
  if (L < ZERO_MAC*ZERO_MAC) return 0.;
  
  L = SQRT(L);
  real_t invL = 1. / L;
  
  *x     *= invL;
  *(x+1) *= invL;
  *(x+2) *= invL;
  
  return L;
}

/// 
/// \param x
/// \param y
/// \return 
template <>
inline real_t
mathvec_t<2>::dot (const real_t* x, const real_t* y) 
{
  return (*x * (*y)) + (*(x+1) * (*(y+1)));
}

/// 
/// \param x
/// \param y
/// \return 
template <>
inline real_t
mathvec_t<3>::dot (const real_t* x, const real_t* y) 
{
  return (*x * (*y)) + (*(x+1) * (*(y+1))) + (*(x+2) * (*(y+2)));
}

// z = x + y
///
/// \param x
/// \param y
/// \param z
/// \return
template <>
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
mathvec_t<2>::sum (InputIterator1 x, InputIterator2 y, InputIterator3 z) 
{
  *z     = *(x)   + *y;
  *(z+1) = *(x+1) + *(y+1);

  return z + 2;
}

// z = ax + y
/// 
/// \param x
/// \param y
/// \param z
/// \return
template <>
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
mathvec_t<3>::sum (InputIterator1 x, InputIterator2 y, InputIterator3 z) 
{
  *z     = *(x)   + *y;
  *(z+1) = *(x+1) + *(y+1);
  *(z+2) = *(x+2) + *(y+2);

  return z + 3;
}

// z = ax + y
/// 
/// \param a
/// \param x
/// \param y
/// \param z
/// \return
template <>
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
mathvec_t<2>::sum (real_t a, InputIterator1 x, InputIterator2 y, InputIterator3 z) 
{
  *z     = *(x)   * a + *y;
  *(z+1) = *(x+1) * a + *(y+1);

  return z + 2;
}

// z = ax + y
/// 
/// \param a
/// \param x
/// \param y
/// \param z
/// \return
template <>
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
mathvec_t<3>::sum (real_t a, InputIterator1 x, InputIterator2 y, InputIterator3 z) 
{
  *z     = *(x)   * a + *y;
  *(z+1) = *(x+1) * a + *(y+1);
  *(z+2) = *(x+2) * a + *(y+2);

  return z + 3;
}

// z = ax + by
/// 
/// \param a
/// \param x
/// \param b
/// \param y
/// \param z
/// \return 
template <>
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
mathvec_t<2>::sum (real_t a, InputIterator1 x, real_t b, InputIterator2 y,  InputIterator3 z) 
{
  *z     = (*x     * a) + (*y     * b);
  *(z+1) = (*(x+1) * a) + (*(y+1) * b);

  return z + 2;
}

// z = ax + by
/// 
/// \param a
/// \param x
/// \param b
/// \param y
/// \param z
/// \return 
template <>
template <typename InputIterator1, typename InputIterator2, typename InputIterator3>
inline InputIterator3
mathvec_t<3>::sum (real_t a, InputIterator1 x, real_t b, InputIterator2 y,  InputIterator3 z) 
{
  *z     = (*x     * a) + (*y     * b);
  *(z+1) = (*(x+1) * a) + (*(y+1) * b);
  *(z+2) = (*(x+2) * a) + (*(y+2) * b);

  return z + 3;
}

/// 
/// \param x
/// \param y
/// \return 
template <>
template <typename InputIterator1, typename InputIterator2>
inline real_t
mathvec_t<2>::L22(InputIterator1 x, InputIterator2 y)
{
  return (*x - *y)*(*x - *y) + (*(x+1) - *(y+1)) * (*(x+1) - *(y+1));
}

/// 
/// \param x
/// \param y
/// \return 
template <>
template <typename InputIterator1, typename InputIterator2>
inline real_t
mathvec_t<3>::L22(InputIterator1 x, InputIterator2 y)
{
  return (*x - *y)*(*x - *y) + (*(x+1) - *(y+1)) * (*(x+1) - *(y+1)) + (*(x+2) - *(y+2)) * (*(x+2) - *(y+2));
}

/// 
/// \param x
/// \param y
/// \return 
template <>
template <typename InputIterator>
inline real_t
mathvec_t<2>::L22(InputIterator x)
{
  return (*x)*(*x) + (*(x+1)) * (*(x+1));
}

/// 
/// \param x
/// \param y
/// \return 
template <>
template <typename InputIterator>
inline real_t
mathvec_t<3>::L22(InputIterator x)
{
  return (*x)*(*x) + (*(x+1)) * (*(x+1)) + (*(x+2)) * (*(x+2));
}


} //NUGA

#endif
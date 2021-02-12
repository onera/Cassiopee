/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __DELAUNAY_INTERPOLATOR_H__
#define __DELAUNAY_INTERPOLATOR_H__

#include "Nuga/include/defs.h"

namespace DELAUNAY
{

template <typename T>
class Interpolator
{
public:

  Interpolator(void)
  {
  }

  virtual ~Interpolator(void)
  {
  }

  virtual T interpolate(const T& H0, const T& H1, E_Float u) const = 0;
};

template <typename T>
class LinearInterpolator : public Interpolator<T>
{
public:
  LinearInterpolator(void)
  {
  }

  ~LinearInterpolator(void)
  {
  }

  inline T interpolate(const T& H0, const T& H1, E_Float u) const {return H0*(1-u) + H1*u;}

};

template <typename T>
class GeometricInterpolator : public Interpolator<T>
{
public:
  GeometricInterpolator(void)
  {
  }

  ~GeometricInterpolator(void)
  {
  }

  inline T interpolate(const T& H0, const T& H1, E_Float u) const {return H0 * ::pow(H1/H0, u);} //iso

};

//aniso
template<> inline
DELAUNAY::Aniso2D
GeometricInterpolator<DELAUNAY::Aniso2D>::interpolate
(const DELAUNAY::Aniso2D& H0, const DELAUNAY::Aniso2D& H1, E_Float u) const
{
  //fixme : in fact Linear currently
  return H0*(1-u) + H1*u;
}

} // namespace DELAUNAY

#endif

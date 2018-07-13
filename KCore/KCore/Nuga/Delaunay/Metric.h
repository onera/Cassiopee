/*    
    Copyright 2013-2018 Onera.

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

#ifndef _DELAUNAY_METRIC_H_
#define _DELAUNAY_METRIC_H_

#include "Def/DefContainers.h"
#include "Def/DefFunction.h"
#include "AnisoMetricType.h"
#include "Interpolator.h"
#include "Linear/DelaunayMath.h"
#include "macros.h"
#include "MeshUtils1D.h"

//#include<list>

namespace DELAUNAY{

  /** Iso/Aniso Variable Metric Class.

  Type T is either a E_Float for an isotropic metric field or an Symmetric matrix (i.e. AnisoMetricType)
  for an anisotropic metric field.
  */
  template <typename T>
  class VarMetric
  {

  public: /** Typedefs */

    typedef K_CONT_DEF::size_type size_type;
    typedef VarMetric             self_type;
    typedef T                     value_type;
    typedef std::vector<T>        field_type;

  public:
    enum eInterpType {LINEAR, GEOMETRIC};

  public: /** Constructors and Destructor */

    VarMetric(K_FLD::FloatArray& pos, E_Float hmin, E_Float hmax, eInterpType interp_type = LINEAR);

    virtual ~VarMetric(void){if (_interpol){delete _interpol; _interpol = 0;} }

  public: /** API methods */

    inline virtual value_type computeMetric(size_type N, size_type Ni, size_type Nj, E_Float r);

    inline virtual void setMetric(E_Int N, const T& m);

    inline E_Float lengthEval (size_type Ni, const T& mi, size_type Nj, const T& mj);

    inline E_Float length (size_type Ni, size_type Nj, const E_Float& threshold, K_CONT_DEF::int_vector_type& tmpNodes);

    inline E_Float getRadius(size_type Ni);

    virtual void init_metric
      (const K_FLD::FloatArray& metric, K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
       const std::vector<E_Int>& hard_nodes);

    inline const value_type& operator[](size_type i) const { assert (i < (size_type)_field.size()); return _field[i];}

    ///inline void smooth(size_type Ni, size_type Nj, E_Float eps);
    
  protected:

    inline void compute_intersection(field_type& metric1, const K_FLD::FloatArray& metric2);

    inline void convertIsoToAniso(const std::vector<E_Float>& isoM, std::vector<T>& anisoM);

    inline E_Bool isValidMetric(const T& mi);

    inline void setUserMetric(const K_FLD::FloatArray& Umetric, field_type& metric);

  protected:

    ///
    E_Float                   _hmin;
    ///
    E_Float                   _hmax;
    /// Metric vector defined at each node in pos.
    field_type                _field;
    /// Coordinate array.
    K_FLD::FloatArray&        _pos;
    ///
    Interpolator<T>*          _interpol;

  };

  // Implementations.

  //Constructors
  template <typename T>
  VarMetric<T>::VarMetric(K_FLD::FloatArray& pos, E_Float hmin, E_Float hmax, eInterpType interp_type)
    :_hmin(hmin), _hmax(hmax), _pos(pos), _interpol(new LinearInterpolator<T>)
  {
  }

   template <typename T>
   void
   VarMetric<T>:: init_metric
   (const K_FLD::FloatArray& metric, K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
    const std::vector<E_Int> & hard_nodes)
   {
     // Compute the isometric metric based on the boundary.
     std::vector<E_Float> m_iso;
     E_Float hmin1, hmax1;
     MeshUtils1D::compute_iso_metric(pos, connectB, hard_nodes, m_iso, hmin1, hmax1);
     
     // Set _hmin and _hmax if not done nor correct.
     if ((_hmin <= 0.) || (_hmin == K_CONST::E_MAX_FLOAT)) // If not (or wrongly) set by the user.
       _hmin = hmin1;
     if ((_hmax <= 0.) || (_hmax == K_CONST::E_MAX_FLOAT)) // If not (or wrongly) set by the user.
       _hmax = hmax1;
     
     if (_hmax != hmax1 || _hmin != hmin1)
     {
       for (size_t i=0; i < m_iso.size(); ++i)
       {
         m_iso[i] = std::min(m_iso[i], _hmax);
         m_iso[i] = std::max(m_iso[i], _hmin);
       }
     }

     // Set _metric by converting m_iso to an aniso type metric.
     // Only relevant in aniso case (in iso just assign it).
     convertIsoToAniso(m_iso, _field);

     // Now take the user metric into account at each node when it's valid.
     // Only the first row of the input matrix is taken into account.
     //compute_intersection(_metric, metric);
     setUserMetric (metric, _field);

   }

  template <typename T>
  E_Float
    VarMetric<T>::length(size_type Ni, size_type Nj, const E_Float & threshold, K_CONT_DEF::int_vector_type& tmpNodes){

      // This computes recursively the length of the input edge by splitting it in
      // a number of bits such all the bits are smaller than the treshold.

      E_Float L = lengthEval(Ni, _field[Ni], Nj, _field[Nj]);
      if (L <= threshold)
        return L;

      // Split
      size_type dim = _pos.rows();
      E_Float* newP = new E_Float[dim];
      E_Float* it   = newP;

      K_FLD::FloatArray::const_iterator it1(_pos.col(Ni)), it2(_pos.col(Nj));

      for (size_type i = 0; i < dim; ++i) // middle point.
        *(it++) = 0.5 * (*(it1++) + *(it2++));

      _pos.pushBack(newP, newP+dim);
      size_type N = _pos.cols()-1;
      delete [] newP;

      tmpNodes.push_back(N);

      this->setMetric(N, this->computeMetric(N, Ni, Nj, 0.5));//fixme

      return (self_type::length(Ni, N, threshold, tmpNodes) + self_type::length(N, Nj, threshold, tmpNodes));
  }

  // Methods specialisation

  static E_Float v[2];
  static E_Float vi[2];
  static E_Float vj[2];
  static E_Float r1;
  static E_Float r2;

  /** Anisotropic Variable field. */
  ///
  template <typename T> inline
    E_Float
    VarMetric<T>::lengthEval (size_type Ni, const T& mi, size_type Nj, const T& mj)
  {
    r1 = r2 = 0.;
    
    K_FUNC::diff<2> (_pos.col(Nj), _pos.col(Ni), v);

    vi[0] = mi[0]*v[0] + mi[1]*v[1];
    vi[1] = mi[1]*v[0] + mi[2]*v[1];
    vj[0] = mj[0]*v[0] + mj[1]*v[1];
    vj[1] = mj[1]*v[0] + mj[2]*v[1];
    /*
    vi = mi * v;
    vj = mj * v;
    */
    for (K_CONT_DEF::size_type i = 0; i < _pos.rows(); ++i)
    {
      r1 += vi[i]*v[i];
      r2 += vj[i]*v[i];
    }

    return 0.5 * (::sqrt(r1) + ::sqrt(r2)); //integral approx.
  }

  /** Isotropic variable field. */
  /*
  ///
  template <>
  E_Float
  VarMetric<E_Float>::length(size_type Ni, size_type Nj)
  {
  E_Float h0 = _metric[Ni], h1 = _metric[Nj];
  E_Float d = ::sqrt(K_FUNC::sqrDistance(_pos.col(Ni), _pos.col(Nj), _pos.rows()));
  E_Float r1 = h1 / h0;
  if (::abs(r1 - 1.) < 0.01)
  return d / std::max(h0,h1);
  E_Float   r = (d + h1) / (d + h0);
  size_type n = size_type(::log(r1)/::log(r));
  r1 = ::pow (r1, 1./E_Float(n));
  return (h1 - r1 * h0)/(r1 - 1);
  }
  */
  
  ///
  template <> inline
    E_Float
    VarMetric<E_Float>::lengthEval (size_type Ni, const E_Float& mi, size_type Nj, const E_Float& mj)
  {
    // Warning : input mi and mj are in fact hi and hj.
    return 0.5 * ((1./mi) + (1./mj)) * ::sqrt(K_FUNC::sqrDistance(_pos.col(Ni), _pos.col(Nj), _pos.rows()));
  }
  
  ///
  template<> inline
    E_Float
    VarMetric<E_Float>::getRadius(size_type Ni)
  {
    return _field[Ni];
  }

  ///
  template<> inline
    E_Float
    VarMetric<Aniso2D>::getRadius(size_type Ni)
  {
    E_Float l1, l2;
    K_LINEAR::DelaunayMath::eigen_values (_field[Ni][0], _field[Ni][2], _field[Ni][1], l1,l2);
    l1 = (l1 < l2) ? l1 : l2;
    return 1./::sqrt(l1);
  }

  ///
  template <typename T> inline
    T
    VarMetric<T>::computeMetric(size_type N, size_type Ni, size_type Nj, E_Float r)
  {
    return _interpol->interpolate(_field[Ni], _field[Nj], r); //fixme
    // need to be in the range [hmin;hmax] if defined...
  }

  ///
  template <typename T> inline
    void
    VarMetric<T>::setMetric(E_Int N, const T& m)
  {
    if (N > (E_Int)_field.size())
      _field.resize(N);
    _field.push_back(m);
  }

/*
  ///
  template<> inline
    void VarMetric<E_Float>::smooth(size_type Ni, size_type Nj, E_Float eps)
  {
    E_Float mi = _field[Ni];
    E_Float mj = _field[Nj];
    E_Float d = ::sqrt(K_FUNC::sqrDistance(_pos.col(Ni), _pos.col(Nj), _pos.rows()));
    _field[Ni] = std::min (mi, mj + eps*d);
    _field[Nj] = std::min (mj, mi + eps*d);
  }

  ///
  template<> inline
    void VarMetric<Aniso2D>::smooth(size_type Ni, size_type Nj, E_Float eps)
  {

    K_FLD::FloatArray M1(2,2),M2(2,2), M1e(2,2), M2e(2,2);
    E_Float d = ::sqrt(K_FUNC::sqrDistance(_pos.col(Ni), _pos.col(Nj), _pos.rows()));

    M1(0,0) = _field[Ni][0];
    M1(1,1) = _field[Ni][2];
    M1(1,0) = M1(0,1) = _field[Ni][1];

    //  std::cout << M1 << std::endl;

    M1e(0,0) = M1(0,0) + eps * d;
    M1e(1,1) = M1(1,1) + eps * d;
    M1e(1,0) = M1e(0,1) = M1(1,0);

    M2(0,0) = _field[Nj][0];
    M2(1,1) = _field[Nj][2];
    M2(1,0) = M2(0,1) = _field[Nj][1];

    M2e(0,0) = M2(0,0) + eps * d;
    M2e(1,1) = M2(1,1) + eps * d;
    M2e(1,0) = M2e(0,1) = M2(1,0);

    //  std::cout << M2 << std::endl;

    K_FLD::FloatArray I;
    K_LINEAR::DelaunayMath::intersect(M1, M2e, I);

    //  std::cout << I << std::endl;

    _field[Ni][0] = I(0,0);
    _field[Ni][1] = I(1,0);
    _field[Ni][2] = I(1,1);

    K_LINEAR::DelaunayMath::intersect(M2, M1e, I);

    //  std::cout << I << std::endl;

    _field[Nj][0] = I(0,0);
    _field[Nj][1] = I(1,0);
    _field[Nj][2] = I(1,1); 
  }
*/

  template<> inline
    void VarMetric<E_Float>::compute_intersection(std::vector<E_Float>& metric1,
    const K_FLD::FloatArray& metric2)
  {
    E_Float m;
    E_Int max = std::min((E_Int)metric1.size(), (E_Int)metric2.cols());
    for (size_type i = 0; i < max; ++i)
    {
      m = metric2(0,i);
      if ((m > 0.) && (m < K_CONST::E_MAX_FLOAT))
        metric1[i] = std::min(m, metric1[i]) ;
    }
  }

  template<> inline
    void VarMetric<Aniso2D>::compute_intersection(std::vector<Aniso2D>& metric1,
    const K_FLD::FloatArray& metric2)
  {  
    // Now set the metric at each node as the minimum(intersection) between the edge length and the input metric.
    // Only the first row of the input matrix is taken into account.
    K_FLD::FloatArray M1(2,2), M2(2,2), I(2,2);
    Aniso2D m;

    E_Int max = std::min((E_Int)metric1.size(), (E_Int)metric2.cols());

    // fixme need to check validity of m2
    for (E_Int i = 0; i < max; ++i)
    {
      M1(0,0) = metric1[i][0];
      M1(1,0) = M1(0,1) = metric1[i][1];
      M1(1,1) = metric1[i][2];

      if (metric2.rows() >= 3)
      {
        M2(0,0) = metric2(0,i);
        M2(1,0) = M2(0,1) = metric2(1,i);
        M2(1,1) = metric2(2,i);
      }
      else // Use only the first row as an iso metric
      {
        M2(0,0) = metric2(0,i);
        M2(1,0) = M2(0,1) = 0.;
        M2(1,1) = metric2(0,i);
      }

      K_LINEAR::DelaunayMath::intersect(M1, M2, I);

      m[0] = I(0,0);
      m[1] = I(1,0);
      m[2] = I(1,1);
      metric1[i] = m;
    }
  }

  template<> inline
  void
  VarMetric<Aniso2D>::convertIsoToAniso
  (const std::vector<E_Float>& isoM, std::vector<Aniso2D>& anisoM)
  {
    anisoM.resize(isoM.size());
    E_Float im;
    for (size_t i = 0; i < isoM.size(); ++i)
    {
      im = isoM[i];
      anisoM[i][0] = anisoM[i][2] = (im > 0.) ? 1./(im*im) : 0.;
      anisoM[i][1] = 0.;
    }
  }

  template<> inline
  void
  VarMetric<E_Float>::convertIsoToAniso
  (const std::vector<E_Float>& isoM, std::vector<E_Float>& anisoM)
  {
    anisoM = isoM;
  }

  template<> inline
  E_Bool
  VarMetric<E_Float>::isValidMetric(const E_Float& mi)
  {
    return ((mi > 0.) && (mi < K_CONST::E_MAX_FLOAT));
  }

  template<> inline
  E_Bool
  VarMetric<Aniso2D>::isValidMetric(const Aniso2D& mi)
  {
    const E_Float& a11 = mi[0];
    const E_Float& a12 = mi[1];
    const E_Float& a22 = mi[2];
    E_Float det = (a11*a22) - (a12*a12);

    return ((a11 > 0.) && (a22 > 0.) && (det > 0.));
  }

  template<> inline
  void
  VarMetric<E_Float>::setUserMetric(const K_FLD::FloatArray& Umetric, field_type& metric)
  {
    E_Int max = std::min(Umetric.cols(), (E_Int)metric.size()); 
    for (E_Int i = 0; i < max; ++i)
    {
      if (isValidMetric(Umetric(0,i)))
        metric[i] = Umetric(0,i);
    }
  }

  template<> inline
  void
  VarMetric<Aniso2D>::setUserMetric(const K_FLD::FloatArray& Umetric, field_type& metric)
  {
    E_Int max = std::min(Umetric.cols(), (E_Int)metric.size()); 
    Aniso2D m;
    for (E_Int i = 0; i < max; ++i)
    {
      m = Umetric.col(i);

      if (isValidMetric(m))
        metric[i] = m;
    }
  }

} // End namespace DELAUNAY

#endif /* DELAUNAY_METRIC_H_ */

/*    
    Copyright 2013-2019 Onera.

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

#ifdef DEBUG_METRIC
#include "Linear/DelaunayMath.h"
#include "iodata.h"
#include "IO/io.h"
#include <sstream>
#endif

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

    inline virtual void computeMetric(size_type N, size_type Ni, size_type Nj, E_Float r);

    inline virtual void setMetric(E_Int N, const T& m);

    inline E_Float lengthEval (size_type Ni, const T& mi, size_type Nj, const T& mj);

    inline E_Float length (size_type Ni, size_type Nj, const E_Float& threshold, K_CONT_DEF::int_vector_type& tmpNodes);

    inline E_Float getRadius(size_type Ni);
    
    /// computes the square of the distance between the center and a point P on the ellipse along direction dir.
    inline E_Float get_h2_along_dir(size_type Ni, const E_Float* dir);
    ///
    inline void directional_metric_reduce(size_type Ni, E_Float h2, const E_Float* normed_dir);
    ///
    inline void metric_reduce(size_type Ni, const E_Float* NiNj, E_Float hjold2, E_Float hjnew2);

    virtual void init_metric
      (const K_FLD::FloatArray& metric, K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
       const std::vector<E_Int>& hard_nodes);

    inline const value_type& operator[](size_type i) const { assert (i < (size_type)_field.size()); return _field[i];}
    
    void __init_refine_points
    (K_FLD::FloatArray& pos, size_type Ni, size_type Nj, E_Float threshold,
     std::vector<std::pair<E_Float, size_type> >& length_to_points, std::vector<size_type>& tmpNodes);
    
    void __compute_refine_points
    (K_FLD::FloatArray& pos, size_type Ni, size_type Nj, E_Float threshold, std::vector<std::pair<E_Float, size_type> >& length_to_points, std::vector<size_type>& tmpNodes);
    
    void smoothing_loop(const std::set<K_MESH::NO_Edge>& edges, E_Float eps, E_Int itermax, E_Int N0 /* threshold for metric changes*/);
    void smoothing_loop(const K_FLD::IntArray& connectB, E_Float eps, E_Int itermax, E_Int N0 /* threshold for metric changes*/);

    inline bool smooth(size_type Ni, size_type Nj, E_Float gr, E_Int N0 /* threshold for metric changes*/);
    
    bool is_valid() {
      bool res= true; E_Int i=0; 
      for (; (i < _field.size()) && res; ++i) res &= isValidMetric(_field[i]); 
      if (i != _field.size()) std::cout << "failed at " << i-1 << std::endl; 
      return res;
    }


#ifdef DEBUG_METRIC
  void append_unity_ellipse(const K_FLD::FloatArray& c, E_Int i, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int Nc = -1);
  void draw_ellipse_field(const char* fname, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const std::vector<bool>* mask = 0);
  void draw_ellipse_field(const char* fname, const K_FLD::FloatArray& crd, const std::vector<size_type>& indices);
  void draw_ellipse(const char* fname, const K_FLD::FloatArray& crd, size_type i);
#endif
    
  protected:

    inline void compute_intersection(field_type& metric1, const K_FLD::FloatArray& metric2);

    inline void convertIsoToAniso(const std::vector<E_Float>& isoM, std::vector<T>& anisoM);

    inline bool isValidMetric(const T& mi);

    inline bool isValidMetric(const E_Float* mi);
    
    inline bool isIsotropic(const T& mi) { return true;}
    
    inline bool isWeakAniso(size_type Ni, E_Float r) { return true;}

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
  public:
    E_Int _N0;
    
    

  };
  
  /// 
  template <> inline
  bool
  VarMetric<Aniso2D>::isValidMetric(const E_Float* mi)
  {
    const E_Float& a11 = mi[0];
    const E_Float& a12 = mi[1];
    const E_Float& a22 = mi[2];
    E_Float det = (a11*a22) - (a12*a12);    
    
    return ((a11 > 0.) && (a22 > 0.) && (det > 0.));
  }
  
    template<> inline
  bool
  VarMetric<E_Float>::isValidMetric(const E_Float& mi)
  {
    return ((mi > 0.) && (mi < K_CONST::E_MAX_FLOAT));
  }

  template<> inline
  bool
  VarMetric<Aniso2D>::isValidMetric(const Aniso2D& mi)
  {
    const E_Float& a11 = mi[0];
    const E_Float& a12 = mi[1];
    const E_Float& a22 = mi[2];
    E_Float det = (a11*a22) - (a12*a12);

    return ((a11 > 0.) && (a22 > 0.) && (det > 0.));
  }
  
  template<> inline
  bool
  VarMetric<Aniso2D>::isIsotropic(const Aniso2D& mi) {
    const E_Float& a11 = mi[0];
    const E_Float& a12 = mi[1];
    const E_Float& a22 = mi[2];
    
    if ( ::fabs(a12) > E_EPSILON) return false;
    if ( ::fabs(a11 - a22) > E_EPSILON) return false;
    
    return true;
  }
  
  /// computes the intersection between an ellipse and a line
  template<> inline
    E_Float
    VarMetric<Aniso2D>::get_h2_along_dir(size_type Ni, const E_Float* dir)
  {
    E_Float L2 = K_FUNC::sqrNorm<2>(dir);
    
    const E_Float& m11 = _field[Ni][0];
    const E_Float& m12 = _field[Ni][1];
    const E_Float& m22 = _field[Ni][2];
    
    E_Float k2 = 1. /(m11*dir[0]*dir[0] + 2.* m12*dir[0]*dir[1] + m22*dir[1]*dir[1]);

    return k2*L2;
  }
  
  template<> inline
  bool
  VarMetric<Aniso2D>::isWeakAniso(size_type Ni, E_Float r) { 
    
    const Aniso2D& mi = _field[Ni];
     
    E_Float lambda0, lambda1, v0[2], v1[2];
    K_LINEAR::DelaunayMath::eigen_vectors (mi[0], mi[2], mi[1], lambda0, lambda1, v0, v1);
    E_Float h1old2 = get_h2_along_dir(Ni, v0);
    E_Float h2old2 = get_h2_along_dir(Ni, v1);
    
    E_Float rr = (std::min(h1old2, h2old2)/std::max(h1old2, h2old2));
    
    return (  rr < r*r );
  }
    
#ifdef DEBUG_METRIC
  
  template <typename T> inline
  void VarMetric<T>::draw_ellipse_field
  (const char* fname, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const std::vector<bool>* mask)
  {
    K_FLD::IntArray cnto;
    K_FLD::FloatArray crdo;
    
    std::vector<E_Int> indices;
    
    if (mask)
    {
      std::set<E_Int> unodes;
      for (size_t i=0; i < cnt.cols(); ++i)
      {
        if (!(*mask)[i]) continue;
        unodes.insert(cnt(0,i));
        unodes.insert(cnt(1,i));
        unodes.insert(cnt(2,i));
      }
      
      indices.insert(indices.end(), unodes.begin(), unodes.end());
    }
    else
      cnt.uniqueVals(indices);
    
    for (size_t i = 0; i < indices.size(); ++i)
      append_unity_ellipse(crd, indices[i], crdo, cnto);
    
    MIO::write(fname, crdo, cnto, "BAR");
  }
  
  template <typename T> inline
  void VarMetric<T>::draw_ellipse
  (const char* fname, const K_FLD::FloatArray& crd, size_type i)
  {
    K_FLD::IntArray cnto;
    K_FLD::FloatArray crdo;
    
    append_unity_ellipse(crd, i, crdo, cnto);
    
    MIO::write(fname, crdo, cnto, "BAR");
  }
  
  template <typename T> inline
  void VarMetric<T>::draw_ellipse_field
  (const char* fname, const K_FLD::FloatArray& crd, const std::vector<size_type>& indices)
  {
    K_FLD::IntArray cnto;
    K_FLD::FloatArray crdo;
    
    for (size_t i = 0; i < indices.size(); ++i)
      append_unity_ellipse(crd, indices[i], crdo, cnto);
    
    MIO::write(fname, crdo, cnto, "BAR");
  }
  
  template<> inline
  void
  VarMetric<Aniso2D>::append_unity_ellipse(const K_FLD::FloatArray& c, E_Int i, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int Nc)
  {
    const Aniso2D& m = this->_field[i];
    
    if (m[0] == 0. || m[2] == 0. ) return;
    
    E_Int SAMPLE = 100;
    
    E_Float alpha = 2. * K_CONST::E_PI /  (E_Float)SAMPLE;
    
    E_Int pos0 = crd.cols();
    
    if (Nc == -1) Nc = i;
    
    for (size_t n=0; n < SAMPLE; ++n)
    {
      E_Float a = alpha * n;
      
      E_Float V[] = {::cos(a), ::sin(a)};

      const E_Float& u=V[0];
      const E_Float& v=V[1];
      const E_Float& m11 = m[0];
      const E_Float& m12 = m[1];
      const E_Float& m22 = m[2];
      
      E_Float k = 1. /::sqrt( u*(m11*u + m12*v) + v*(m12*u + m22*v) );
      
      E_Float Pt[] = {k*u , k*v};
      
      K_FUNC::sum<2>(Pt, c.col(Nc), Pt); //center it at node Nc
      
      crd.pushBack(Pt, Pt+2);
      
    }
    
    for (size_t n=0; n < SAMPLE; ++n)
    {
      E_Int e[] = {pos0+n, pos0+(n+1)%SAMPLE};
      cnt.pushBack(e, e+2);
    }   
  }
  
  template<> inline
  void
  VarMetric<E_Float>::append_unity_ellipse(const K_FLD::FloatArray& c, E_Int i, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, E_Int Nc)
  {    
    E_Int SAMPLE = 50;
    
    const E_Float& h = this->_field[i];
    
    E_Float alpha = 2. * K_CONST::E_PI /  (E_Float)SAMPLE;
    
    E_Int pos0 = crd.cols();
    
    if (Nc == -1) Nc = i;
    
    for (size_t n=0; n < SAMPLE; ++n)
    {
      E_Float a = alpha * n;
      E_Float V[] = {::cos(a), ::sin(a)};
      
      E_Float Pt[] = {h*V[0] , h*V[1]};
      
      K_FUNC::sum<2>(Pt, c.col(Nc), Pt); //center it at node Nc
      
      crd.pushBack(Pt, Pt+2);
      
    }
    
    for (size_t n=0; n < SAMPLE; ++n)
    {
      E_Int e[] = {pos0+n, pos0+(n+1)%SAMPLE};
      cnt.pushBack(e, e+2);
    }   
  }
  
#endif

  // Implementations.

  //Constructors
  template <typename T>
  VarMetric<T>::VarMetric(K_FLD::FloatArray& pos, E_Float hmin, E_Float hmax, eInterpType interp_type)
    :_hmin(hmin), _hmax(hmax), _pos(pos)
  {
    if (interp_type == LINEAR)
      _interpol = new LinearInterpolator<T>;
    else //GEOMETRIC
      _interpol = new GeometricInterpolator<T>;
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
     if (_hmax <= 0.)
       _hmax = K_CONST::E_MAX_FLOAT;
     if ((_hmin <= 0.) || (_hmin > _hmax) || (_hmin == K_CONST::E_MAX_FLOAT) )
       _hmin = hmin1;
     
     if (_hmax != K_CONST::E_MAX_FLOAT)
     {
       for (size_t i=0; i < m_iso.size(); ++i)
         m_iso[i] = std::min(m_iso[i], _hmax);
     }
     if (_hmin != hmin1)
     {
       for (size_t i=0; i < m_iso.size(); ++i)
         m_iso[i] = std::max(m_iso[i], _hmin);
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

      this->computeMetric(N, Ni, Nj, 0.5);

      return (self_type::length(Ni, N, threshold, tmpNodes) + self_type::length(N, Nj, threshold, tmpNodes));
  }

  // Methods specialisation

 

  /** Anisotropic Variable field. */
  ///
  template <typename T> inline
    E_Float
    VarMetric<T>::lengthEval (size_type Ni, const T& mi, size_type Nj, const T& mj)
  {
    E_Float v[2];
    E_Float vi[2];
    E_Float vj[2];
    E_Float r1;
    E_Float r2;
    r1 = r2 = 0.;
    
#ifdef DEBUG_METRIC
    assert (isValidMetric(mi));
    assert (isValidMetric(mj));
#endif
    
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
    
    E_Float res = 0.5 * (::sqrt(r1) + ::sqrt(r2)); //integral approx.

#ifdef DEBUG_METRIC
    assert (res > E_EPSILON);
#endif

    return res; 
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
    E_Float lmax,lmin;
    _field[Ni].eigen_values(lmax, lmin);
    return 1./::sqrt(lmin);
  }
  
  ///
  template<> inline
  void
  VarMetric<Aniso2D>::directional_metric_reduce(size_type Ni, E_Float hdir2, const E_Float* normed_dir)
  {
    // Aim : modifying the ellipse at Ni such d(Ni,X)*d(Ni,X) = h2 
    //       where X is the intersection point of the ellipse with the line (Ni, dir)
    
    // (h2, dir) => X => hnew, the one between h1new and h2new that maximize hnew/hold
    
    E_Float hdir = ::sqrt(hdir2);
        
    E_Float NiX[] = {normed_dir[0]*hdir, normed_dir[1]*hdir};
    
    
    // Diagonalization
    E_Float lambda0, lambda1, v0[2], v1[2];
    K_LINEAR::DelaunayMath::eigen_vectors (_field[Ni][0], _field[Ni][2], _field[Ni][1], lambda0, lambda1, v0, v1);
    K_FLD::FloatArray D(2,2, 0.);
    D(0,0) = lambda0;
    D(1,1) = lambda1;
    
    K_FUNC::normalize<2>(v0);
    K_FUNC::normalize<2>(v1);
    
    // transformation Matrix : diagonalizing base (BD) : Main axis -> (i,j)
    K_FLD::FloatArray P(2,2, 0.), tP(2,2,0.);
    tP(0,0) = P(0,0) = v0[0];
    tP(0,1) = P(1,0) = v0[1];
    tP(1,0) = P(0,1) = v1[0];
    tP(1,1) = P(1,1) = v1[1];
    
    // X expressed in BD
    E_Float Xbd[2];
    Xbd[0] = tP(0,0) * NiX[0] + tP(0,1) * NiX[1];
    Xbd[1] = tP(1,0) * NiX[0] + tP(1,1) * NiX[1];
    
    // Computing h1old and h2old to get first fundamental form coefficient (to keep them in the transfo)
    // E = lambda0 * h1old2
    // G = lambda1 * h2old2
    E_Float h1old2 = get_h2_along_dir(Ni, v0);
    E_Float h2old2 = get_h2_along_dir(Ni, v1);
    
    // Computing h1new and h2new
    E_Float xbd2 = Xbd[0]*Xbd[0];
    E_Float ybd2 = Xbd[1]*Xbd[1];
    
    if (xbd2 < E_EPSILON*E_EPSILON)       // along second axis
      D(1,1) *= (h2old2 / ybd2);
    else if (ybd2 < E_EPSILON*E_EPSILON) // along first axis
      D(0,0) *= (h1old2 / xbd2);
//    else if (::fabs(lambda0-lambda1) < E_EPSILON) //isotropic reduce
//    {
//      //assert (::fabs(xbd2 - ybd2) < E_EPSILON*E_EPSILON);
//      D(0,0) = 1. / xbd2;
//      D(1,1) = 1. / ybd2;
//    }
    else
    {
      bool first_axis = (ybd2/xbd2 < ::fabs(lambda0/lambda1));
//      E_Float psu2 = K_FUNC::dot<2>(normed_dir, v0);
//      psu2 *= psu2;
//      E_Float psv2 = K_FUNC::dot<2>(normed_dir, v1);
//      psv2 *= psv2;
//      bool first_axis = (psu2 >= psv2);
//      E_Float l0new = ::fabs((1. - lambda1*(Xbd[1]*Xbd[1])) / (Xbd[0]*Xbd[0]));
//      E_Float l1new = ::fabs((1. - lambda0*(Xbd[0]*Xbd[0])) / (Xbd[1]*Xbd[1]));

      // Replace the eigen value to reflect the transfo
      if (first_axis)
        D(0,0) *= (h1old2 * ::fabs(1. - lambda1 * ybd2) / xbd2);
      else
        D(1,1) *= (h2old2 * ::fabs(1. - lambda0 * xbd2) / ybd2);
    }

    K_FLD::FloatArray new_metric(2,2, 0.);
    new_metric = P * D * tP;
    
#ifdef DEBUG_METRIC
    E_Float sym_check = ::fabs(new_metric(0,1) - new_metric(1,0));
    if (sym_check >=E_EPSILON)
    {
      sym_check = sym_check / std::max(new_metric(0,1), new_metric(1,0));
    }
    //std::cout << sym_check << std::endl;
    assert (sym_check < E_EPSILON);
//if (Ni == 112)
//    {
//    K_FLD::FloatArray crdo;
//    K_FLD::IntArray cnto;
//    append_unity_ellipse(_pos, Ni, crdo, cnto);
//
//    std::ostringstream o;
//    o << "metric_at_Ni_before" << Ni << ".mesh";
//    MIO::write(o.str().c_str(), crdo, cnto, "BAR");
//    
//    }
#endif

    _field[Ni][0] = new_metric(0,0);
    _field[Ni][1] = new_metric(0,1);
    _field[Ni][2] = new_metric(1,1);
    
#ifdef DEBUG_METRIC
    
//    if (Ni == 112)
//    {
//    std::ostringstream o;
//    o << "metric_at_Ni_after" << Ni << ".mesh";
//    draw_ellipse(o.str().c_str(), _pos, Ni);
//    }
    assert (isValidMetric( _field[Ni]));
#endif
  }
  
  template<> inline
  void
  VarMetric<Aniso2D>::metric_reduce(size_type Nj, const E_Float* normed_dir, E_Float hjold2, E_Float hjnew2)
  {
//    if (isValidMetric(_field[Nj]) && !isWeakAniso(Nj, 0.25))
//      directional_metric_reduce(Nj, hjnew2, normed_dir); //orientation of NiNj does not matter
//    else //isotropic reduction
    {
      _field[Nj][0] *= ( hjold2 / hjnew2) ;
      _field[Nj][2] *= ( hjold2 / hjnew2);
      _field[Nj][1] *= ( hjold2 / hjnew2) ;
    }
  }

  ///
  template <typename T> inline
  void
    VarMetric<T>::computeMetric(size_type N, size_type Ni, size_type Nj, E_Float r)
  {
    setMetric(N, _interpol->interpolate(_field[Ni], _field[Nj], r));
  }

  ///
  template <typename T> inline
    void
    VarMetric<T>::setMetric(E_Int N, const T& m)
  {
    if (N >= (E_Int)_field.size())
      _field.resize(N+1);
    _field[N] = m;
  }
  
  template <typename T> inline
  void
  VarMetric<T>::smoothing_loop
  (const std::set<K_MESH::NO_Edge>& edges, E_Float gr, E_Int itermax, E_Int N0 /* threshold for metric changes*/)
  {
    E_Int iter(0);
    bool has_changed = false;
    //
    do
    {
      has_changed = false;
      for (const auto& Ei : edges)
        has_changed |= this->smooth(Ei.node(0), Ei.node(1), gr, N0);
    }
    while ( has_changed && (++iter < itermax) );
  }
  
  template <typename T> inline
  void
  VarMetric<T>::smoothing_loop
  (const K_FLD::IntArray& connectB, E_Float gr, E_Int itermax, E_Int N0 /* threshold for metric changes*/)
  {
    E_Int iter(0);
    bool has_changed = false;
    //
    do
    {
      has_changed = false;
      for (E_Int i=0; i < connectB.cols(); ++i)
        has_changed |= this->smooth(connectB(0,i), connectB(1,i), gr, N0);
    }
    while ( has_changed && (++iter < itermax) );
  }
  
  ///
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
  
  ///
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
  
  ///
  template<> inline
    bool VarMetric<E_Float>::smooth(size_type Ni, size_type Nj, E_Float gr, E_Int N0 /* threshold for metric changes*/)
  {
    //WARNING : assume growth ratio gr > 1.
    
    // The smaller might smooth the bigger
    
    E_Float hi0 = _field[Ni];
    E_Float hj0 = _field[Nj];

    if (::fabs(hj0 - hi0) < E_EPSILON) return false; //same metric so nothing to smooth
    if (hj0 < hi0)
    {
      std::swap(Ni,Nj);
      std::swap(hi0, hj0);
    }
    
    // (Ni,mi) is the finer now, Nj is the one to smooth
    
    if (Nj <= N0) return false; //do not touch the hard nodes

    E_Float dij = lengthEval(Ni, hi0, Nj, hi0);
    
    if (hi0 >= dij && hj0 >= dij) return false; //fixme : ??

    E_Float hs = hi0 *(1. + (gr-1.) * dij); // extrapolated h at Nj with respect to growth ratio

    hs = std::max(_hmin, hs);
    hs = std::min(_hmax, hs);
    
    if (hj0 <= hs) return false; //false if metric at Nj is smaller than the replacement one.
    
    _field[Nj] = hs;
            
    return true;
    
  }
  
  ///
  template<> inline
    bool VarMetric<Aniso2D>::smooth(size_type Ni, size_type Nj, E_Float gr, E_Int N0 /* threshold for metric changes*/)
  {
    //WARNING : assume growth ratio gr > 1.
    
    // The smaller might smooth the bigger. Discuss on the spectral radius
    
    const Aniso2D& mi0 = _field[Ni];
    const Aniso2D& mj0 = _field[Nj];

    E_Float NiNj[2];
    K_FUNC::diff<2>(_pos.col(Nj), _pos.col(Ni), NiNj);
       
    E_Float hi02 = get_h2_along_dir(Ni, NiNj); // trace on NiNj of the ellipse centered at Ni
    E_Float hj02 = get_h2_along_dir(Nj, NiNj); // trace on NiNj of the ellipse centered at Nj
    
    const Aniso2D* pmi0 = &mi0;
    
    if (::fabs(hj02 - hi02) < E_EPSILON*E_EPSILON) return false; //same metric so nothing to smooth
    if (hj02 < hi02)
    {
      std::swap(Ni,Nj);
      std::swap(hi02, hj02);
      pmi0 = &mj0;
    }
    
    // (Ni,mi) is the finer now, Nj is the one to smooth
    
    if (Nj <= N0) return false; //do not touch the hard nodes

    E_Float dij = lengthEval(Ni, *pmi0, Nj, *pmi0);
    
    E_Float hs2 = hi02 * (1. + (gr-1.) * dij) * (1. + (gr-1.) * dij); // extrapolated h at Nj with respect to growth ratio

    metric_reduce(Nj, NiNj, hj02, hs2);
    
    return true;
  }

  ///
  template <typename T>
  void
    VarMetric<T>::__init_refine_points
    (K_FLD::FloatArray& pos, size_type Ni, size_type Nj, E_Float threshold, 
    std::vector<std::pair<E_Float, size_type> >& length_to_points, std::vector<size_type>& tmpNodes)
  {
    size_type dim(pos.rows());
    
    E_Float newP[3];
    
    K_FLD::FloatArray::const_iterator it1(pos.col(Ni)), it2(pos.col(Nj));
    for (size_type i = 0; i < dim; ++i) // middle point.
      newP[i] = 0.5 * (*(it1++) + *(it2++));

    pos.pushBack(newP, newP+dim);
    size_type N = pos.cols()-1;
    
    E_Float   d = ::sqrt(K_FUNC::sqrDistance(pos.col(Ni), pos.col(Nj), pos.rows()));
    
    // set the metric as a large valid value : half of the initial edge length
    E_Float factor = 0.5;
    T m(factor*d);
    setMetric(N, m);
           
    length_to_points.push_back(std::make_pair(-d, N));
 
  }
  
  ///
  template <typename T>
  void
  VarMetric<T>::__compute_refine_points
  (K_FLD::FloatArray& pos, size_type Ni, size_type Nj, E_Float threshold, 
   std::vector<std::pair<E_Float, size_type> >& length_to_points, std::vector<size_type>& tmpNodes)
  {
    tmpNodes.clear();
    E_Float   d = length(Ni, Nj, threshold, tmpNodes);
    size_type n = std::max(size_type(d), 1);

    if ((n * (n + 1)) < (d * d))
      ++n;

    if (n == 1) // Saturated
      return;

    size_type ni = 0, ii = 0, dim(pos.rows()), nb_nodes;
    E_Float l = 0.;
    size_type Nstart = Ni, Nk = Ni, Nl = Ni;
    E_Float x = 1.;//fixme
    std::vector<std::pair<E_Float, size_type> > length_to_point;
    length_to_point.push_back(std::make_pair(0., Ni));
    E_Float* pNi = pos.col(Ni);

    nb_nodes = (size_type)tmpNodes.size();
    for (size_type i = 0; i < nb_nodes; ++i)
    {
      Nk = tmpNodes[i];
      x = K_FUNC::sqrDistance(pNi, pos.col(Nk), dim);
      length_to_point.push_back(std::make_pair(x, Nk));
    }
    length_to_point.push_back(std::make_pair(K_CONST::E_MAX_FLOAT, Nj));

    std::sort(length_to_point.begin(), length_to_point.end());

    E_Float newP[2], r, xx;
    while (ni++ < n)
    {
      l = 0.;
      while (l < 1.)
      {
        Nk = (l == 0.) ? Nstart : length_to_point[ii].second;
        Nl = length_to_point[++ii].second;
        if (Nl == Nj)
          break;
        x = lengthEval(Nk, /*metric*/(*this)[Nk], Nl, /*metric*/(*this)[Nl]);
        if ((l+x) < 1.)
          l += x;
        else
          break;
      }

      if (Nl == Nj)
        break;

      r = ((1 - l)/x);
      
      K_FUNC::diff<2>(pos.col(Nl), pos.col(Nk), newP);

      for (size_type j = 0; j < 2; ++j)
      {
        xx = pos(j, Nk);
        newP[j] = xx + r * newP[j];
      }

      pos.pushBack(newP, newP+2);
      Nstart = pos.cols()-1;

      computeMetric(Nstart, Nk, Nl, r);

      length_to_points.push_back(std::make_pair(-n, Nstart));
    }
  }




} // End namespace DELAUNAY

#endif /* DELAUNAY_METRIC_H_ */

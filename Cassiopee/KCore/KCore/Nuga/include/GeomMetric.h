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

#ifndef __GENERATOR_GEOM_METRIC_H__
#define __GENERATOR_GEOM_METRIC_H__

#include "Metric.h"

namespace DELAUNAY
{

  template <typename T, typename SurfaceType>
  class GeomMetric : public VarMetric<T>
  {
  public:

    enum GMmode
    {
      ISO_CST, ///< A constant size is specified to mesh the surface.
      ISO_RHO, ///< local minimum curvature radius is used to compute the metric.
      ANISO    ///< both local principal curvature radii are used to compute the metric.
    };

    typedef  VarMetric<T>     parent_type;
    typedef  NUGA::size_type  size_type;

  public:

    /// Constructor for ISO mode : mesh size is specified.
    GeomMetric(K_FLD::FloatArray& pos, const SurfaceType& surface, E_Float h0, E_Float gr)
      : parent_type (pos, h0, h0),
      _mode(ISO_CST), _surface(surface), _hmax2(h0*h0), _unbounded_h(false),
      _h0(h0), _chordal_error(1.), _gr(gr)
    {}

    /// Constructor for adaptive mode (ISO_RHO, ANISO) : mode and relative chordal error is specified.
    GeomMetric(K_FLD::FloatArray& pos, const SurfaceType& surface, GMmode mode, 
               E_Float chordal_error, E_Float hmin, E_Float hmax, E_Float gr)
      : parent_type (pos, hmin, hmax), _mode(mode), _surface(surface),
      _h0(NUGA::FLOAT_MAX), _chordal_error(chordal_error),
      _alpha2(4. * chordal_error*(2. - chordal_error)), _gr(gr), _unbounded_h(false)  
    {}
    
    void set_pos2D(const K_FLD::FloatArray& pos2D){_pos2D = &pos2D;}

    ~GeomMetric(void){}

    inline virtual void computeMetric(size_type N, size_type Ni, size_type Nj, E_Float /*dummy*/);

    void init_metric
      (const K_FLD::FloatArray& metric, K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
       const std::vector<E_Int>& hard_nodes);
    
    void reduce_by_checking_real_space(size_type Ni);

  private:
    void __update_boundary_metric_with_surface(const K_FLD::IntArray& connectB);

    void __compute_1st_fundamental_form(size_type N0, E_Float& E, E_Float& F, E_Float& G);

    void __computeMetric(size_type N, K_FLD::FloatArray& Mout, E_Float hmax2);

  private:
    GMmode             _mode;
    const SurfaceType& _surface;
    E_Float            _hmax2;
    E_Float            _h0;
    E_Float            _chordal_error;
    E_Float            _alpha2;
    E_Float            _gr;
    const K_FLD::FloatArray* _pos2D; //hack to avoid to pass a dummy argument for Metric::init_metric as pos2D is only required for GeomMetric
    
    //T _boundary_metric_max;
    //E_Float _humax2, _hvmax2;
    bool _unbounded_h;
  };

  ///
  template <typename T, typename SurfaceType>
  void
  GeomMetric<T, SurfaceType>::init_metric
  (const K_FLD::FloatArray& metric, K_FLD::FloatArray& pos3D, const K_FLD::IntArray& connectB,
   const std::vector<E_Int>& hard_nodes)
  {
    std::vector<E_Int> BNodes;
    //
    connectB.uniqueVals(BNodes);
    E_Int maxID = *std::max_element(ALL(BNodes));
    if (!hard_nodes.empty())
      maxID = std::max(maxID, *std::max_element(ALL(hard_nodes)));

    T m; // invalid by default
    parent_type::_field.resize(parent_type::_pos->cols(), m);

    E_Int max = std::min(maxID+1, metric.cols());

    // Set the input metric
    if (metric.rows() >= 3)
    {
      for (E_Int i = 0; i < max; ++i)
      {
        m = metric.col(i);
        
        if (parent_type::isValidMetric(m))
          parent_type::_field[i] = m;
      }
    }   

    // Compute the isometric metric based on the boundary.
    std::vector<E_Float> m_iso;
    E_Float hmin1, hmax1;
    MeshUtils1D::compute_iso_metric(pos3D, connectB, hard_nodes, m_iso, hmin1, hmax1);

    // Set _hmin and _hmax if not done nor correct.
     if (parent_type::_hmax <= 0.)
     {
       parent_type::_hmax = hmax1 * connectB.cols();
       _unbounded_h = true;
     }
     if ((parent_type::_hmin <= 0.) || (parent_type::_hmin > parent_type::_hmax) || (parent_type::_hmin == NUGA::FLOAT_MAX) )
       parent_type::_hmin = std::min(hmin1, parent_type::_hmax);
 
     if (parent_type::_hmax < hmax1)
     {
       for (size_t i=0; i < m_iso.size(); ++i)
         m_iso[i] = std::min(m_iso[i], parent_type::_hmax);
     }
     if (parent_type::_hmin > hmin1)
     {
       for (size_t i=0; i < m_iso.size(); ++i)
         m_iso[i] = std::max(m_iso[i], parent_type::_hmin);
     }
     
    _hmax2 = (parent_type::_hmax != NUGA::FLOAT_MAX) ? parent_type::_hmax * parent_type::_hmax : NUGA::FLOAT_MAX; //fixme : important to be done before __update_boundary_metric_with_surface

    // Set _metric by converting m_iso to an aniso type metric.
    // Only relevant in aniso case (in iso just assign it) //fixme
    this->convertIsoToAniso(m_iso, parent_type::_field);
    
    // Transform the metric in the parameter space.
    E_Float E,F,G;
    for (size_t Ni = 0; Ni < parent_type::_field.size(); ++Ni)
    {
      __compute_1st_fundamental_form(Ni, E, F, G);
      
      // Do the transform
      parent_type::_field[Ni][0] *= E;
      //parent_type::_field[Ni][1] *= F; //here m12 is 0
      parent_type::_field[Ni][2] *= G;
      
#ifdef DEBUG_METRIC
      assert (parent_type::isValidMetric(parent_type::_field[Ni]));
#endif
  }

#ifdef DEBUG_METRIC
      {
        //parent_type::draw_ellipse_field("ellipse_raw_param_space.mesh", *_pos2D, connectB);
      }
#endif
  
  
//    // LIMITERS
//    
//    // 1. compute the max ellipse base on the max du and dv over the boundary
//    E_Float hu_max(0.), hv_max(0.);
//    E_Int nbe = connectB.cols();
//    const K_FLD::FloatArray &crd = *_pos2D;
//    for (E_Int i=0; i < nbe; ++i)
//    {
//      E_Float du = ::fabs(crd(0,connectB(0,i)) - crd(0,connectB(1,i)));
//      E_Float dv = ::fabs(crd(1,connectB(0,i)) - crd(1,connectB(1,i)));
//      
//      hu_max = std::max(hu_max, du);
//      hv_max = std::max(hv_max, dv);
//    }
//    
////    E_Float huvmax = std::max(hu_max, hv_max);
////
////    _boundary_metric_max[1] = 0.;
//    _humax2 = hu_max*hu_max;
//    _hvmax2 = hv_max*hv_max;
//    E_Float mu_min  =  1./ _humax2;
//    E_Float mv_min = 1./ _hvmax2;
//    
//    for (size_t Ni = 0; Ni < parent_type::_field.size(); ++Ni)
//    {
//      parent_type::_field[Ni][0] = std::max(parent_type::_field[Ni][0], mu_min);
//      parent_type::_field[Ni][2] = std::max(parent_type::_field[Ni][2], mv_min);
//    }
//    
//#ifdef DEBUG_METRIC
//    parent_type::draw_ellipse_field("ellipse_cap_param_space.mesh", *_pos2D, connectB); //fixme : wrong pos, shoule be pos2D
//#endif
//    
//    // 2. optimal reduce : ensure that the ellipses intersection with the bounday does not contain inner nodes.
//    //    APPROX : check only connected nodes => not true for erratic boundary
//    for (E_Int i=0; i < nbe; ++i)
//    {
//      const E_Int& Ni = connectB(0,i);
//      const E_Int& Nj = connectB(1,i);
//      
//      E_Float NiNj[2];
//      NUGA::diff<2>(crd.col(Nj),crd.col(Ni), NiNj);
//      E_Float dij2 = NUGA::normalize<2>(NiNj);
//      dij2 *= dij2;
//      
//      E_Float hi02 = parent_type::get_h2_along_dir(Ni, NiNj); // trace on NiNj of the ellipse centered at Ni
//      
//      if (hi02> dij2) //need to be reduced along that direction : the ellipse will pass through Nj after transfo
//        parent_type::metric_reduce(Ni, NiNj, hi02, dij2);
//
//      E_Float hj02 = parent_type::get_h2_along_dir(Nj, NiNj); // trace on NiNj of the ellipse centered at Nj
//      
//      if (hj02 > dij2) //need to be reduced along that direction
//        parent_type::metric_reduce(Nj, NiNj, hj02, dij2);
//    }
//
//#ifdef DEBUG_METRIC
//    parent_type::draw_ellipse_field("ellipse_opt_reduce_param_space.mesh", *_pos2D, connectB); //fixme : wrong pos, shoule be pos2D
//#endif
//    
//#ifdef DEBUG_METRIC
//    std:: cout << "Is metric valid ? " << parent_type::is_valid() << std::endl;
//#endif
//    
//    // 3. real space checking
//    if (!_unbounded_h)
//    {
//      for (size_t Ni = 0; Ni < parent_type::_field.size(); ++Ni)
//        reduce_by_checking_real_space(Ni);
//    }
// 
//#ifdef DEBUG_METRIC
//    parent_type::draw_ellipse_field("ellipse_real_reduce_param_space.mesh", *_pos2D, connectB); //fixme : wrong pos, shoule be pos2D
//#endif

    // Now take the user metric into account at each node when it's valid.
    // Only the first row of the input matrix is taken into account.
    //compute_intersection(_metric, metric);
    
    //__update_boundary_metric_with_surface(connectB);

#ifdef DEBUG_METRIC
    std:: cout << "Is metric valid ? " << parent_type::is_valid() << std::endl;
#endif
  }

  template <typename T, typename SurfaceType>
  void GeomMetric<T, SurfaceType>::__computeMetric(size_type N0, K_FLD::FloatArray& Mout, E_Float hmax2)
  {
    Mout.resize(2,2);
    Mout(0,0) = Mout(1,1) = 1./(parent_type::_hmin*parent_type::_hmin);
    Mout(0,1) = Mout(1,0) = 0.; // Iso by default.

    E_Float* Q = parent_type::_pos->col(N0);
    E_Float  u = Q[0];
    E_Float  v = Q[1];
    E_Float  E, F, G; // Tangential plane's first fundamental form coefficients.
    E_Float  L, M, N; // Tangential plane's second fundamental form coefficients.
    //E_Float  dU1[3], dU2[3], dV1[3], dV2[3], dUV[3];
    E_Float dA[15];
    E_Float* dU1 = dA; E_Float* dV1 = dA+3;
    E_Float* dU2 = dA+6; E_Float* dV2 = dA+9; E_Float* dUV = dA+12;
    E_Float n[3];
    E_Float hmin2 = parent_type::_hmin*parent_type::_hmin;
    
    // Ask the surface.
    //_surface.DU1(u, v, dU1); // First base vector of the tangential plane.
    //_surface.DU2(u, v, dU2);
    //_surface.DV1(u, v, dV1); // Second base vector of the tangential plane.
    //_surface.DV2(u, v, dV2);
    //_surface.DUV(u, v, dUV);
    _surface.DAUV2(u, v, dA);
    
    NUGA::crossProduct<3> (dU1, dV1, n); // Normal to the tangential plane.
    E_Float ln = NUGA::normalize<3>(n);
    
    bool singular = (::fabs(ln) < EPSILON); //undefined plane : dU1 and dV2 are colinear !

    if (!singular)
    {
      // First form. = metric tensor
      E = NUGA::sqrNorm<3> (dU1);
      F = NUGA::dot<3> (dU1, dV1);
      G = NUGA::sqrNorm<3> (dV1);
      
      singular = ((E < EPSILON || G < EPSILON));
    }
    
    if (singular) return;

    if (_mode == ISO_CST) // impose hmin
    {
      Mout(0,0) *= E;
      Mout(1,1) *= G;
      Mout(1,0) = Mout(0,1) *= F;
      return;
    }

    // Second form. = shape tensor
    L = NUGA::dot<3>(n, dU2);
    M = NUGA::dot<3>(n, dUV);
    N = NUGA::dot<3>(n, dV2);

    bool locally_iso = ((::fabs((F*L)-(E*M)) < EPSILON) && 
                        (::fabs((G*L)-(E*N)) < EPSILON) && 
                        (::fabs((G*M)-(F*N)) < EPSILON));

    if (locally_iso)
    {
      E_Float R2 = E/L; // pourquoi cette valeur ? 
      R2 *= R2;
      E_Float h2 = std::min(hmax2, _alpha2*R2); // alpha2 prend en compte hausd
      h2 = std::max(h2, hmin2); //fixme!!
      h2 = 1./h2;

      Mout(0,0) = E*h2;
      Mout(1,1) = G*h2;
      Mout(1,0) = Mout(0,1) = F*h2;
      
      return;
    }

    K_FLD::FloatArray M1(2,2), M2(2,2), NN(2,2);

    // matrice du tenseur metrique
    M1(0,0) = E;
    M1(0,1) = M1(1,0) = F;
    M1(1,1) = G;
    // matrice du tenseur de forme
    M2(0,0) = L;
    M2(0,1) = M2(1,0) = M;
    M2(1,1) = N;

    // vecteurs propres de inv(M1) * M2
    E_Float V1[2], V2[2];
    
    // Simultaneous reduction. NN = inv(M1) * M2
    K_LINEAR::DelaunayMath::simultaneous_reduction(M1, M2, V1, V2);

    // CB : gaussian curvature pour verification (K = K1*K2)
    //E_Float K = (L*N - M*M) / (E*G - M*M);

    // principal curvatures
    E_Float rho1 = (E*V1[0]*V1[0] + 2*F*V1[0]*V1[1] + G*V1[1]*V1[1]) /
      (L*V1[0]*V1[0] + 2*M*V1[0]*V1[1] + N*V1[1]*V1[1]);
    E_Float rho2 = (E*V2[0]*V2[0] + 2*F*V2[0]*V2[1] + G*V2[1]*V2[1]) /
      (L*V2[0]*V2[0] + 2*M*V2[0]*V2[1] + N*V2[1]*V2[1]);

    rho1 = (rho1 < 0.) ? -rho1 : rho1; // valeur absolue
    rho2 = (rho2 < 0.) ? -rho2 : rho2;

    if (rho2 < rho1)
    {
      std::swap(rho1, rho2);
      std::swap(V1, V2);
    }
    // rho1, V1 est maintenant le plus petit rayon de courbure

    E_Float rho1_2 = rho1 * rho1;

    if (_mode == ISO_RHO) // use min curvature in all directions + impose hmin
    {
      E_Float h2 = std::min(hmax2, _alpha2*rho1_2);
      h2 = std::max(h2, hmin2);
      h2 = 1./h2;

      Mout(0,0) = E*h2;
      Mout(1,1) = G*h2;
      Mout(1,0) = Mout(0,1) = F*h2;
      
#ifdef DEBUG_METRIC
      {
      //assert (parent_type::isValidMetric(Mout));
      
      const E_Float& a11 = Mout(0,0);
      const E_Float& a12 = Mout(0,1);
      const E_Float& a22 = Mout(1,1);
      E_Float det = (a11*a22) - (a12*a12);
      assert ((a11 > 0.) && (a22 > 0.) && (det > 0.)); //i.e. isValidMetric
      }
#endif

      return;
    }

    E_Float rho2_2 = rho2 * rho2; // plus grand rayon de courbure

    rho1_2 = std::min(hmax2/_alpha2, rho1_2); // impose hmax et hmin
    rho2_2 = std::min(hmax2/_alpha2, rho2_2);
    rho1_2 = std::max(hmin2/_alpha2, rho1_2); //fixme
    rho2_2 = std::max(hmin2/_alpha2, rho2_2); //fixme

    E_Float q = 1. - sqrt(rho1_2/rho2_2);
    E_Float rho2_2c = rho2_2*(1.-q*q); // reduce the largest according to q

    E_Float h1_2 = _alpha2*rho1_2;
    E_Float h2_2 = _alpha2*rho2_2c;

    // Use an interpolation for the third metric
    //  E_Float k = 0.;
    //  E_Float maxh = std::max(h1_2, h2_2);
    E_Float minh = std::min(h1_2, h2_2);
    E_Float h3_2 = minh;//(((1.-k)*minh) + (k * maxh));

    K_FLD::FloatArray Base(3,3, 0.);
    E_Float zero = 0.;
    Mout.resize(3,3,&zero);

    Base(0,0) = V1[0];
    Base(1,0) = V1[1];
    Base(0,1) = V2[0];
    Base(1,1) = V2[1];
    Base(2,2) = 1.;

    Mout(0,0) = 1. / h1_2;
    Mout(1,1) = 1. / h2_2;
    Mout(2,2) = 1. / h3_2;

    // Express the metric in (dU1, dV1).
    Mout = Mout * Base;
    Mout = Base.transpose() * Mout;

    Base(0,0) = dU1[0];
    Base(1,0) = dU1[1];
    Base(2,0) = dU1[2];
    Base(0,1) = dV1[0];
    Base(1,1) = dV1[1];
    Base(2,1) = dV1[2];
    Base(0,2) = n[0];
    Base(1,2) = n[1];
    Base(2,2) = n[2];

    // Express the metric in the canonical base.
    Mout = Mout * Base;
    Mout = Base.transpose() * Mout;
  }

  // N0 : calul la metrique au point N0
  // Ni, Nj: points interpolants
  template <typename T, typename SurfaceType>
  void GeomMetric<T, SurfaceType>::computeMetric
    (size_type N0, size_type /*dummy*/Ni, size_type /*dummy*/Nj, E_Float /*dummy*/r)
  {
    
    if ((E_Int)parent_type::_field.size() > N0)//fixme : work around to avoid to set more than once
      return;

    K_FLD::FloatArray M(2,2);
    __computeMetric(N0, M, _hmax2);

    T m;
    m[0] = M(0,0);
    m[1] = M(1,0);
    m[2] = M(1,1);

    if (!parent_type::isValidMetric(m)) // e.g. surface is locally planar
      m = parent_type::_interpol->interpolate(parent_type::_field[Ni], parent_type::_field[Nj], r);

    parent_type::setMetric(N0, m);
    
    if (_gr > 1.)
    {
      parent_type::smooth(Ni, N0, _gr, Ni);
      parent_type::smooth(Nj, N0, _gr, Nj);
    }
    
    assert (parent_type::isValidMetric(parent_type::_field[N0]));
    
    //parent_type::draw_ellipse("ellipse_step1.mesh", *_pos2D, N0);

    //if (!_unbounded_h)
    //  reduce_by_checking_real_space(N0);
    //parent_type::draw_ellipse("ellipse_real.mesh", *_pos2D, N0);
    
    assert (parent_type::isValidMetric(parent_type::_field[N0]));
  }

  template <typename T, typename SurfaceType>
  void
  GeomMetric<T, SurfaceType>::__compute_1st_fundamental_form(size_type N0, E_Float& E, E_Float& F, E_Float& G)
  {
    E_Float* Q = parent_type::_pos->col(N0);
    E_Float  u = Q[0];
    E_Float  v = Q[1];
    //E_Float  dU1[3], dV1[3];
    E_Float dA[6];
    E_Float* dU1 = dA; E_Float* dV1 = dA+3;
    
    // Ask the surface.
    //_surface.DU1(u,v, dU1); // First base vector of the tangential plane.
    //_surface.DV1(u,v, dV1);
    _surface.DAUV1(u,v, dA);

    // First form.
    E = NUGA::sqrNorm<3> (dU1);
    F = NUGA::dot<3> (dU1, dV1);
    G = NUGA::sqrNorm<3> (dV1);
  }

  template <typename T, typename SurfaceType>
  void
    GeomMetric<T, SurfaceType>::__update_boundary_metric_with_surface
    (const K_FLD::IntArray& connectB)
  {
    std::vector<E_Int> Bnodes;
    connectB.uniqueVals(Bnodes);

    K_FLD::FloatArray M1(2,2), M2(2,2), ID(2,2);
    E_Int Ni;
    std::vector<E_Float> s(3);
    
    for (size_t i = 0; i < Bnodes.size(); ++i)
    {
      Ni = Bnodes[i];
      M1(0,0) = parent_type::_metric[Ni][0];
      M1(1,0) = M1(0,1) = parent_type::_metric[Ni][1];
      M1(1,1) = parent_type::_metric[Ni][2];

      __computeMetric(Ni, M2, 1./parent_type::_metric[Ni][0]); // Metric obtained by the geometry (reduced by the edge length.fixme)

      //fixme : should update somehow hmin and hmax here. ?
      //K_LINEAR::DelaunayMath::intersect(M1, M2, I);
      //std::cout << I << std::endl;

      s[0] = M2(0,0);
      s[1] = M2(1,0);
      s[2] = M2(1,1);
      parent_type::_metric[i] = Aniso2D(&s[0]);
    }
  }

  template <typename T, typename SurfaceType>
  void
  GeomMetric<T, SurfaceType>::reduce_by_checking_real_space(size_type Ni)
  {
    if (_hmax2 == NUGA::FLOAT_MAX || _hmax2 <= 0.) return;
    
    // Diagonalization
    E_Float lambda0, lambda1, v0[2], v1[2];
    K_LINEAR::DelaunayMath::eigen_vectors(parent_type::_field[Ni][0], parent_type::_field[Ni][2], parent_type::_field[Ni][1], lambda0, lambda1, v0, v1);
    K_FLD::FloatArray D(2,2, 0.);
    D(0,0) = lambda0;
    D(1,1) = lambda1;
    
    NUGA::normalize<2>(v0);
    NUGA::normalize<2>(v1);
    
    // transformation Matrix : BD : Main axis -> (i,j)
    K_FLD::FloatArray P(2,2, 0.), tP(2,2,0.);
    tP(0,0) = P(0,0) = v0[0];
    tP(0,1) = P(1,0) = v0[1];
    tP(1,0) = P(0,1) = v1[0];
    tP(1,1) = P(1,1) = v1[1];
    
    
    E_Float P0[2], P1[2];
    
    E_Float h0 = sqrt(1./lambda0);
    E_Float h1 = sqrt(1./lambda1);
    //E_Float Pix = parent_type::_pos(0,Ni);
    //E_Float Piy = parent_type::_pos(1,Ni);

    // P0 = Pi +/- h0*v0
    NUGA::sum<3>(1., parent_type::_pos->col(Ni), h0, v0, P0);
    bool P0_is_inside = _surface.in_bounds(P0[0], P0[1]);
    if (!P0_is_inside)
    {
      NUGA::sum<3>(1., parent_type::_pos->col(Ni), -h0, v0, P0); //try with the opposite point
      P0_is_inside = _surface.in_bounds(P0[0], P0[1]);
    }

    // P1 = Pi +/- h1*v1
    NUGA::sum<3>(1., parent_type::_pos->col(Ni), h1, v1, P1);
    bool P1_is_inside = _surface.in_bounds(P1[0], P1[1]);
    if (!P1_is_inside)
    {
      NUGA::sum<3>(1., parent_type::_pos->col(Ni), -h1, v1, P1); //try with the opposite point
      P1_is_inside = _surface.in_bounds(P1[0], P1[1]);
    }
    
    if (!P0_is_inside && !P1_is_inside) return; 
    
    // Transform in the real space
    E_Float P0r[3], Pi[3], P1r[3];
    _surface.point(parent_type::_pos(0,Ni), parent_type::_pos(1,Ni), Pi); //center of the ellipse
    
    E_Float k02 = 1.;
    if (P0_is_inside)
    {
      _surface.point(P0[0], P0[1], P0r);
      E_Float dPiP02 = NUGA::sqrDistance(Pi, P0r, 3);
      k02 = dPiP02 / _hmax2;
    }
    
    E_Float k12 = 1.;
    if (P1_is_inside)
    {
      _surface.point(P1[0], P1[1], P1r);
      E_Float dPiP12 = NUGA::sqrDistance(Pi, P1r, 3);
      k12 =  dPiP12 / _hmax2;
    }

    if (k02 <= 1. && k12 <= 1.) return; //smaller than hmax
    
    E_Float KMAX2 = (1. + sqrt(2.));
    KMAX2 *= KMAX2;
        
    if (k02 > 1.1 && k02 < KMAX2) D(0,0) *= k02;
    if (k12 > 1.1 && k12 < KMAX2) D(1,1) *= k12;
        
//    if (k12 > KMAX2 || k02 > KMAX2)
//    {
//      std::cout << "wthfck Ni : " << Ni << " k02 :" << k02 << " k12 : " << k12 << std::endl;
//      std::ostringstream o;
//      o << "error_real_" << Ni << ".mesh";
//      //parent_type::draw_ellipse(o.str().c_str(), parent_type::_pos, Ni);
//    }
    
    K_FLD::FloatArray M = P * D * tP;
    
    parent_type::_field[Ni][0] = M(0,0);
    parent_type::_field[Ni][1] = M(1,0);
    parent_type::_field[Ni][2] = M(1,1);
  }
}

#endif

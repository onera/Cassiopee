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

    typedef     VarMetric<T>            parent_type;
    typedef     K_CONT_DEF::size_type   size_type;

  public:

    /// Constructor for ISO mode : mesh size is specified.
    GeomMetric(K_FLD::FloatArray& pos, const SurfaceType& surface, E_Float h0)
      :parent_type (pos, h0, h0),
      _mode(ISO_CST), _surface(surface), _hmax2(h0*h0),
      _h0(h0), _chordal_error(1.)
    {}

    /// Constructor for adaptive mode (ISO_RHO, ANISO) : mode and relative chordal error is specified.
    GeomMetric(K_FLD::FloatArray& pos, const SurfaceType& surface, GMmode mode, 
               E_Float chordal_error, E_Float hmin, E_Float hmax)
      :parent_type (pos, hmin, hmax), _mode(mode), _surface(surface),
      _h0(K_CONST::E_MAX_FLOAT), _chordal_error(chordal_error),
      _alpha2(4. * chordal_error*(2. - chordal_error))  
    {}

    ~GeomMetric(void){}

    inline virtual void computeMetric(size_type N, size_type Ni, size_type Nj, E_Float /*dummy*/);

    inline virtual void setMetric(E_Int N, const T& m);

    void init_metric
      (const K_FLD::FloatArray& metric, K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
       const std::vector<E_Int>& hard_nodes);
    
    void __init_refine_points
    (K_FLD::FloatArray& pos, size_type Ni, size_type Nj, E_Float threshold,
    std::vector<std::pair<E_Float, size_type> >& length_to_points, std::vector<size_type>& tmpNodes);

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
  };

  ///
  template <typename T, typename SurfaceType>
  void
  GeomMetric<T, SurfaceType>::init_metric
  (const K_FLD::FloatArray& metric, K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
   const std::vector<E_Int>& hard_nodes)
  {
    std::vector<E_Int> BNodes;
    //
    connectB.uniqueVals(BNodes);
    E_Int maxID = *std::max_element(ALL(BNodes));
    if (!hard_nodes.empty())
      maxID = std::max(maxID, *std::max_element(ALL(hard_nodes)));

    T m;// invalid by default
    parent_type::_field.resize(parent_type::_pos.cols(), m);

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
    MeshUtils1D::compute_iso_metric(pos, connectB, hard_nodes, m_iso, hmin1, hmax1);
    
    // Set _hmin and _hmax if not done nor correct.
    if ((parent_type::_hmin <= 0.) || (parent_type::_hmin == K_CONST::E_MAX_FLOAT)) // If not (or wrongly) set by the user.
      parent_type::_hmin = /*0.25 **/ hmin1;
    if ((parent_type::_hmax <= 0.) || (parent_type::_hmax == K_CONST::E_MAX_FLOAT)) // If not (or wrongly) set by the user.
      parent_type::_hmax = hmax1;
    
    if (parent_type::_hmax != hmax1) //i.e. user-defined
    {
      for (size_t i=0; i < m_iso.size(); ++i)
        m_iso[i] = std::min(m_iso[i], parent_type::_hmax);
    }

    if (parent_type::_hmin > hmin1)
    {
      for (size_t i=0; i < m_iso.size(); ++i)
        m_iso[i] = std::max(m_iso[i], parent_type::_hmin);
    }
      
    _hmax2 = (parent_type::_hmax != K_CONST::E_MAX_FLOAT) ? parent_type::_hmax * parent_type::_hmax : K_CONST::E_MAX_FLOAT; //fixme : important to be done before __update_boundary_metric_with_surface

    // Set _metric by converting m_iso to an aniso type metric.
    // Only relevant in aniso case (in iso just assign it) //fixme
    this->convertIsoToAniso(m_iso, parent_type::_field);
    
#ifdef DEBUG_METRIC
      {
        parent_type::draw_ellipse_field("ellipse_init_metric.mesh", pos, connectB);
      }
#endif

    // Transform the metric in the parameter space.
    E_Float E,F,G;
    for (size_t Ni = 0; Ni < parent_type::_field.size(); ++Ni)
    {
      __compute_1st_fundamental_form(Ni, E, F, G);

      if ((E < E_EPSILON) || (G < E_EPSILON)) //singular point.
        continue;

      // Do the transform if and only if the ellispse is contained in the initial iso circle
      if (E > (1. / (parent_type::_field[Ni][0] * hmax1 * hmax1) ) )
        parent_type::_field[Ni][0] *= E;
      
      //parent_type::_field[Ni][1] *= F; //here m12 is 0
      
      // Do the transform if and only if the ellispse is contained in the initial iso circle
      if (G > (1. / (parent_type::_field[Ni][2] * hmax1 * hmax1) ) )
        parent_type::_field[Ni][2] *= G;
      
#ifdef DEBUG_METRIC
      assert (isValidMetric(parent_type::_field[Ni]));
#endif
    }
    
#ifdef DEBUG_METRIC
      {
        parent_type::draw_ellipse_field("ellipse_param_space.mesh", pos, connectB);
      }
#endif

    // Now take the user metric into account at each node when it's valid.
    // Only the first row of the input matrix is taken into account.
    //compute_intersection(_metric, metric);
    
/*
    __update_boundary_metric_with_surface(connectB);
    */
  }

  template <typename T, typename SurfaceType>
  void GeomMetric<T, SurfaceType>::__computeMetric(size_type N0, K_FLD::FloatArray& Mout, E_Float hmax2)
  {
    Mout.resize(2,2);
    Mout(0,0) = Mout(1,1) = 1./(parent_type::_hmin*parent_type::_hmin);
    Mout(0,1) = Mout(1,0) = 0.; // Iso by default.

    E_Float* Q = parent_type::_pos.col(N0);
    E_Float  u = Q[0];
    E_Float  v = Q[1];
    E_Float  E, F, G; // Tangential plane's first fundamental form coefficients.
    E_Float  L, M, N; // Tangential plane's second fundamental form coefficients.
    E_Float  dU1[3], dU2[3], dV1[3], dV2[3], dUV[3], n[3];
    E_Float  hmin2 = parent_type::_hmin*parent_type::_hmin;
    
    // Ask the surface.
    _surface.DU1(u,v, dU1); // First base vector of the tangential plane.
    _surface.DU2(u,v, dU2); // Second base vector of the tangential plane.
    _surface.DV1(u,v, dV1);
    _surface.DV2(u,v, dV2);
    _surface.DUV(u,v, dUV);

    K_FUNC::crossProduct<3> (dU1, dV1, n); // Normal to the tangential plane.
    E_Float l = K_FUNC::normalize<3>(n);
    
    bool singular = (::fabs(l) < E_EPSILON); //undefined plane : dU1 and dV2 are colinear !

    if (!singular)
    {
      // First form.
      E = K_FUNC::sqrNorm<3> (dU1);
      F = K_FUNC::dot<3> (dU1, dV1);
      G = K_FUNC::sqrNorm<3> (dV1);
      
      singular = ((E < E_EPSILON || G < E_EPSILON));
    }
    
    if (singular) return;

    if (_mode == ISO_CST) // but not singular
    {
      Mout(0,0) *= E;
      Mout(1,1) *= G;
      Mout(1,0) = Mout(0,1) *= F;
      
      return;
    }

    // Second form.
    L = K_FUNC::dot<3>(n, dU2);
    M = K_FUNC::dot<3>(n, dUV);
    N = K_FUNC::dot<3>(n, dV2);

    bool locally_iso = ((::fabs((F*L)-(E*M)) < E_EPSILON) && 
                        (::fabs((G*L)-(E*N)) < E_EPSILON) && 
                        (::fabs((G*M)-(F*N)) < E_EPSILON));

    if (locally_iso)
    {
      E_Float R2 = E/L;
      R2 *= R2;
      E_Float h2 = std::min( hmax2, _alpha2*R2);
      h2 = std::max( h2, hmin2);//fixme!!
      h2 = 1./h2;

      Mout(0,0) = E*h2;
      Mout(1,1) = G*h2;
      Mout(1,0) = Mout(0,1) = F*h2;
      
      return;
    }

    K_FLD::FloatArray M1(2,2), M2(2,2), NN(2,2);

    M1(0,0) = E;
    M1(0,1) = M1(1,0) = F;
    M1(1,1) = G;
    M2(0,0) = L;
    M2(0,1) = M2(1,0) = M;
    M2(1,1) = N;

    E_Float v1[2], v2[2];
    E_Float *V1(v1), *V2(v2);

    // Simultaneous reduction. NN = inv(M1) * M2
    K_LINEAR::DelaunayMath::simultaneous_reduction (M1, M2, V1, V2);

    E_Float rho1 = (E*V1[0]*V1[0] + 2*F*V1[0]*V1[1] + G*V1[1]*V1[1]) /
      (L*V1[0]*V1[0] + 2*M*V1[0]*V1[1] + N*V1[1]*V1[1]);
    E_Float rho2 = (E*V2[0]*V2[0] + 2*F*V2[0]*V2[1] + G*V2[1]*V2[1]) /
      (L*V2[0]*V2[0] + 2*M*V2[0]*V2[1] + N*V2[1]*V2[1]);

    rho1 = (rho1 < 0.) ? -rho1 : rho1;
    rho2 = (rho2 < 0.) ? -rho2 : rho2;

    if (rho2 < rho1)
    {
      std::swap(rho1, rho2);
      std::swap(V1, V2);
    }

    E_Float rho1_2 = rho1 * rho1;

    if (_mode == ISO_RHO) //use min rho
    {
      E_Float h2 = std::min( hmax2, _alpha2*rho1_2);
      h2 = std::max( h2, hmin2);//fixme!!
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

    E_Float rho2_2 = rho2 * rho2;

    rho1_2 = std::min (hmax2/_alpha2, rho1_2);
    rho2_2 = std::min (hmax2/_alpha2, rho2_2);
    rho1_2 = std::max (hmin2/_alpha2, rho1_2);//fixme
    rho2_2 = std::max (hmin2/_alpha2, rho2_2);//fixme

    E_Float q = 1. - ::sqrt(rho1_2/rho2_2);
    E_Float rho2_2c = rho2_2*(1.-(q*q));// reduce the largest according to q

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

  template <typename T, typename SurfaceType>
  void
    GeomMetric<T, SurfaceType>::computeMetric
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

    parent_type::setMetric(N0, m);
  }

  //fixme : implementation required ?
   template <typename T, typename SurfaceType>
  inline 
  void GeomMetric<T, SurfaceType>::setMetric(E_Int N, const T& m)
  {
    //if (isValidMetric(m)) // relates to the above work around.
    if ((E_Int)parent_type::_field.size() > N)//fixme : work around to avoid to set more than once
      return;
    parent_type::setMetric(N, m);
  }

  template <typename T, typename SurfaceType>
  void
  GeomMetric<T, SurfaceType>::__compute_1st_fundamental_form(size_type N0, E_Float& E, E_Float& F, E_Float& G)
  {
    E_Float* Q = parent_type::_pos.col(N0);
    E_Float  u = Q[0];
    E_Float  v = Q[1];
    E_Float  dU1[3], dV1[3];
    
    // Ask the surface.
    _surface.DU1(u,v, dU1); // First base vector of the tangential plane.
    _surface.DV1(u,v, dV1);
    
    // First form.
    E = K_FUNC::sqrNorm<3> (dU1);
    F = K_FUNC::dot<3> (dU1, dV1);
    G = K_FUNC::sqrNorm<3> (dV1);
  }

  template <typename T, typename SurfaceType>
  void
    GeomMetric<T, SurfaceType>::__update_boundary_metric_with_surface
    (const K_FLD::IntArray& connectB)
  {
    std::vector<E_Int> Bnodes;
    connectB.uniqueVals(Bnodes);

    K_FLD::FloatArray M1(2,2), M2(2,2), I(2,2);
    E_Int Ni ;
    std::vector<E_Float> s(3);
    
    for (size_t i = 0; i < Bnodes.size(); ++i)
    {
      Ni = Bnodes[i];
      M1(0,0) = parent_type::_metric[Ni][0];
      M1(1,0) = M1(0,1) = parent_type::_metric[Ni][1];
      M1(1,1) = parent_type::_metric[Ni][2];

      __computeMetric(Ni, M2, 1./parent_type::_metric[Ni][0]); // Metric obtained by the geometry (reduced by the edge length.fixme)

      //fixme : should update somehow hmin and hmax here. ?

//      std::cout << M1 << std::endl;
 //     std::cout << M2 << std::endl;

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
  GeomMetric<T, SurfaceType>::__init_refine_points
  (K_FLD::FloatArray& pos, size_type Ni, size_type Nj, E_Float threshold,
   std::vector<std::pair<E_Float, size_type> >& length_to_points, std::vector<size_type>& tmpNodes)
  {
    // For a Geom Mesh, it doesn't make sense to symetrize the param mesh, so
   // do regular refinement
    parent_type::__compute_refine_points(pos, Ni, Nj, threshold, length_to_points, tmpNodes);
  }

}

#endif

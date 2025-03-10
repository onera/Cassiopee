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

#ifndef NUGA_XSENSOR2_HXX
#define NUGA_XSENSOR2_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/mesh_t.hxx"
#include "V1_smoother.hxx"
#include "shell_smoother.hxx"
#include "Nuga/include/localizer.hxx"
#include "Nuga/include/collider.hxx"

namespace NUGA
{

  enum eMetricPolicy {MIN=0, MEAN=1, MAX=2, MIN_OR_MAX=3};

/// X geometric sensor 2 : recursive-collision test with surrounding metric
template <typename mesh_t> //ngu for surfacic (PGs) or ngon_t for volumic
class xsensor2 : public sensor<mesh_t, pg_smesh_t>
{  
  public:
    using parent_t = sensor<mesh_t, pg_smesh_t>;
    using sensor_input_t = pg_smesh_t;
    using output_t = typename mesh_t::output_t;
    
    //
    xsensor2(mesh_t& mesh, eSmoother smoo_type, eMetricPolicy metric_policy, E_Int itermax = -1) :
      parent_t(mesh, nullptr), _metric_policy(metric_policy), _iter_max((itermax <= 0) ? INT_MAX : itermax), _iter(0), _done(false)
    {
      // smoother
      if (smoo_type == eSmoother::V1_NEIGH)
        parent_t::_smoother = new V1_smoother<mesh_t>();
      else if (smoo_type == eSmoother::SHELL)
        parent_t::_smoother = new shell_smoother<mesh_t>();
    }

    E_Int assign_data(const sensor_input_t& data) override;

    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override;

    bool update() override;

    bool stop() override;

    virtual ~xsensor2() {}

private:
  eMetricPolicy _metric_policy;
  E_Int _iter_max, _iter;
  bool _done;
  std::map<E_Int, std::vector<E_Int>> _candidates;
  static constexpr double RTOL = 1.e-15;
  //static constexpr NUGA::eMetricType MTYPE = NUGA::ISO_MIN;
};

///
template <typename mesh_t>
E_Int xsensor2<mesh_t>::assign_data(const sensor_input_t& data)
{
  parent_t::assign_data(data);

  _candidates.clear();
  _done = false;
  _iter = 0;

  return 0;
}

///
template <typename mesh_t>
bool xsensor2<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  bool filled{ false };

  E_Int nb_cells = parent_t::_hmesh._ng.PHs.size();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  
  adap_incr.face_adap_incr.clear();
  adap_incr.cell_adap_incr.clear();

  using cell_incr_t = typename output_t::cell_incr_t;
  using face_incr_t = typename output_t::face_incr_t;

  adap_incr.cell_adap_incr.resize(nb_cells, cell_incr_t(0));
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  ph_mesh_t m1;//fixme
  m1.crd = parent_t::_hmesh._crd;
  m1.cnt = parent_t::_hmesh._ng;
  m1.get_nodal_metric2(NUGA::ISO_MAX);//fixme : just here to construct ae. ISO_MAX to avoid useless KdTree construc in build_nodal_metric2


  if (_candidates.empty()) // first time with this parent_t::_data
  { 
    pg_smesh_t& src_mesh = parent_t::_data;
    auto loc2 = src_mesh.get_localizer();
    if (loc2 == nullptr) return filled;

    std::vector<E_Int> cands;
    for (E_Int i = 0; i < nb_cells; ++i)
    {
      if (!parent_t::_hmesh._PHtree.is_enabled(i)) continue;

      auto ae1 = m1.aelement(i);
      //double Lref21 = ae1.Lref2(ae1.m_crd, MTYPE);

      cands.clear();
      loc2->get_candidates(ae1, ae1.m_crd, cands, 0, RTOL); //return as 0-based (fixme for volumic, was 1-based)
      if (cands.empty()) continue;

      bool is_x = NUGA::COLLIDE::get_colliding(ae1, src_mesh, cands, 0, RTOL, false/*all of them*/);
      if (!is_x) continue;
      if (cands.empty()) continue;

      _candidates[i] = cands;
    }

  }

  pg_smesh_t& src_mesh = parent_t::_data;
  double K = 2.;

  for (auto it : _candidates)
  {
    E_Int PHi = it.first;
    const auto& cands = it.second;

    auto ph = m1.element(PHi);

    // criterion based on ISO MIN
    if ( (_metric_policy == NUGA::MIN) || (_metric_policy == NUGA::MIN_OR_MAX) )
    {
      double phlr2 = ph.Lref2(m1.crd, NUGA::ISO_MIN); //evaluated on this element
      double LREF2_AMBIANT = NUGA::FLOAT_MAX;
      for (auto c : cands)
      {
        auto pg = src_mesh.element(c);
        double lr2 = pg.Lref2(src_mesh.crd, NUGA::ISO_MIN); //evaluated on this element
        LREF2_AMBIANT = std::min(LREF2_AMBIANT, lr2);
      }
      if (phlr2 > K* LREF2_AMBIANT)
      {
        filled = true;
        adap_incr.cell_adap_incr[PHi] = 1;
        continue;
      }
    }

    // criterion based on ISO MAX
    if ((_metric_policy == NUGA::MAX) || (_metric_policy == NUGA::MIN_OR_MAX))
    {
      double phlr2 = ph.Lref2(m1.crd, NUGA::ISO_MAX); //evaluated on this element
      double LREF2_AMBIANT = -1.;
      for (auto c : cands)
      {
        auto pg = src_mesh.element(c);
        double lr2 = pg.Lref2(src_mesh.crd, NUGA::ISO_MAX); //evaluated on this element
        LREF2_AMBIANT = std::max(LREF2_AMBIANT, lr2);
      }
      if (phlr2 > K* LREF2_AMBIANT) {
        filled = true;
        adap_incr.cell_adap_incr[PHi] = 1;
        continue;
      }
    }
    
    // criterion based on ISO MEAN
    if (_metric_policy == NUGA::MEAN)
    {
      double phlr2 = ph.Lref2(m1.crd, NUGA::ISO_MEAN); //evaluated on this element
      double LREF2_AMBIANT = -1.;
      for (auto c : cands)
      {
        auto pg = src_mesh.element(c);
        double lr2 = pg.Lref2(src_mesh.crd, NUGA::ISO_MEAN); //evaluated on this element
        LREF2_AMBIANT = std::max(LREF2_AMBIANT, lr2);
      }
      if (phlr2 > K* LREF2_AMBIANT) {
        filled = true;
        adap_incr.cell_adap_incr[PHi] = 1;
        continue;
      }
    }
  }

  return filled;
}

///
template <typename mesh_t>
bool xsensor2<mesh_t>::update()
{
  std::map<E_Int, std::vector<E_Int>> new_candidates;
  
  ph_mesh_t m1;//fixme
  m1.crd = parent_t::_hmesh._crd;
  m1.cnt = parent_t::_hmesh._ng;
  m1.get_nodal_metric2(NUGA::ISO_MAX);//fixme : just here to construct ae. ISO_MAX to avoid useless KdTree construc in build_nodal_metric2
  
  pg_smesh_t& src_mesh = parent_t::_data;
  auto loc2 = src_mesh.get_localizer();

  for (auto it : _candidates)
  {
    E_Int PHi = it.first;
    auto& cands = it.second;

    E_Int nbc = parent_t::_hmesh._PHtree.nb_children(PHi);
    const E_Int* children = parent_t::_hmesh._PHtree.children(PHi);
    for (E_Int c = 0; c < nbc; ++c)
    {
      E_Int cphi = children[c];

      auto ae1 = m1.aelement(cphi);
      //double Lref21 = ae1.Lref2(m1.crd, MTYPE);

      //medith::write<ngon_type>("ae", m1.crd, m1.cnt, cphi);

      cands.clear();
      loc2->get_candidates(ae1, ae1.m_crd, cands, 0, RTOL); //return as 0-based (fixme for volumic, was 1-based)
      if (cands.empty()) continue;

      bool is_x = NUGA::COLLIDE::get_colliding(ae1, src_mesh, cands, 0, RTOL, false/*all of them*/);
      if (!is_x) continue;
      if (cands.empty()) continue;

      new_candidates[cphi] = cands;
    }

  }

  _candidates = new_candidates;
  _done = _candidates.empty();

  return !_done;
}

///
template <typename mesh_t>
bool xsensor2<mesh_t>::stop()
{
  //std::cout << "stop ? : " << _iter << "/" << _iter_max << std::endl;
  if (++_iter > _iter_max) return true;
  return _done;
}

}


#endif


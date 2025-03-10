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

#include "Nuga/include/localizer.hxx"
#include "Nuga/include/collider.hxx"
#include <vector>

#ifndef NUGA_ESTIMATOR_HXX
#define NUGA_ESTIMATOR_HXX


namespace NUGA
{

  template <typename mesh_t1, typename mesh_t2, typename cell_incr_t>
  bool estimate_adap_req(mesh_t1& m1, mesh_t2& m2, NUGA::eMetricType M, double RTOL, std::vector<cell_incr_t>& data, E_Int MINVAL = 0, E_Int MAXVAL = 10)
  {
    data.clear();
    data.resize(m1.ncells(), cell_incr_t(0));

    m1.get_nodal_metric2(M); // ISO_MAX, ISO_MIN or ISO_MEAN
    m2.get_nodal_metric2(M);

    auto loc2 = m2.get_localizer();

    std::vector<E_Int> cands;
    //E_Int maxv = 0;
    //E_Int count{ 0 };

    /*std::cout << "estimate_adap_req : M : " << M <<std::endl;
    std::cout << "estimate_adap_req : RTOL : " << RTOL <<std::endl;
    std::cout << "estimate_adap_req : MINVAL : " << MINVAL <<std::endl;
    std::cout << "estimate_adap_req : MAXVAL : " << MAXVAL <<std::endl;*/

    for (E_Int i = 0; i < m1.ncells(); ++i)
    {
      auto ae1 = m1.aelement(i);
      double Lref21 = ae1.Lref2();

      cands.clear();
      loc2->get_candidates(ae1, ae1.m_crd, cands, 1, RTOL); //return as 0-based (fixme for volumic, was 1-based)
      if (cands.empty()) continue;

      bool is_x = NUGA::COLLIDE::get_colliding(ae1, m2, cands, 1, RTOL, false/*all of them*/);
      if (!is_x) continue;
      if (cands.empty()) continue;

      double cands_lref2 = 0.;
      if (M == NUGA::ISO_MIN) cands_lref2 = NUGA::FLOAT_MAX;

      for (size_t k = 0; k < cands.size(); ++k)
      {
        E_Int PGi = cands[k] - 1;
        auto as2 = m2.aelement(PGi);
        double Lref22 = as2.Lref2();

        if (M == NUGA::ISO_MEAN)
          cands_lref2 += Lref22;
        else if (M == NUGA::ISO_MAX)
          cands_lref2 = std::max(cands_lref2, Lref22);
        else if (M == NUGA::ISO_MIN)
          cands_lref2 = std::min(cands_lref2, Lref22);
      }

      if (M == NUGA::ISO_MEAN)
        cands_lref2 /= cands.size();

      /*bool caught = false;
      for (size_t uu = 0; uu < ae1.m_crd.cols(); ++uu)
      {
      if (ae1.m_crd(1, uu) != 0.) continue;
      if (ae1.m_crd(2, uu) != 0.) continue;

      if (ae1.m_crd(0, uu) + 0.150719 < 1.e-5)
      caught = true;
      }

      if (caught)
      {
      std::ostringstream o;
      o << "ae1_" << i;
      medith::write(o.str().c_str(), ae1);
      std::ostringstream o1;
      o1 << "s2pg_" << i;
      medith::write(o1.str().c_str(), s2.crd, s2.cnt, &cands, 1);
      }*/

      double r = 0.5*::log2(Lref21 / cands_lref2); // half because Lrefs are squares

      if (r < 1.) continue;
      E_Int nsub = E_Int(r);
      if (nsub < 0) continue;

      data[i] = nsub;

      //maxv = std::max(nsub, maxv);

      data[i] = NUGA::min(cell_incr_t(MAXVAL), data[i]);
      data[i] = NUGA::max(cell_incr_t(MINVAL), data[i]);

    }

    //std::cout << "maxval computed : " << maxv << std::endl;

    //std::cout << count << std::endl;
    return false;
  }

}

#endif

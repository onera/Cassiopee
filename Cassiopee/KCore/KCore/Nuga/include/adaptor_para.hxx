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

#ifndef NUGA_ADAPTOR_MPI_HXX
#define NUGA_ADAPTOR_MPI_HXX

#include "Nuga/include/mpi_stl.hxx"
#include "Nuga/include/adaptor.hxx"
#include "Nuga/include/join_sensor.hxx"
#ifdef _MPI
#include "Nuga/include/mpi_messages.hxx"
#include "Nuga/include/hybrid_para_algo.hxx"
#endif
#include "Nuga/include/omp_algo.hxx"
#ifdef ADAPT_TIMER
#include "Nuga/include/chrono.h"
#endif

namespace NUGA
{

  //
  template <typename para_algo_t, typename hmesh_t, typename sensor_t>
  class adaptor_para : public para_algo_t
  {
  public:

    using parent_t = para_algo_t;
    using zid_to_rid_to_ptlist_t = typename parent_t::zid_to_rid_to_ptlist_t;
    using id_to_PG_to_plan_t = typename parent_t::id_to_PG_to_plan_t;
    using rid_to_zones_t = typename parent_t::rid_to_zones_t;

    std::vector<sensor_t*> sensors;
    bool do_agglo;

    inline int get_data_stride() override { return 4; } //fixme : assume ISO


    ///
    bool prepare_data_to_send
    (
      const hmesh_t & mesh,
      const std::map<int, std::vector<E_Int>>& rid_to_list,
      id_to_PG_to_plan_t & rid_to_PG_to_plan
    ) override
    {
      bool has_packs{ false };
      //std::cout << "pack : loop on joins : " << join << std::endl;
      for (auto& it : rid_to_list)
      {
        int rid = it.first;
        auto& ptlist = it.second;

        for (size_t i = 0; i < ptlist.size(); ++i)
        {
          K_FLD::IntArray p;
          E_Int PGi = ptlist[i] - 1;
          //std::cout << "i/PGi/sz : " << i << "/" << PGi << "/" << ptlist.size() << std::endl;
          mesh.extract_plan(PGi, true/*reverse*/, 0/*because previous sync*/, p);
          //std::cout << "after extract_plan" << std::endl;
          if (p.getSize())
          {
            rid_to_PG_to_plan[rid][i/*PGi*/] = p;
            has_packs = true;
          }
        }
      }

      return has_packs;

    }

    void autonomous_run(const std::vector<hmesh_t*>& hmeshes, int i) override
    {
      using adaptor_t = NUGA::adaptor<hmesh_t, sensor_t>;
      adaptor_t::run(*hmeshes[i], *sensors[i], do_agglo);
    }

    bool run_with_data
    (const std::vector<hmesh_t*>& hmeshes,
     const std::vector<int>& zids,//useless for hmesh as it has a zid atribute
     const id_to_PG_to_plan_t & zid_to_PG_to_plan) override
    {
      E_Int NBZ{ E_Int(hmeshes.size()) };

      bool has_omp_changes{false};

      using join_sensor_t = NUGA::join_sensor<hmesh_t>;

      ePara PARA = COARSE_OMP;
 
#pragma omp parallel for reduction ( || : has_omp_changes) if(PARA == COARSE_OMP)         
      for (E_Int i = 0; i < NBZ; ++i)
      {
        join_sensor_t jsensor(*hmeshes[i]);
        auto it = zid_to_PG_to_plan.find(hmeshes[i]->zid);
        if (it == zid_to_PG_to_plan.end()) continue;
        
        jsensor.assign_data(it->second);

        E_Int npgs0 = hmeshes[i]->_ng.PGs.size();
        E_Int nphs0 = hmeshes[i]->_ng.PHs.size();

        NUGA::adaptor<hmesh_t, join_sensor_t>::run(*hmeshes[i], jsensor, do_agglo); //fixme : agglo

        has_omp_changes |= (hmeshes[i]->_ng.PGs.size() != npgs0) || (hmeshes[i]->_ng.PHs.size() != nphs0);
      }

      return has_omp_changes;
    }

  };

}

#endif

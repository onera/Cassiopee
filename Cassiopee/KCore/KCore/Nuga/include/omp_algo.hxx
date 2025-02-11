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

#include "Nuga/include/macros.h"
#include <map>
#include <vector>
#include <assert.h>

#ifndef NUGA_OMP_ALGO_HXX
#define NUGA_OMP_ALGO_HXX

namespace NUGA
{
  
  // zid and rid will probably never need to be long int
  // PG id or poition in pointlist must be E_Int

  enum ePara { SEQ, FINE_OMP, COARSE_OMP , DISTRIB };

  template <typename mesh_t, typename T>
  class omp_algo
  {
  public:

    using zid_to_rid_to_ptlist_t = std::map<int, std::map<int, std::vector<E_Int>>>;
    using id_to_PG_to_plan_t = std::map<int, std::map<E_Int, K_FLD::DynArray<T>>>;
    using rid_to_zones_t = std::map<int, std::pair<int, int>>;

    virtual bool prepare_data_to_send
    (
      const mesh_t & mesh,
      const std::map<int, std::vector<E_Int>>& rid_to_list,
      id_to_PG_to_plan_t & rid_to_PG_to_plan
    ) = 0;

    virtual void autonomous_run(const std::vector<mesh_t*>& mesh, int i) = 0;

    virtual bool run_with_data(const std::vector<mesh_t*>& meshes,
                               const std::vector<int>& zids,
                               const id_to_PG_to_plan_t & zid_to_PG_to_plan) = 0;

    virtual int get_data_stride() = 0;

    E_Int get_opp_zone(const rid_to_zones_t& rid_to_zones, E_Int rid, E_Int zid)
    {
      auto it = rid_to_zones.find(rid);
      if (it == rid_to_zones.end()){ assert(false); return IDX_NONE;}
      if (it->second.first == zid) return it->second.second;
      return it->second.first;
    }

    ///
    void run
    (
      std::vector<mesh_t*>& meshes,
      std::vector<int>& zids,
      const zid_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones
    );

    bool exchange_and_run
    (
      const std::vector<mesh_t*>& hmeshes,
      const std::vector<int>& zids,
      const zid_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones,
      id_to_PG_to_plan_t& zid_to_PG_to_plan
    );

    void exchange_omp_data
    (
      const std::vector<mesh_t*>& hmeshes,
      const std::vector<int>& zids,
      const zid_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones,
      id_to_PG_to_plan_t& zid_to_PG_to_plan
    );

  };

  ///
  template <typename mesh_t, typename T>
  void omp_algo<mesh_t, T>::run
  (
    std::vector<mesh_t*>& meshes,
    std::vector<int>& zids,
    const zid_to_rid_to_ptlist_t& zone_to_rid_to_list,
    const rid_to_zones_t& rid_to_zones
  )
  {
    int NBZ{ int(meshes.size()) };

    //1. autonomous runs
    NUGA::ePara PARA = COARSE_OMP;
#pragma omp parallel for if(PARA == COARSE_OMP)
    for (int i = 0; i < NBZ; ++i)
      this->autonomous_run(meshes, i);

    //2. exchange and run untill convergence
    id_to_PG_to_plan_t zid_to_PG_to_plan;
    exchange_and_run(meshes, zids, zone_to_rid_to_list, rid_to_zones, zid_to_PG_to_plan);

  }


  template <typename mesh_t, typename T>
  bool omp_algo<mesh_t, T>::exchange_and_run
  (
    const std::vector<mesh_t*>& meshes,
    const std::vector<int>& zids,
    const zid_to_rid_to_ptlist_t& zid_to_rid_to_list,
    const rid_to_zones_t& rid_to_zones,
    id_to_PG_to_plan_t& zid_to_PG_to_plan)
  {
    bool has_omp_changes{ true }, has_local_changes{ false };
    int omp_iter = -1;
    

    while (has_omp_changes)
    {
      ++omp_iter;
      //std::cout << "rank : " << rank << " : C : omp iter : " << omp_iter << std::endl;

      exchange_omp_data(meshes, zids, zid_to_rid_to_list, rid_to_zones, zid_to_PG_to_plan); //zid_to_PG_to_plan is appended with local contributions (might have distant contrib upon entry)

      has_omp_changes = this->run_with_data(meshes, zids, zid_to_PG_to_plan); // OVERLOADED
      has_local_changes |= has_omp_changes;
    }

    return has_local_changes;
  }


  ///
  template <typename mesh_t, typename T>
  void omp_algo<mesh_t, T>::exchange_omp_data
  (
    const std::vector<mesh_t*>& meshes,
    const std::vector<int>& zids,
    const zid_to_rid_to_ptlist_t& zone_to_rid_to_list,
    const rid_to_zones_t& rid_to_zones,
    id_to_PG_to_plan_t & zid_to_PG_to_plan
  )
  {
    bool has_packs{ false };

    for (size_t i = 0; i < meshes.size(); ++i)
    {
      int zid = zids[i];
      const auto it = zone_to_rid_to_list.find(zid);

      if (it == zone_to_rid_to_list.end()) continue; // current zone has no OMP joins

      const auto & rid_to_list = it->second;

      std::map<int, std::map<E_Int, K_FLD::DynArray<T>>> rid_to_PGrk_to_plan; //PGrk is the position in the pointlist
      has_packs |= this->prepare_data_to_send(*meshes[i], rid_to_list, rid_to_PGrk_to_plan);

      // convert to sensor data
      for (auto& r : rid_to_PGrk_to_plan)
      {
        int rid = r.first;
        int jzid = get_opp_zone(rid_to_zones, rid, zid);
        //get the OPP ptlist
        const auto itopp = zone_to_rid_to_list.find(jzid);
        assert(itopp != zone_to_rid_to_list.end());

        const auto itptl = itopp->second.find(rid);
        assert(itptl != itopp->second.end());
        const auto& ptlist = itptl->second;

        const auto & PG_to_plan = r.second;

        if (jzid != zid)
        {
          for (const auto & k : PG_to_plan)
          {
            const E_Int& j = k.first;
            const auto& plan = k.second;
            E_Int PGi = ptlist[j] - 1;
            zid_to_PG_to_plan[jzid][PGi] = plan; // we pass the plan associated to j-th face of zid to the correponding face in the joined zone
          }
        }
        else // auto-join
        {
          E_Int sz = ptlist.size();
          assert (sz % 2 == 0); // must always be even
          E_Int sz2 = sz / 2;

          for (const auto & k : PG_to_plan)
          {
            const E_Int& j = k.first;
            const auto& plan = k.second;

            E_Int j2 = (j + sz2) % sz; // j2 is the rank in the appropriate half of ptlist associated with j-th face in second half

            E_Int PGi = ptlist[j2] - 1;
            zid_to_PG_to_plan[zid][PGi] = plan;
          }
        }
      }
    }
  }
  
}

#endif

/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#include "Nuga/include/macros.h"
#include <map>
#include <vector>
#include <assert.h>

#ifndef NUGA_OMP_ALGO_HXX
#define NUGA_OMP_ALGO_HXX

namespace NUGA
{
  



  template <typename mesh_t, typename T>
  class omp_algo
  {
  public:

    using zone_to_rid_to_ptlist_t = std::map<int, std::map<int, std::vector<E_Int>>>;
    using rid_to_zones_t = std::map<int, std::pair<int, int>>;

    virtual bool prepare_data_to_send
    (
      const mesh_t & mesh,
      const std::map<int, std::vector<E_Int>>& rid_to_list,
      std::map<int, std::map<int, K_FLD::DynArray<T>>> & rid_to_PG_to_plan
    ) = 0;

    virtual void autonomous_run(const std::vector<mesh_t*>& mesh, int i) = 0;

    virtual bool run_with_data(const std::vector<mesh_t*>& meshes, const std::map<int, std::map<int, K_FLD::DynArray<T>>> & zid_to_data) = 0;

    virtual int get_data_stride() = 0;

    E_Int get_opp_zone(const std::map<E_Int, std::pair<E_Int,E_Int>>& rid_to_zones, E_Int rid, E_Int zid)
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
      const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_rid_to_list,
      const std::map<E_Int, std::pair<int, int>>& rid_to_zones
    );

    void exchange_and_run
    (
      const std::vector<mesh_t*>& hmeshes,
      const std::vector<int>& zids,
      const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones
    );

    void exchange_omp_data
    (
      const std::vector<mesh_t*>& hmeshes,
      const std::vector<int>& zids,
      const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones,
      std::map<int, std::map<int, K_FLD::DynArray<T>>> & zid_to_data
    );

  };

  ///
  template <typename mesh_t, typename T>
  void omp_algo<mesh_t, T>::run
  (
    std::vector<mesh_t*>& meshes,
    std::vector<int>& zids,
    const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_rid_to_list,
    const std::map<E_Int, std::pair<int, int>>& rid_to_zones
  )
  {
    E_Int err(0);
    E_Int NBZ{ E_Int(meshes.size()) };

    //1. autonomous runs
    ePara PARA = COARSE_OMP;
#pragma omp parallel for if(PARA == COARSE_OMP)
    for (E_Int i = 0; i < NBZ; ++i)
      this->autonomous_run(meshes, i);

    //2. exchange and run untill convergence
    exchange_and_run(meshes, zids, zone_to_rid_to_list, rid_to_zones);

  }


  template <typename mesh_t, typename T>
  void omp_algo<mesh_t, T>::exchange_and_run
  (
    const std::vector<mesh_t*>& meshes,
    const std::vector<int>& zids,
    const zone_to_rid_to_ptlist_t& zone_to_rid_to_list_omp,
    const rid_to_zones_t& rid_to_zones
  )
  {
    bool has_omp_changes{ true }, has_local_changes{ false };
    int omp_iter = -1;
    std::map<int, std::map<int, K_FLD::DynArray<T>>> zid_to_jdata;

    while (has_omp_changes)
    {
      ++omp_iter;
        //std::cout << "rank : " << rank << " : C : omp iter : " << omp_iter << std::endl;

      exchange_omp_data(meshes, zids, zone_to_rid_to_list_omp, rid_to_zones, zid_to_jdata);

      has_omp_changes = this->run_with_data(meshes, zid_to_jdata); // OVERLOADED
      has_local_changes |= has_omp_changes;
    }
  }


  ///
  template <typename mesh_t, typename T>
  void omp_algo<mesh_t, T>::exchange_omp_data
  (
    const std::vector<mesh_t*>& meshes,
    const std::vector<int>& zids,
    const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
    const rid_to_zones_t& rid_to_zones,
    std::map<int, std::map<int, K_FLD::DynArray<T>>> & zid_to_data
  )
  {
    //zid_to_data.clear();
    bool has_packs{ false };

    for (size_t i = 0; i < meshes.size(); ++i)
    {
      int zid = zids[i];
      const auto it = zone_to_rid_to_list.find(zid);

      if (it == zone_to_rid_to_list.end()) continue; // current zone has no OMP joins

      const auto & rid_to_list = it->second;

      std::map<int, std::map<int, K_FLD::DynArray<T>>> rid_to_PG_to_plan;
      has_packs |= this->prepare_data_to_send(*meshes[i], rid_to_list, rid_to_PG_to_plan);

      // convert to sensor data
      for (auto& r : rid_to_PG_to_plan)
      {
        int rid = r.first;
        int jzid = get_opp_zone(rid_to_zones, rid, zid);
        //get the OPP ptlist
        const auto itopp = zone_to_rid_to_list.find(jzid);
        assert(itopp != zone_to_rid_to_list.end());

        const auto itptl = itopp->second.find(rid);
        assert(itptl != itopp->second.end());
        const auto& ptlist = itptl->second;

        auto & PG_to_plan = r.second;

        if (jzid != zid)
        {
          for (auto & k : PG_to_plan)
          {
            E_Int PGi = ptlist[k.first] - 1;
            zid_to_data[jzid][PGi] = k.second;
          }
        }
        else // auto-join
        {
          int sz = ptlist.size() / 2; // outside of the loop to earn some calculation time
          int stock_Plan_size = PG_to_plan.size();

          for (auto & k : PG_to_plan) // copy refinement plan to apply it to adequate face
          {
            // keys we work with, each associated to a Plan
            int key1 = k.first; // 1st Plan (can be Left or Right, doesn't matter)
            int key2 = (k.first + sz) % ptlist.size(); // 2nd Plan // modulo as k.first may begin on List Right side (namely above ptlist.size()/2)

            E_Int PGi = ptlist[key1] - 1; // Face on which we test plan 1

            int PGj = ptlist[key2] - 1;
            zid_to_data[jzid][PGj] = k.second;       // zone to refine - Plan 1

            const auto Plan2 = PG_to_plan.find(key2);
            if (Plan2 == PG_to_plan.end()) continue;
            zid_to_data[jzid][PGi] = Plan2->second;  // zone to refine - Plan 2         
          }
        }
      }
    }
  }

  
}

#endif

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

    bool exchange_and_run
    (
      const std::vector<mesh_t*>& hmeshes,
      const std::vector<int>& zids,
      const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones,
      std::map<int, std::map<int, K_FLD::DynArray<T>>>& zid_to_jdata
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
    std::map<int, std::map<int, K_FLD::DynArray<T>>> zid_to_jdata;
    exchange_and_run(meshes, zids, zone_to_rid_to_list, rid_to_zones, zid_to_jdata);

  }


  template <typename mesh_t, typename T>
  bool omp_algo<mesh_t, T>::exchange_and_run
  (
    const std::vector<mesh_t*>& meshes,
    const std::vector<int>& zids,
    const zone_to_rid_to_ptlist_t& zone_to_rid_to_list_omp,
    const rid_to_zones_t& rid_to_zones,
    std::map<int, std::map<int, K_FLD::DynArray<T>>>& zid_to_jdata)
  {
    bool has_omp_changes{ true }, has_local_changes{ false };
    int omp_iter = -1;
    

    while (has_omp_changes)
    {
      ++omp_iter;
      //std::cout << "rank : " << rank << " : C : omp iter : " << omp_iter << std::endl;

      exchange_omp_data(meshes, zids, zone_to_rid_to_list_omp, rid_to_zones, zid_to_jdata); //zid_to_jdata is appended with local contributions (might have distant contrib upon entry)

      has_omp_changes = this->run_with_data(meshes, zid_to_jdata); // OVERLOADED
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
    const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
    const rid_to_zones_t& rid_to_zones,
    std::map<int, std::map<int, K_FLD::DynArray<T>>> & zid_to_data
  )
  {
    bool has_packs{ false };

    for (size_t i = 0; i < meshes.size(); ++i)
    {
      int zid = zids[i];
      const auto it = zone_to_rid_to_list.find(zid);

      if (it == zone_to_rid_to_list.end()) continue; // current zone has no OMP joins

      const auto & rid_to_list = it->second;

      std::map<int, std::map<E_Int, K_FLD::DynArray<T>>> rid_to_PG_to_plan;
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

        const auto & PG_to_plan = r.second;

        if (jzid != zid)
        {
          for (const auto & k : PG_to_plan)
          {
            const E_Int& j = k.first;
            const auto& plan = k.second;
            E_Int PGi = ptlist[j] - 1;
            zid_to_data[jzid][PGi] = plan; // we pass the plan associated to j-th face of zid to the correponding face in the joined zone
          }
        }
        else // auto-join
        {
          E_Int sz = ptlist.size();
          E_Int sz2 = sz / 2;

          for (const auto & k : PG_to_plan)
          {
            const E_Int& j = k.first;
            const auto& plan = k.second;

            E_Int j2 = (j + sz2) % sz; // j2 is the rank in the appropriate half of ptlist associated with j-th face in second half

            E_Int PGi = ptlist[j2] - 1;
            zid_to_data[zid][PGi] = k.second;
          }
        }
      }
    }
  }
  
}

#endif

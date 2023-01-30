/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_ADAPTOR_MPI_HXX
#define NUGA_ADAPTOR_MPI_HXX

#include "Nuga/include/mpi_stl.hxx"
#include "Nuga/include/mpi_messages.hxx"
#include "Nuga/include/adaptor.hxx"
#include "Nuga/include/join_sensor.hxx"
#include "Nuga/include/hybrid_para_algo.hxx"
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

    std::vector<sensor_t*> sensors;
    bool do_agglo;

    inline int get_data_stride() override { return 4; } //fixme : assume ISO


    ///
    bool prepare_data_to_send
    (
      const hmesh_t & mesh,
      const std::map<int, std::vector<E_Int>>& rid_to_list,
      std::map<int, std::map<E_Int, K_FLD::IntArray>> & rid_to_PG_to_plan
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
    (const std::vector<hmesh_t*>& hmeshes, const std::map<int, std::map<E_Int, K_FLD::IntArray>> & zid_to_PG_to_plan) override
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

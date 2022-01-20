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


namespace NUGA
{
  

  template <typename hmesh_t, typename sensor_t>
  class adaptor_mpi
  {
  public:

    ///
    template <typename data_t>
    static bool prepare_data_to_send
    (
      const hmesh_t & mesh,
      const std::map<E_Int, std::vector<E_Int>>& zone2jlist,
      std::map<E_Int, data_t> & raczone_to_PG_to_plan
    )
    {
      bool has_packs{ false };
      //std::cout << "pack : loop on joins : " << join << std::endl;
      for (auto& it : zone2jlist)
      {
        E_Int jzid = it.first;
        auto& ptlist = it.second;

        for (size_t i = 0; i < ptlist.size(); ++i)
        {
          typename hmesh_t::pg_arr_t p; //only implemented for IntArray
          E_Int PGi = ptlist[i] - 1;
          //std::cout << "i/PGi/sz : " << i << "/" << PGi << "/" << ptlist.size() << std::endl;
          mesh.extract_plan(PGi, true/*reverse*/, 0/*because previous sync*/, p);
          //std::cout << "after extract_plan" << std::endl;
          if (p.getSize())
          {
            raczone_to_PG_to_plan[jzid][i/*PGi*/] = p;
            has_packs = true;
          }
        }
      }

      return has_packs;

    }

    //
    template <typename data_t>
    static void exchange_mpi_data
    (
      const std::vector<hmesh_t*>& hmeshes,
      const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone2jlists,
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      int rank, int nranks,
      std::map<int, data_t>& zone_to_sensor_data
    )
    {
      if (zone_to_zone2jlists.empty()) return;
      if (hmeshes.empty()) return;

      // Each zone builds its MPI-data-to-send by zone
      std::map<int, std::map<int, std::map<int, K_FLD::IntArray>>> sz_to_jz_to_PG_to_plan;
      bool has_packs{ false };
      for (size_t i = 0; i < hmeshes.size(); ++i)
      {
        //if (rank == 2) std::cout << "rank : " << rank << " : hmeshes[" << i << "] : " << hmeshes[i] << std::endl;
        const auto it = zone_to_zone2jlists.find(hmeshes[i]->zid);
        
        if (it == zone_to_zone2jlists.end()) continue; // current zone has no MPI joins

        //if (rank == 2) std::cout << "rank : " << rank << " found zid : " << hmeshes[i]->zid << std::endl;

        const auto & zone2jlists = it->second;
        
        std::map<int, std::map<int, K_FLD::IntArray>> jz_to_PG_to_plan;
        has_packs |= prepare_data_to_send(*hmeshes[i], zone2jlists, jz_to_PG_to_plan);

        sz_to_jz_to_PG_to_plan[hmeshes[i]->zid] = jz_to_PG_to_plan;

      }

      // Convert for sending and gather by rank : rank/zone/face/plan
      //if (rank == 2) std::cout << "rank : " << rank << " convert_to_MPI_exchange_format ..." << std::endl;
      std::map<int, plan_msg_type> rank_to_mpi_data;
      plan_msg_type::convert_to_MPI_exchange_format(sz_to_jz_to_PG_to_plan, zonerank, rank_to_mpi_data);

      // Send MPI data   
      //if (rank == 2) std::cout << "rank : " << rank << " send_data ..." << std::endl;
      plan_msg_type::send_data(rank, nranks, COM, rank_to_mpi_data);

      // Receive MPI data and build sensor data by zone
      //if (rank == 2) std::cout << "rank : " << rank << " receive_data ..." << std::endl;
      plan_msg_type::receive_data(rank, nranks, COM, zone_to_zone2jlists, zone_to_sensor_data);
      //if (rank == 2) std::cout << "rank : " << rank << " DONE. exit exch mpi" << std::endl;
    }

    //
    template <typename data_t>
    static void exchange_omp_data
    (
      const std::vector<hmesh_t*>& hmeshes,
      const std::map<int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone2jlists,
      std::map<int, data_t>& sensor_data
    )
    {
      
      bool has_packs{ false };
      for (size_t i = 0; i < hmeshes.size(); ++i)
      {
        const auto it = zone_to_zone2jlists.find(hmeshes[i]->zid);

        if (it == zone_to_zone2jlists.end()) continue; // current zone has no OMP joins

        const auto & zone2jlists = it->second;

        std::map<int, data_t> jz_to_PG_to_plan;
        has_packs |= prepare_data_to_send(*hmeshes[i], zone2jlists, jz_to_PG_to_plan);

        // convert to sensor data
        for (auto& j : jz_to_PG_to_plan)
        {
          int jzid = j.first;

          //get the ptlist
          const auto it2 = zone_to_zone2jlists.find(jzid);
          if (it2 == zone_to_zone2jlists.end()) continue;
          const auto & zone2jlists2 = it2->second;
          const auto& itPtList = zone2jlists2.find(hmeshes[i]->zid);
          assert(itPtList != zone2jlists2.end());
          const auto& ptlist = itPtList->second;
          
          auto & PG_to_plan = j.second;

          for (auto & k : PG_to_plan)
          {
            E_Int PGi = ptlist[k.first] - 1;
            sensor_data[jzid][PGi] = k.second;
          }
        }
      }
    }

    ///
    static void run
    (
      std::vector<hmesh_t*>& hmeshes, std::vector<sensor_t*>& sensors,
      const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone2jlists,
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      bool do_agglo
    )
    {
      E_Int err(0);
      E_Int NBZ{ E_Int(hmeshes.size()) };

      assert(NBZ == sensors.size());

      using adaptor_t = NUGA::adaptor<hmesh_t, sensor_t>;
      using join_sensor_t = NUGA::join_sensor<hmesh_t>;
      using data_t = typename join_sensor_t::input_t;// std::map<int, K_FLD::IntArray>;

      ePara PARA = COARSE_OMP;
#pragma omp parallel for if(PARA == COARSE_OMP)
      for (E_Int i = 0; i < NBZ; ++i)
      {
        adaptor_t::run(*hmeshes[i], *sensors[i], do_agglo);
      }

      int nranks{ 0 };
      MPI_Comm_size(COM, &nranks);
      int rank{ 0 };
      MPI_Comm_rank(COM, &rank);

      for (E_Int i = 0; i < NBZ; ++i)
        assert(rank == zonerank[hmeshes[i]->zid]);

      //separate MPI/OMP joins
      std::map<int, std::map<E_Int, std::vector<E_Int>>> zone_to_zone2jlists_mpi, zone_to_zone2jlists_omp;
      NUGA::split_mpi_omp_joins(zone_to_zone2jlists, rank, zonerank, zone_to_zone2jlists_omp, zone_to_zone2jlists_mpi);

      bool has_mpi_changes{ true };

      int mpi_iter = -1;

      while (has_mpi_changes)
      {
        ++mpi_iter;
        //std::cout << "rank : " << rank << " mpi iter : " << mpi_iter << std::endl;
        //if (mpi_iter == 1) return;
        has_mpi_changes = false;

        std::map<int, data_t> zone_to_sensor_data;
        ///*if (rank ==2)*/ std::cout << "rank : " << rank << " before exch mpi" << std::endl;
        exchange_mpi_data(hmeshes, zone_to_zone2jlists_mpi, zonerank, COM, rank, nranks, zone_to_sensor_data);
        
        bool has_omp_changes{ true }, has_local_changes{ false };
        int omp_iter = -1;

        while (has_omp_changes)
        {
          ++omp_iter;
          //std::cout << "rank : " << rank << " omp iter : " << mpi_iter << std::endl;

          has_omp_changes = false;
          
          exchange_omp_data(hmeshes, zone_to_zone2jlists_omp, zone_to_sensor_data);

          // Adapt each zone
          //std::cout << "rank : " << rank << " : adapt omp loop" << mpi_iter << std::endl;
#pragma omp parallel for reduction ( || : has_omp_changes, has_local_changes) if(PARA == COARSE_OMP)
          for (E_Int i = 0; i < NBZ; ++i)
          {
            join_sensor_t jsensor(*hmeshes[i]);
            jsensor.assign_data(zone_to_sensor_data[hmeshes[i]->zid]);

            E_Int npgs0 = hmeshes[i]->_ng.PGs.size();
            E_Int nphs0 = hmeshes[i]->_ng.PHs.size();

            NUGA::adaptor<hmesh_t, join_sensor_t>::run(*hmeshes[i], jsensor, false/*do_agglo*/); //fixme : agglo

            has_omp_changes |= (hmeshes[i]->_ng.PGs.size() != npgs0) || (hmeshes[i]->_ng.PHs.size() != nphs0);
            has_local_changes |= has_omp_changes;
          }
        }
        MPI_Barrier(COM);
        //if (rank == 0) std::cout << "rank : " << rank << " MPI_Barrier " << std::endl;
        has_mpi_changes = exchange_status(has_local_changes, rank, nranks, COM);
      }
    }



    ///
    static bool exchange_status (bool local_changes, int rank, int nranks, MPI_Comm COM)
    {
      bool has_changes{ false };

      STACK_ARRAY(bool, nranks, all_status);
      MPI_Request req;
      if (rank != 0) {
        MPI_Isend(&local_changes, 1, MPI_C_BOOL, 0, TAG_MPI_EXCH_STATUS, COM, &req);
        //if (rank == 1) std::cout << "rank : " << rank << " is sendind its change status : " << has_changes << std::endl;
      }
      else
      {
        all_status[0] = local_changes;
        MPI_Status status;
        for (size_t r = 1; r < nranks; ++r) {
          MPI_Recv(&all_status[r], 1, MPI_C_BOOL, r, TAG_MPI_EXCH_STATUS, COM, &status);
          //std::cout << "rank 0 has received status from rank " << r << std::endl;
        }
      }

      //if (rank == 1) std::cout << "rank : " << rank << " is before barrier" << std::endl;

      MPI_Barrier(COM);

      //std::cout << "after barrier" << std::endl;

      if (rank == 0)
      {
        has_changes = false;
        for (size_t k = 0; (k < nranks) && !has_changes; ++k)
          has_changes |= all_status[k];

        //std::cout << "rank : 0 check the change status among received: " << has_changes << std::endl;

        STACK_ARRAY(MPI_Request, nranks - 1, req);
        STACK_ARRAY(MPI_Status, nranks - 1, status);

        for (size_t r = 1; r < nranks; ++r)
          MPI_Isend(&has_changes, 1, MPI_C_BOOL, r, TAG_MPI_EXCH_STATUS, COM, &req[r - 1]);
        MPI_Waitall(nranks - 1, req.get(), status.get());
        //std::cout << "rank : 0 is telling to all change status : " << has_changes << std::endl;
      }
      else
      {
        MPI_Status status;
        MPI_Recv(&has_changes, 1, MPI_C_BOOL, 0, TAG_MPI_EXCH_STATUS, COM, &status);
        //std::cout << "rank : " << rank << "has been told change status : " << has_changes << std::endl;
      }

      return has_changes;
    }

  };

}

#endif

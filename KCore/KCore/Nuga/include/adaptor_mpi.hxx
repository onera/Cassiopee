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

#ifdef ADAPT_TIMER
#include "Nuga/include/chrono.h"
#endif

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
      const std::map<E_Int, std::vector<E_Int>>& rid_to_list,
      std::map<E_Int, data_t> & rid_to_PG_to_plan
    )
    {
      bool has_packs{ false };
      //std::cout << "pack : loop on joins : " << join << std::endl;
      for (auto& it : rid_to_list)
      {
        E_Int rid = it.first;
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
            rid_to_PG_to_plan[rid][i/*PGi*/] = p;
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
      const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_rid_to_list,
      const std::map<int, std::pair<int, int>>& rid_to_zones,
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      int rank, int nranks,
      std::map<int, data_t>& zone_to_sensor_data
    )
    {
      if (zone_to_rid_to_list.empty()) return;
      if (hmeshes.empty()) return;
      if (nranks == 1) return;

      // Each zone builds its MPI-data-to-send by sending zone
      std::map<int, std::map<int, std::map<int, K_FLD::IntArray>>> sz_to_rid_to_PG_to_plan;
      bool has_packs{ false };
      for (size_t i = 0; i < hmeshes.size(); ++i)
      {
        //if (rank == 2) std::cout << "rank : " << rank << " : hmeshes[" << i << "] : " << hmeshes[i] << std::endl;
        const auto it = zone_to_rid_to_list.find(hmeshes[i]->zid);
        
        if (it == zone_to_rid_to_list.end()) continue; // current zone has no MPI joins

        //if (rank == 2) std::cout << "rank : " << rank << " found zid : " << hmeshes[i]->zid << std::endl;

        const auto & rid_to_list = it->second;
        
        std::map<int, std::map<int, K_FLD::IntArray>> rid_to_PG_to_plan;
        has_packs |= prepare_data_to_send(*hmeshes[i], rid_to_list, rid_to_PG_to_plan);

        sz_to_rid_to_PG_to_plan[hmeshes[i]->zid] = rid_to_PG_to_plan;

      }

      // Convert for sending and gather by rank : rank/zone/face/plan
      //if (rank == 2) std::cout << "rank : " << rank << " convert_to_MPI_exchange_format ..." << std::endl;
      std::map<int, plan_msg_type> rank_to_mpi_data;
      plan_msg_type::convert_to_MPI_exchange_format(sz_to_rid_to_PG_to_plan, rid_to_zones, zonerank, rank_to_mpi_data);

      // Send MPI data   
      //if (rank == 2) std::cout << "rank : " << rank << " send_data ..." << std::endl;

      int nsranks = rank_to_mpi_data.size(); // nb of ranks to send to
      int NB_TOT_REQS = nranks + (7 * nsranks);    // 'has_sent' to all ranks (nranks)= + 7 vectors for each connected rank : data/datarange/pgs/szone/szonerange/jzone/jzonerange.  2 req per vec. => 14

      std::vector<MPI_Request> sreqs(NB_TOT_REQS);
      STACK_ARRAY(bool, nranks, has_sent);
      for (size_t n = 0; n < nranks; ++n) has_sent[n] = false;

      plan_msg_type::isend_data(rank, nranks, COM, rank_to_mpi_data, has_sent.get(), sreqs);

      // Receive MPI data and build sensor data by zone
      //if (rank == 2) std::cout << "rank : " << rank << " receive_data ..." << std::endl;
      plan_msg_type::receive_data(rank, nranks, COM, rid_to_zones, zone_to_rid_to_list, sreqs, zone_to_sensor_data);
      //if (rank == 2) std::cout << "rank : " << rank << " DONE. exit exch mpi" << std::endl;
    }

    //
    template <typename data_t>
    static void exchange_omp_data
    (
      const std::vector<hmesh_t*>& hmeshes,
      const std::map<int, std::map<E_Int, std::vector<E_Int>>>& zone_to_rid_to_list,
      const std::map<int, std::pair<int,int>> & rid_to_zones,
      std::map<int, data_t>& sensor_data
    )
    {
      
      bool has_packs{ false };
      for (size_t i = 0; i < hmeshes.size(); ++i)
      {
        int zid = hmeshes[i]->zid;
        const auto it = zone_to_rid_to_list.find(zid);

        if (it == zone_to_rid_to_list.end()) continue; // current zone has no OMP joins

        const auto & rid_to_list = it->second;

        std::map<int, data_t> rid_to_PG_to_plan;
        has_packs |= prepare_data_to_send(*hmeshes[i], rid_to_list, rid_to_PG_to_plan);

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

          if (jzid != hmeshes[i]->zid)
          {
            for (auto & k : PG_to_plan)
            {
              E_Int PGi = ptlist[k.first] -1;
              sensor_data[jzid][PGi] = k.second;
            }
          }
          
          else // auto-join
          {
            int sz = ptlist.size() / 2; // outside of the loop to earn some calculation time
            int stock_Plan_size = PG_to_plan.size();
            
            for (auto & k : PG_to_plan) // copy refinement plan to apply it to adequate face
            {
              // keys we work with, each associated to a Plan
              int key1 =  k.first; // 1st Plan (can be Left or Right, doesn't matter)
              int key2 = (k.first + sz) % ptlist.size(); // 2nd Plan // modulo as k.first may begin on List Right side (namely above ptlist.size()/2)
              
              // Plan 1
              E_Int PGi = ptlist[key1] -1; // Face on which we test plan 1
              
              // ------------------------------------- //
              
              // Plan 2
              int PGj = ptlist[key2] -1;
              sensor_data[jzid][PGj] = k.second;       // zone to refine - Plan 1

              const auto Plan2 = PG_to_plan.find(key2);
              if (Plan2 == PG_to_plan.end()) continue;
              sensor_data[jzid][PGi] = Plan2->second;  // zone to refine - Plan 2         
            }
          }
        }
      }
    }

    ///
    static void run
    (
      std::vector<hmesh_t*>& hmeshes, std::vector<sensor_t*>& sensors,
      const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_rid_to_list,
      const std::map<E_Int, std::pair<int, int>>& rid_to_zones,
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      bool do_agglo
    )
    {
      E_Int err(0);
      E_Int NBZ{ E_Int(hmeshes.size()) };

      int nranks{ 0 };
      MPI_Comm_size(COM, &nranks);
      int rank{ 0 };
      MPI_Comm_rank(COM, &rank);

      for (E_Int i = 0; i < NBZ; ++i)
        assert(rank == zonerank[hmeshes[i]->zid]);

      assert(NBZ == sensors.size());

      using adaptor_t = NUGA::adaptor<hmesh_t, sensor_t>;
      using join_sensor_t = NUGA::join_sensor<hmesh_t>;
      using data_t = typename join_sensor_t::input_t;// std::map<int, K_FLD::IntArray>;

#ifdef ADAPT_TIMER
      NUGA::chrono c;
      c.start();
#endif

      ePara PARA = COARSE_OMP;
#pragma omp parallel for if(PARA == COARSE_OMP)
      for (E_Int i = 0; i < NBZ; ++i)
      {
        adaptor_t::run(*hmeshes[i], *sensors[i], do_agglo);
      }

#ifdef ADAPT_TIMER
      //if (rank == 0)
      //std::cout << "rank : " << rank << " : C : time independant adapt : " << c.elapsed()<< std::endl;
#endif

      //separate MPI/OMP joins
      std::map<int, std::map<E_Int, std::vector<E_Int>>> zone_to_rid_to_list_mpi, zone_to_rid_to_list_omp;
      NUGA::split_mpi_omp_joins(zone_to_rid_to_list, rank, rid_to_zones, zonerank, zone_to_rid_to_list_omp, zone_to_rid_to_list_mpi);

      bool has_mpi_changes{ true };

      int mpi_iter = -1;

      while (has_mpi_changes)
      {
        ++mpi_iter;
        //std::cout << "rank : " << rank << " mpi iter : " << mpi_iter << std::endl;
        //if (mpi_iter == 1) return;
        has_mpi_changes = false;

        std::map<int, data_t> zone_to_sensor_data;

#ifdef ADAPT_TIMER
        //c.start();
#endif
        exchange_mpi_data(hmeshes, zone_to_rid_to_list_mpi, rid_to_zones, zonerank, COM, rank, nranks, zone_to_sensor_data);
        
#ifdef ADAPT_TIMER
        //if (rank == 0)
        //std::cout << "rank : " << rank << " : C : time exchange_mpi_data : " << c.elapsed() << std::endl;
#endif

        bool has_omp_changes{ true }, has_local_changes{ false };
        int omp_iter = -1;

#ifdef ADAPT_TIMER
        /*c.start();
        NUGA::chrono c2;
        double dexch=0.;
        double dadapt=0.;
        double dsensor=0.;*/
#endif

        while (has_omp_changes)
        {
          ++omp_iter;
          //std::cout << "rank : " << rank << " : C : omp iter : " << omp_iter << std::endl;

          has_omp_changes = false;

#ifdef ADAPT_TIMER
          //c2.start();
#endif

          exchange_omp_data(hmeshes, zone_to_rid_to_list_omp, rid_to_zones, zone_to_sensor_data);
          
#ifdef ADAPT_TIMER
          //dexch += c2.elapsed();
          // Adapt each zone
          //std::cout << "rank : " << rank << " : adapt omp loop" << mpi_iter << std::endl;
#endif


//no reduction mode  STACK_ARRAY(bool, NBZ, has_local_changes);
//no reduction mode  for (size_t kkk = 0; kkk < NBZ; ++kkk)has_local_changes[kkk]=false;

#pragma omp parallel for reduction ( || : has_omp_changes, has_local_changes) if(PARA == COARSE_OMP)
//no reduction mode #pragma omp parallel for if(PARA == COARSE_OMP)          
          for (E_Int i = 0; i < NBZ; ++i)
          {

#ifdef ADAPT_TIMER
            //c2.start();
#endif
            join_sensor_t jsensor(*hmeshes[i]);
            jsensor.assign_data(zone_to_sensor_data[hmeshes[i]->zid]);
            
#ifdef ADAPT_TIMER
            //dsensor += c2.elapsed();
#endif

            E_Int npgs0 = hmeshes[i]->_ng.PGs.size();
            E_Int nphs0 = hmeshes[i]->_ng.PHs.size();

#ifdef ADAPT_TIMER
            //c2.start();
#endif

            NUGA::adaptor<hmesh_t, join_sensor_t>::run(*hmeshes[i], jsensor, false/*do_agglo*/); //fixme : agglo

#ifdef ADAPT_TIMER
            //dadapt += c2.elapsed();
#endif

            has_omp_changes |= (hmeshes[i]->_ng.PGs.size() != npgs0) || (hmeshes[i]->_ng.PHs.size() != nphs0);
            has_local_changes |= has_omp_changes;

//no reduction mode            has_local_changes[i] = (hmeshes[i]->_ng.PGs.size() != npgs0) || (hmeshes[i]->_ng.PHs.size() != nphs0);
          }

 //no reduction mode         for (size_t kkk=0; kkk < NBZ; ++kkk)
 //no reduction mode           has_omp_changes |= has_local_changes[kkk];

          //zone_to_sensor_data.clear();

        }

#ifdef ADAPT_TIMER
        //if (rank == 0)
        /*{
          //std::cout << "rank : " << rank << " : C : OMP LOOP : time creation data sensor : " << dexch << std::endl;
          std::cout << "rank : " << rank << " : C : OMP LOOP : time adapt sync : " << dadapt << std::endl;

          //std::cout << "rank : " << rank << " : C : OMP LOOP : time TOTAL : " << c.elapsed() << std::endl;
        }*/
#endif
        
        MPI_Barrier(COM);

#ifdef ADAPT_TIMER
        //c.start();
#endif

        //if (rank == 0) std::cout << "rank : " << rank << " MPI_Barrier " << std::endl;
        has_mpi_changes = exchange_status(has_local_changes, rank, nranks, COM);

#ifdef ADAPT_TIMER
        //if (rank == 0)
        //std::cout << "rank : " << rank << " : C : OMP LOOP : time exchange_status : " << c.elapsed() << std::endl;
#endif
      }
    }



    ///
    static bool exchange_status (bool local_changes, int rank, int nranks, MPI_Comm COM)
    {
      bool has_changes{ false };
      if (nranks==1) return false;

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
          //std::cout << " has received status from rank " << r << std::endl;
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

        //std::cout << " check the change status among received: " << has_changes << std::endl;

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

/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#include "Nuga/include/omp_algo.hxx"
#include "Nuga/include/mpi_stl.hxx"
#include "Nuga/include/mpi_messages.hxx"
#include "Nuga/include/macros.h"
#include <map>
#include <vector>
#include <assert.h>

#ifndef NUGA_HYBRID_PARA_ALGO_HXX
#define NUGA_HYBRID_PARA_ALGO_HXX

namespace NUGA
{
  
  template <typename mesh_t, typename T>
  class hybrid_para_algo : public omp_algo<mesh_t, T>
  {
  public:

    using parent_t = omp_algo<mesh_t, T>;

    using zone_to_rid_to_ptlist_t = std::map<int, std::map<int, std::vector<E_Int>>>;
    using rid_to_zones_t = std::map<int, std::pair<int, int>>;

    ///
    void run
    (
      std::vector<mesh_t*>& meshes,
      std::vector<int>& zids,
      const std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list,
      const std::map<int, std::pair<int, int>>& rid_to_zones,
      const std::vector<int>& zonerank,
      MPI_Comm COM
    );

    void exchange_and_run
    (
      const std::vector<mesh_t*>& hmeshes,
      const std::vector<int>& zids,
      const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones,
      const std::vector<int>& zonerank,
      MPI_Comm COM
    );

    void exchange_mpi_data
    (
      const std::vector<mesh_t*>& meshes,
      const std::vector<int>& zids,
      const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
      const rid_to_zones_t& rid_to_zones,
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      int rank,
      int nranks,
      std::map<int, std::map<E_Int, K_FLD::DynArray<T>>> & zid_to_data
    );

    ///
    bool exchange_status(bool local_changes, int rank, int nranks, MPI_Comm COM);
    
  };

  ///
  template <typename mesh_t, typename T>
  void hybrid_para_algo<mesh_t, T>::run
  (
    std::vector<mesh_t*>& meshes,
    std::vector<int>& zids,
    const std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list,
    const std::map<int, std::pair<int, int>>& rid_to_zones,
    const std::vector<int>& zonerank,
    MPI_Comm COM
  )
  {
    int err(0);
    int NBZ{ E_Int(meshes.size()) };

    //1. autonomous runs
    ePara PARA = COARSE_OMP;
#pragma omp parallel for if(PARA == COARSE_OMP)
    for (int i = 0; i < NBZ; ++i)
      this->autonomous_run(meshes, i);

    //2. exchange and run untill convergence
    exchange_and_run(meshes, zids, zone_to_rid_to_list, rid_to_zones, zonerank, COM);

  }


  template <typename mesh_t, typename T>
  void hybrid_para_algo<mesh_t, T>::exchange_and_run
  (
    const std::vector<mesh_t*>& meshes,
    const std::vector<int>& zids,
    const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
    const rid_to_zones_t& rid_to_zones,
    const std::vector<int>& zonerank,
    MPI_Comm COM
  )
  {
    int err(0);

    int nranks{ 0 };
    MPI_Comm_size(COM, &nranks);

    int rank{ 0 };
    MPI_Comm_rank(COM, &rank);

    //separate MPI/OMP joins
    zone_to_rid_to_ptlist_t zone_to_rid_to_list_mpi, zone_to_rid_to_list_omp;
    split_mpi_omp_joins(zone_to_rid_to_list, rank, rid_to_zones, zonerank, zone_to_rid_to_list_omp, zone_to_rid_to_list_mpi);

    bool has_mpi_changes{ true };

    int mpi_iter = -1;

    while (has_mpi_changes)
    {
      ++mpi_iter;
      has_mpi_changes = false;

      std::map<int, std::map<E_Int, K_FLD::DynArray<T>>> zid_to_jdata;

      exchange_mpi_data(meshes, zids, zone_to_rid_to_list_mpi, rid_to_zones, zonerank, COM, rank, nranks, zid_to_jdata);

      bool has_local_changes = parent_t::exchange_and_run(meshes, zids, zone_to_rid_to_list_omp, rid_to_zones, zid_to_jdata);

      MPI_Barrier(COM);

      //if (rank == 0) std::cout << "rank : " << rank << " MPI_Barrier " << std::endl;
      has_mpi_changes = exchange_status(has_local_changes, rank, nranks, COM);
    }
  }

  ///
  template <typename mesh_t, typename T>
  void hybrid_para_algo<mesh_t, T>::exchange_mpi_data
  (
    const std::vector<mesh_t*>& meshes,
    const std::vector<int>& zids,
    const zone_to_rid_to_ptlist_t& zone_to_rid_to_list,
    const rid_to_zones_t& rid_to_zones,
    const std::vector<int>& zonerank,
    MPI_Comm COM,
    int rank,
    int nranks,
    std::map<int, std::map<E_Int, K_FLD::DynArray<T>>> & zid_to_data
  )
  {
    if (zone_to_rid_to_list.empty()) return;
    if (nranks == 1) return;

    // Each zone builds its MPI-data-to-send by sending zone
    std::map<int, std::map<int, std::map<E_Int, K_FLD::DynArray<T>>>> sz_to_rid_to_PG_to_data;
    bool has_packs{ false };
    for (size_t i = 0; i < meshes.size(); ++i)
    {
      //i-th mesh has zid[i] as id

      //if (rank == 2) std::cout << "rank : " << rank << " : hmeshes[" << i << "] : " << hmeshes[i] << std::endl;
      const auto it = zone_to_rid_to_list.find(zids[i]);

      if (it == zone_to_rid_to_list.end()) continue; // current zone has no MPI joins

      //if (rank == 2) std::cout << "rank : " << rank << " found zid : " << hmeshes[i]->zid << std::endl;
      const auto & rid_to_list = it->second;

      std::map<int, std::map<E_Int, K_FLD::DynArray<T>>> rid_to_PG_to_data;
      has_packs |= this->prepare_data_to_send(*meshes[i], rid_to_list, rid_to_PG_to_data);

      sz_to_rid_to_PG_to_data[zids[i]] = rid_to_PG_to_data;

    }

    // Convert for sending and gather by rank : rank/zone/face/data
    //if (rank == 2) std::cout << "rank : " << rank << " convert_to_MPI_exchange_format ..." << std::endl;
    using plan_type = NUGA::plan_msg_type<T>;
    std::map<int, plan_type> rank_to_mpi_data;
    plan_type::convert_to_MPI_exchange_format(sz_to_rid_to_PG_to_data, rid_to_zones, zonerank, rank_to_mpi_data);

    // Send MPI data   
    //if (rank == 2) std::cout << "rank : " << rank << " send_data ..." << std::endl;

    int nsranks = rank_to_mpi_data.size(); // nb of ranks to send to
    int NB_TOT_REQS = nranks + (7 * nsranks);    // 'has_sent' to all ranks (nranks)= + 7 vectors for each connected rank : data/datarange/pgs/szone/szonerange/jzone/jzonerange.  2 req per vec. => 14

    std::vector<MPI_Request> sreqs(NB_TOT_REQS);
    STACK_ARRAY(bool, nranks, has_sent);
    for (size_t n = 0; n < nranks; ++n) has_sent[n] = false;

    plan_type::isend_data(rank, nranks, COM, rank_to_mpi_data, has_sent.get(), sreqs);

    // Receive MPI data and build sensor data by zone
    //if (rank == 2) std::cout << "rank : " << rank << " receive_data ..." << std::endl;
    plan_type::receive_data(rank, nranks, COM, this->get_data_stride(), rid_to_zones, zone_to_rid_to_list, sreqs, zid_to_data);
    //if (rank == 2) std::cout << "rank : " << rank << " DONE. exit exch mpi" << std::endl;
  }

  ///
  template <typename mesh_t, typename T>
  bool hybrid_para_algo<mesh_t, T>::exchange_status
  (bool local_changes, int rank, int nranks, MPI_Comm COM)
  {
    bool has_changes{ false };
    if (nranks == 1) return false;

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
}

#endif

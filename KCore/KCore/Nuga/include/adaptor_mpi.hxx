/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#ifndef NUGA_ADAPTOR_MPI_HXX
#define NUGA_ADAPTOR_MPI_HXX

#include "mpi.h"
#include "NUGA/include/adaptor.hxx"
#include "Nuga/include/join_sensor.hxx"


namespace NUGA
{
  template <typename mesh_t, typename sensor_t>
  class adaptor_mpi
  {
  public:

    enum eMPI_Tag { TAG_DATA_SZ = 2, TAG_DATA, TAG_XRANGE_SZ, TAG_XRANGE, TAG_HASSENT };

    template <typename hmesh_t>
    static bool prepare_data_to_send
    (
      const hmesh_t & mesh,
      const std::vector<std::pair<E_Int, std::vector<int>>>& zone2jlist,
      std::map<int, std::map<int, K_FLD::IntArray>> & raczone_to_PG_to_plan
    )
    {
      raczone_to_PG_to_plan.clear();

      bool has_packs{ false };
      //std::cout << "pack : loop on joins : " << join << std::endl;
      for (auto& i : zone2jlist)
      {
        E_Int jid = i.first;
        auto& ptlist = i.second;

        for (size_t i = 0; i < ptlist.size(); ++i)
        {
          K_FLD::IntArray p; //only implemented for IntArray
          E_Int PGi = ptlist[i] - 1;
          //std::cout << "i/PGi : " << i << "/" << PGi << std::endl;
          mesh.extract_plan(PGi, true/*reverse*/, 0/*because previous sync*/, p);
          //std::cout << "after extract_plan" << std::endl;
          if (p.getSize())
          {
            raczone_to_PG_to_plan[jid][PGi] = p;
            has_packs = true;
          }
        }
      }

      return has_packs;

    }

    template <typename hmesh_t, typename sensor_t>
    static void run
    (
      hmesh_t& hmesh, sensor_t& sensor, 
      const std::vector<std::pair<E_Int, std::vector<E_Int>>>& zone2jlist, 
      const std::vector<int>& zonerank, 
      MPI_Comm COM, 
      bool do_agglo
    )
    {
      using adaptor_t = NUGA::adaptor<hmesh_t, sensor_t>;

      adaptor_t::run(hmesh, sensor);

      int rank{ 0 };
      int nranks{ 0 };
      MPI_Comm_size(COM, &nranks);
      MPI_Comm_rank(COM, &rank);

      std::vector<int> jzones;
      jzones.reserve(zone2jlist.size());
      for (auto z : zone2jlist)
        jzones.push_back(z.first);

      bool has_changes{ true };

      while (has_changes)
      {
        has_changes = false;

        //if (rank == 1) std::cout << "rank : " << rank << " prepare_data_to_send... " << std::endl;

        // prepare & send
        std::map<int, std::map<int, K_FLD::IntArray>> jz_to_PG_to_plan;
        bool has_packs = prepare_data_to_send(hmesh, zone2jlist, jz_to_PG_to_plan);

        //if (rank == 1) std::cout << "rank : " << rank << " convert_to_MPI_exchange_format... " << std::endl;

        std::map<int, std::vector<int>> jzrank_to_data, jzrank_to_xrange;
        convert_to_MPI_exchange_format(jz_to_PG_to_plan, zonerank, jzrank_to_data, jzrank_to_xrange);

        //if (rank == 1) std::cout << "rank : " << rank << " send_data... " << std::endl;

        send_data(rank, nranks, COM, jzrank_to_data, jzrank_to_xrange);

        //if (rank == 1) std::cout << "rank : " << rank << " receive_data... " << std::endl;

        MPI_Barrier(COM);

        // receive and process
        receive_data(rank, nranks, COM, jzrank_to_data, jzrank_to_xrange);

        //if (rank == 1) std::cout << "rank : " << rank << " convert_MPI_format_to_sensor_format... " << std::endl;

        std::map<int, K_FLD::IntArray> sensor_data;
        convert_MPI_format_to_sensor_format(jzrank_to_data, jzrank_to_xrange, zone2jlist, sensor_data);

        using join_sensor_t = NUGA::join_sensor<hmesh_t>;
        join_sensor_t jsensor(hmesh);
        jsensor.assign_data(sensor_data);

        //if (rank == 1) std::cout << "rank : " << rank << " has started adapt" << std::endl;

        E_Int npgs0 = hmesh._ng.PGs.size();
        E_Int nphs0 = hmesh._ng.PHs.size();
        NUGA::adaptor<mesh_t, join_sensor_t>::run(hmesh, jsensor, false/*do_agglo*/); //fixme : agglo

        bool local_changes = (hmesh._ng.PGs.size() != npgs0) || (hmesh._ng.PHs.size() != nphs0);

        has_changes = exchange_status(local_changes, rank, nranks, COM);
      }
    }

    ///
    static void convert_to_MPI_exchange_format
    (
      const std::map<int, std::map<int, K_FLD::IntArray>> & jz_to_PG_to_plan,
      const std::vector<int>& zonerank,
      std::map<int, std::vector<int>> & jzrank_to_data,
      std::map<int, std::vector<int>> & jzrank_to_xrange
    )
    {
      jzrank_to_data.clear();
      jzrank_to_xrange.clear();

      for (auto i : jz_to_PG_to_plan)
      {
        std::vector<int> data, xrange;

        auto& dat = i.second;

        for (auto k : dat)
        {
          //int PGi = k.first; //fixme : useless because same sorting (?)
          auto p = k.second; //copy to ensure contiguous
                             //std::cout << "Int arr : r/c : " << p.rows() << "/" << p.cols() << std::endl;

          xrange.push_back(data.size());
          data.insert(data.end(), p.begin(), p.end());
        }
        xrange.push_back(data.size());//one-pass-the-end

        int jzid = i.first;
        assert(jzid > -1 && jzid < zonerank.size());
        int jzrank = zonerank[jzid];

        jzrank_to_data[jzrank] = data;
        jzrank_to_xrange[jzrank] = xrange;
      }
    }

    ///
    static void send_data
    (
      int rank, int nrank, MPI_Comm COM,
      const std::map<int, std::vector<int>>& szone_to_data,
      const std::map<int, std::vector<int>>& szone_to_xrange
    )
    {
      int nszone = szone_to_data.size();
      assert(nszone == szone_to_xrange.size());

      STACK_ARRAY(MPI_Request, nszone, sreqs0); // for buffer size
      STACK_ARRAY(MPI_Request, nszone, sreqs1); // for data
      STACK_ARRAY(MPI_Request, nszone, sreqs2); // for xrange size
      STACK_ARRAY(MPI_Request, nszone, sreqs3); // for xrange

      STACK_ARRAY(bool, nrank, has_sent);
      for (size_t n = 0; n < nrank; ++n) has_sent[n] = false;

      int i = 0;
      for (auto& d : szone_to_data) // WARNING : reference is manadatory ! otherwise ISend might not finish before iteration d copy is deleted
      {
        int zrankid = d.first;
        auto& data = d.second;

        has_sent[zrankid] = true;
 
        auto itxr = szone_to_xrange.find(zrankid);
        assert(itxr != szone_to_xrange.end());

        //std::cout << "envoie de la data sz : " << data.size() << std::endl;
        int data_sz = data.size();
        int err = MPI_Isend(&data_sz, 1, MPI_INT, zrankid, int(TAG_DATA_SZ), COM, &sreqs0[i]);
       /* if (rank == 3 && zid == 2)
        {
          std::cout << "nb de donnees de 3 envoyees a 2 : " << data_sz << std::endl;
          for (size_t k = 0; k < data.size(); ++k)std::cout << data[k] << "/";
          std::cout << std::endl;
        }*/
        
        err = MPI_Isend(&data[0], data.size(), MPI_INT, zrankid, int(TAG_DATA), COM, &sreqs1[i]);
        

        auto& xrange = itxr->second;
        int xrange_sz = xrange.size();
        //std::cout << "envoie de la xrange sz : " << xrange_sz << std::endl;
        err = MPI_Isend(&xrange_sz, 1, MPI_INT, zrankid, int(TAG_XRANGE_SZ), COM, &sreqs2[i]);

        err = MPI_Isend(&xrange[0], xrange_sz, MPI_INT, zrankid, int(TAG_XRANGE), COM, &sreqs3[i]);

        ++i;
      }

      {
        STACK_ARRAY(MPI_Status, nszone, status);
        MPI_Waitall(i, sreqs0.get(), status.get());
        MPI_Waitall(i, sreqs1.get(), status.get());
        MPI_Waitall(i, sreqs2.get(), status.get());
        MPI_Waitall(i, sreqs3.get(), status.get());
      }

      //std::cout << "ALL SENT " << std::endl;

      // tell if sent or not to all ranks (to avoid deadlock when receiving)
      STACK_ARRAY(MPI_Request, nrank, sreq4);
      for (size_t r = 0; r < nrank; ++r) {
        if (r == rank) continue;
        MPI_Isend(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreq4[r]);
      }
      STACK_ARRAY(MPI_Status, nrank, status);
      MPI_Waitall(nrank, sreq4.get(), status.get());

      /*for (size_t r = 0; r < nrank; ++r)
      {
        std::cout << "rank : " << rank << " has sent to " << r << "? : " << has_sent[r] << std::endl;
      }*/

      //std::cout << "HAS SENT OK" << std::endl;

    }
     
    ///
    static void receive_data
    (
      int rank, int nrank, MPI_Comm COM,
      std::map<int, std::vector<int>>& jzrank_to_data,
      std::map<int, std::vector<int>>& jzrank_to_xrange
    )
    {
      jzrank_to_data.clear();
      jzrank_to_xrange.clear();

      STACK_ARRAY(bool, nrank, has_sent);
      has_sent[rank] = false;
      
      STACK_ARRAY(MPI_Request, nrank, sreqs);
      for (size_t r = 0; r < nrank; ++r)
      {
        if (r == rank) continue;
        int err = MPI_Irecv(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreqs[r]);
      }

      {
        STACK_ARRAY(MPI_Status, nrank, status);
        MPI_Waitall(nrank, sreqs.get(), status.get());
      }

      // 1. catch all data sizes
      STACK_ARRAY(int, nrank, data_sz);
      STACK_ARRAY(MPI_Request, nrank, sreqs0); // for buffer size

      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        int err = MPI_Irecv(&data_sz[r], 1, MPI_INT, r, int(TAG_DATA_SZ), COM, &sreqs0[r]);
      }

      {
        STACK_ARRAY(MPI_Status, nrank, status);
        MPI_Waitall(nrank, sreqs0.get(), status.get());
      }

      /*if (rank == 2)
      for (size_t r = 0; r < nrank; ++r)
      {
        std::cout <<r << " has sent to " << rank << " : "  << has_sent[r] << std::endl;
        std::cout << "data_sz : " << data_sz[r] << std::endl;
      }*/

      // 2. get data
      std::vector<int*> data_buffers(nrank, nullptr);
      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        data_buffers[r] = new int[data_sz[r]];
      }

      STACK_ARRAY(MPI_Request, nrank, sreqs1); // for data
      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        int err = MPI_Irecv(data_buffers[r], data_sz[r], MPI_INT, r, int(TAG_DATA), COM, &sreqs1[r]);
      }

      {
        STACK_ARRAY(MPI_Status, nrank, status);
        MPI_Waitall(nrank, sreqs1.get(), status.get());
      }

      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        auto & entry = jzrank_to_data[r];
        entry.insert(entry.end(), data_buffers[r], data_buffers[r] + data_sz[r]);
      }

      // 3. catch all xrange sizes

      STACK_ARRAY(int, nrank, xrange_sz);
      STACK_ARRAY(MPI_Request, nrank, sreqs2); // for buffer size
      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        int err = MPI_Irecv(&xrange_sz[r], 1, MPI_INT, r, int(TAG_XRANGE_SZ), COM, &sreqs2[r]);
      }

      {
        STACK_ARRAY(MPI_Status, nrank, status);
        MPI_Waitall(nrank, sreqs2.get(), status.get());
      }

      // 2. get xrange
      std::vector<int*> xrange_buffers(nrank, nullptr);
      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        xrange_buffers[r] = new int[xrange_sz[r]];
      }

      STACK_ARRAY(MPI_Request, nrank, sreqs3);
      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        int err = MPI_Irecv(xrange_buffers[r], xrange_sz[r], MPI_INT, r, int(TAG_XRANGE), COM, &sreqs3[r]);
      }

      {
        STACK_ARRAY(MPI_Status, nrank, status);
        MPI_Waitall(nrank, sreqs3.get(), status.get());
      }

      for (size_t r = 0; r < nrank; ++r)
      {
        if (has_sent[r] == false) continue;
        auto & entry = jzrank_to_xrange[r];
        entry.insert(entry.end(), xrange_buffers[r], xrange_buffers[r] + xrange_sz[r]);
      }

      for (size_t r = 0; r < nrank; ++r)
      {
        delete data_buffers[r];
        delete xrange_buffers[r];
      }
    }

    static void convert_MPI_format_to_sensor_format
    (
      const std::map<int, std::vector<int>>& rzone_to_data,
      const std::map<int, std::vector<int>>& rzone_to_xrange,
      const std::vector<std::pair<E_Int, std::vector<int>>>& zone2jlist,
      std::map<int, K_FLD::IntArray>& sensor_data
    )
    {
      sensor_data.clear();

      for (auto r2d : rzone_to_data)
      {
        int zid = r2d.first;
        auto& data = r2d.second;

        const std::vector<int>* ptList{ nullptr };
        //std::cout << "zid : " << zid << std::endl;
        //std::cout << "nb of entries : " << zone2jlist.size() << std::endl;
        for (size_t i = 0; i < zone2jlist.size(); ++i)
          if (zone2jlist[i].first == zid)
          {
            ptList = &zone2jlist[i].second;
          }
        assert(ptList != nullptr);

        auto itxr = rzone_to_xrange.find(zid);
        assert(itxr != rzone_to_xrange.end());
        auto & xrange = itxr->second;

        int npgs = xrange.size() - 1;

        assert(npgs == ptList->size());


        for (size_t n = 0; n < npgs; ++n)
        {
          int sz = xrange[n + 1] - xrange[n];
          int cols = sz / 4;
          //std::cout << "sz : " << sz << std::endl;
          //std::cout << "cols : " << cols << std::endl;
          K_FLD::IntArray a;
          a.resize(4, cols);
          int* p = &data[0];
          for (size_t c = 0; c < cols; ++c)
          {
            a(0, c) = *(p++);
            a(1, c) = *(p++);
            a(2, c) = *(p++);
            a(3, c) = *(p++);
          }

          int PGi = (*ptList)[n] - 1;
          //std::cout << "PGI : " << PGi << std::endl;
          sensor_data[PGi] = a;

        }

      }
    }

    static bool exchange_status (bool local_changes, int rank, int nranks, MPI_Comm COM)
    {
      bool has_changes{ false };

      STACK_ARRAY(bool, nranks, all_status);
      MPI_Request req;
      if (rank != 0) {
        MPI_Isend(&local_changes, 1, MPI_C_BOOL, 0, 7, COM, &req);
        //if (rank == 1) std::cout << "rank : " << rank << " is sendind its change status : " << has_changes << std::endl;
      }
      else
      {
        all_status[0] = local_changes;
        MPI_Status status;
        for (size_t r = 1; r < nranks; ++r) {
          MPI_Recv(&all_status[r], 1, MPI_C_BOOL, r, 7, COM, &status);
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
          MPI_Isend(&has_changes, 1, MPI_C_BOOL, r, 7, COM, &req[r - 1]);
        MPI_Waitall(nranks - 1, req.get(), status.get());
        //std::cout << "rank : 0 is telling to all change status : " << has_changes << std::endl;
      }
      else
      {
        MPI_Status status;
        MPI_Recv(&has_changes, 1, MPI_C_BOOL, 0, 7, COM, &status);
        //std::cout << "rank : " << rank << "has been told change status : " << has_changes << std::endl;
      }

      return has_changes;
    }

  };

}

#endif
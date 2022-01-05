/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_ADAPTOR_MPI_HXX
#define NUGA_ADAPTOR_MPI_HXX

#include "Nuga/include/mpi_stl.hxx"
#include "Nuga/include/adaptor.hxx"
#include "Nuga/include/join_sensor.hxx"


namespace NUGA
{
  struct join_msg_type
  {
    std::vector<int> data, datarange, pgs, szone, szonerange, jzone, jzonerange;

    void clear() { data.clear(); datarange.clear(); pgs.clear(); szone.clear(); szonerange.clear(); jzone.clear(); jzonerange.clear(); }
  };

  template <typename mesh_t, typename sensor_t>
  class adaptor_mpi
  {
  public:

    enum eMPI_Tag { TAG_DATA_SZ = 2, TAG_DATA=3, TAG_XRANGE_SZ=4, TAG_XRANGE=5, TAG_PGS_SZ=6, TAG_PGS=7, TAG_SZONE_SZ=8, TAG_SZONERANGE_SZ=10, TAG_JZONE_SZ = 12, TAG_JZONERANGE_SZ = 14, TAG_MPI_EXCH_STATUS=16, TAG_HASSENT=99};

    template <typename hmesh_t>
    static bool prepare_data_to_send
    (
      const hmesh_t & mesh,
      const std::vector<std::pair<E_Int, std::vector<int>>>& zone2jlist,
      std::map<int, std::map<int, K_FLD::IntArray>> & raczone_to_PG_to_plan
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
          K_FLD::IntArray p; //only implemented for IntArray
          E_Int PGi = ptlist[i] - 1;
          //std::cout << "i/PGi : " << i << "/" << PGi << std::endl;
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
    template <typename hmesh_t, typename data_t>
    static void exchange_mpi_data
    (
      const std::vector<hmesh_t>& hmeshes,
      const std::map<int, std::vector<std::pair<E_Int, std::vector<E_Int>>>>& zone_to_zone2jlists,
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
        const auto it = zone_to_zone2jlists.find(hmeshes[i].zid);
        
        if (it == zone_to_zone2jlists.end()) continue; // current zone has no MPI joins

        const auto & zone2jlists = it->second;
        
        std::map<int, std::map<int, K_FLD::IntArray>> jz_to_PG_to_plan;
        has_packs |= prepare_data_to_send(hmeshes[i], zone2jlists, jz_to_PG_to_plan);

        sz_to_jz_to_PG_to_plan[hmeshes[i].zid] = jz_to_PG_to_plan;

      }

      // Convert for sending and gather by rank : rank/zone/face/plan
      std::map<int, join_msg_type> rank_to_mpi_data;
      convert_to_MPI_exchange_format(sz_to_jz_to_PG_to_plan, zonerank, rank_to_mpi_data);

      // Send MPI data   
      send_data(rank, nranks, COM, rank_to_mpi_data);

      // Receive MPI data and build sensor data by zone
      receive_data(rank, nranks, COM, zone_to_zone2jlists, zone_to_sensor_data);
    }

    //
    template <typename hmesh_t, typename data_t>
    static void exchange_omp_data
    (
      const std::vector<hmesh_t>& hmeshes,
      const std::map<int, std::vector<std::pair<E_Int, std::vector<E_Int>>>>& zone_to_zone2jlists,
      std::map<int, data_t>& sensor_data
    )
    {
      
      bool has_packs{ false };
      for (size_t i = 0; i < hmeshes.size(); ++i)
      {
        const auto it = zone_to_zone2jlists.find(hmeshes[i].zid);

        if (it == zone_to_zone2jlists.end()) continue; // current zone has no MPI joins

        const auto & zone2jlists = it->second;

        std::map<int, std::map<int, K_FLD::IntArray>> jz_to_PG_to_plan;
        has_packs |= prepare_data_to_send(hmeshes[i], zone2jlists, jz_to_PG_to_plan);

        // convert to sensor data
        for (auto& j : jz_to_PG_to_plan)
        {
          int jzid = j.first;

          //get the ptlist
          const auto it2 = zone_to_zone2jlists.find(jzid);
          if (it2 == zone_to_zone2jlists.end()) continue;
          const auto & zone2jlists2 = it2->second;
          const std::vector<int>* ptList{ nullptr };
          for (size_t z = 0; z < zone2jlists2.size(); ++z)
            if (zone2jlists2[z].first == hmeshes[i].zid)
            {
              ptList = &zone2jlists2[z].second;
            }
          assert(ptList != nullptr);

          auto & PG_to_plan = j.second;

          for (auto & k : PG_to_plan)
          {
            E_Int PGi = (*ptList)[k.first] - 1;
            sensor_data[jzid][PGi] = k.second;
          }
        }
      }
    }

    template <typename hmesh_t, typename sensor_t>
    static void run
    (
      std::vector<hmesh_t>& hmeshes, std::vector<sensor_t*>& sensors,
      const std::map<int, std::vector<std::pair<E_Int, std::vector<E_Int>>>>& zone_to_zone2jlists,
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      bool do_agglo
    )
    {
      E_Int err(0);
      E_Int NBZ{ E_Int(hmeshes.size()) };

      assert(NBZ == sensors.size());

      using adaptor_t = NUGA::adaptor<hmesh_t, sensor_t>;
      using data_t = std::map<int, K_FLD::IntArray>;
      using join_sensor_t = NUGA::join_sensor<hmesh_t>;

      ePara PARA = COARSE_OMP;
#pragma omp parallel for if(PARA == COARSE_OMP)
      for (E_Int i = 0; i < NBZ; ++i)
      {
        NUGA::adaptor<mesh_t, sensor_t>::run(hmeshes[i], *sensors[i], do_agglo);
      }

      int nranks{ 0 };
      MPI_Comm_size(COM, &nranks);
      int rank{ 0 };
      MPI_Comm_rank(COM, &rank);

      for (E_Int i = 0; i < NBZ; ++i)
        assert(rank == zonerank[hmeshes[i].zid]);

      //separate MPI/OMP joins
      std::map<int, std::vector<std::pair<E_Int, std::vector<E_Int>>>> zone_to_zone2jlists_mpi, zone_to_zone2jlists_omp;
      for (auto & it : zone_to_zone2jlists)
      {
        int zid = it.first;
        auto & zones2jlist = it.second;
        for (auto & it2 : zones2jlist)
        {
          int jzid = it2.first;
          auto& ptlist = it2.second;
          int rk = zonerank[jzid];
          if (rk == rank) // local join
            zone_to_zone2jlists_omp[zid].push_back(std::make_pair(jzid, ptlist));
          else // distant
            zone_to_zone2jlists_mpi[zid].push_back(std::make_pair(jzid, ptlist));
        }
      }

      bool has_mpi_changes{ true };

      int mpi_iter = -1;

      while (has_mpi_changes)
      {
        ++mpi_iter;
        //std::cout << "rank : " << rank << " mpi iter : " << mpi_iter << std::endl;
        //if (mpi_iter == 1) return;
        has_mpi_changes = false;

        std::map<int, data_t> zone_to_sensor_data;
        exchange_mpi_data(hmeshes, zone_to_zone2jlists_mpi, zonerank, COM, rank, nranks, zone_to_sensor_data);
        
        bool has_omp_changes{ true }, has_local_changes{ false };
        int omp_iter = -1;

        while (has_omp_changes)
        {
          ++omp_iter;
          has_omp_changes = false;

          exchange_omp_data(hmeshes, zone_to_zone2jlists_omp, zone_to_sensor_data);

          // Adapt each zone
#pragma omp parallel for reduction ( || : has_omp_changes, has_local_changes) if(PARA == COARSE_OMP)
          for (E_Int i = 0; i < NBZ; ++i)
          {
            join_sensor_t jsensor(hmeshes[i]);
            jsensor.assign_data(zone_to_sensor_data[hmeshes[i].zid]);

            E_Int npgs0 = hmeshes[i]._ng.PGs.size();
            E_Int nphs0 = hmeshes[i]._ng.PHs.size();

            NUGA::adaptor<mesh_t, join_sensor_t>::run(hmeshes[i], jsensor, false/*do_agglo*/); //fixme : agglo

            has_omp_changes |= (hmeshes[i]._ng.PGs.size() != npgs0) || (hmeshes[i]._ng.PHs.size() != nphs0);
            has_local_changes |= has_omp_changes;
          }
        }
        MPI_Barrier(COM);
        has_mpi_changes = exchange_status(has_local_changes, rank, nranks, COM);
      }
    }

    ///
    static void convert_to_MPI_exchange_format
    (
      const std::map<int, std::map<int, std::map<int, K_FLD::IntArray>>> & sz_to_jz_to_PG_to_plan,
      const std::vector<int>& zonerank,
      std::map<int,join_msg_type> & rank_to_data
    )
    {
      rank_to_data.clear();

      // Gather data by rank
      std::map<int, std::map<int, std::map<int, std::map<int, K_FLD::IntArray>>>> rank_to_sz_to_jz_to_PG_to_plan;
      for (auto & it : sz_to_jz_to_PG_to_plan)
      {
        int szid = it.first;
        const auto& jz_to_PG_to_plan = it.second;

        for (auto& it2 : jz_to_PG_to_plan)
        {

          int jzid = it2.first;
          auto& PG_to_plan = it2.second;
          assert(jzid > -1 && jzid < zonerank.size());
          int jrk = zonerank[jzid];

          rank_to_sz_to_jz_to_PG_to_plan[jrk][szid][jzid] = PG_to_plan;
        }
      }
      
      //
      for (auto& i : rank_to_sz_to_jz_to_PG_to_plan)
      {
        int rank = i.first;
        auto & sz_to_jz_to_PG_to_plan_PG_to_plan = i.second;

        auto & rankdata = rank_to_data[rank];

        for (auto & kk : sz_to_jz_to_PG_to_plan_PG_to_plan)
        {
          int szid = kk.first;
          auto & jz_to_PG_to_plan = kk.second;

          rankdata.szonerange.push_back(rankdata.jzone.size());
          rankdata.szone.push_back(szid);

          for (auto& j : jz_to_PG_to_plan)
          {
            int jzid = j.first;
            auto& PG_to_plan = j.second;

            rankdata.jzonerange.push_back(rankdata.pgs.size());
            rankdata.jzone.push_back(jzid);

            for (auto& k : PG_to_plan)
            {
              int iloc/*PGi*/ = k.first;
              auto p = k.second; //copy to ensure contiguous

              assert(p.cols() != 0);
              int sz0 = rankdata.data.size();
              rankdata.datarange.push_back(rankdata.data.size());
              rankdata.data.insert(rankdata.data.end(), p.begin(), p.end());
              assert(rankdata.data.size() > sz0);
              rankdata.pgs.push_back(iloc);
            }
          }
        }

        //one-pass-the-end
        rankdata.datarange.push_back(rankdata.data.size());
        rankdata.jzonerange.push_back(rankdata.pgs.size());
        rankdata.szonerange.push_back(rankdata.jzone.size());
      }
    }

    ///
    static void send_data
    (
      int rank, int nrank,
      MPI_Comm COM,
      const std::map<int, join_msg_type> & rank_to_data
    )
    {
      int nsranks = rank_to_data.size(); // nb of ranks to send to
      int NB_TOT_REQS = 14 * nsranks;    // for 7 vectors : data/datarange/pgs/szone/szonerange/jzone/jzonerange.  2 req per vec. => 14

      STACK_ARRAY(MPI_Request, NB_TOT_REQS, sreqs);

      STACK_ARRAY(bool, nrank/*all ranks*/, has_sent);
      for (size_t n = 0; n < nrank; ++n) has_sent[n] = false;
      
      int count_req{ -1 };
      for (auto& d : rank_to_data) // WARNING : reference is manadatory ! otherwise ISend might not finish before iteration d copy is deleted
      {
        int rankid = d.first;
        auto& data = d.second;

        has_sent[rankid] = true;

        NUGA::MPI::Isend(data.data, rankid, TAG_DATA_SZ, COM, &sreqs[++count_req], &sreqs[++count_req]);

        NUGA::MPI::Isend(data.datarange, rankid, TAG_XRANGE_SZ, COM, &sreqs[++count_req], &sreqs[++count_req]);

        NUGA::MPI::Isend(data.pgs, rankid, TAG_PGS_SZ, COM, &sreqs[++count_req], &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szone, rankid, TAG_SZONE_SZ, COM, &sreqs[++count_req], &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szonerange, rankid, TAG_SZONERANGE_SZ, COM, &sreqs[++count_req], &sreqs[++count_req]);

        NUGA::MPI::Isend(data.jzone, rankid, TAG_JZONE_SZ, COM, &sreqs[++count_req], &sreqs[++count_req]);

        NUGA::MPI::Isend(data.jzonerange, rankid, TAG_JZONERANGE_SZ, COM, &sreqs[++count_req], &sreqs[++count_req]);
      }

      ++count_req;
      if (count_req > 0)
      {
        STACK_ARRAY(MPI_Status, count_req, status);
        MPI_Waitall(count_req, sreqs.get(), status.get());
      }

      //std::cout << "ALL SENT " << std::endl;

      // tell if sent or not to all ranks (to avoid deadlock when receiving)
      STACK_ARRAY(MPI_Request, nrank, sreq);
      count_req = -1;
      for (size_t r = 0; r < nrank; ++r) {
        //std::cout << "rank : " << r << " has sent ? : " << has_sent[r] << std::endl;
        if (r == rank) continue;
        MPI_Isend(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreq[++count_req]);
      }

      STACK_ARRAY(MPI_Status, nrank-1, stats);
      MPI_Waitall(nrank-1, sreq.get(), stats.get());
    }
     
   
    ///
    template <typename data_t>
    static int receive_data
    (
      int rank, int nranks, MPI_Comm COM,
      const std::map<int, std::vector<std::pair<E_Int, std::vector<E_Int>>>>& zone_to_zone2jlists,
      std::map<int, data_t>& zone_to_sensor_data
    )
    {
      zone_to_sensor_data.clear();

      STACK_ARRAY(bool, nranks, has_sent);
      has_sent[rank] = false;

      int NB_REQ_FOR_HASSENT = nranks - 1;
      {
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_HASSENT, sreqs);
        int req_count = -1;
        for (size_t r = 0; r < nranks; ++r)
        {
          if (r == rank) continue;
          int err = MPI_Irecv(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreqs[++req_count]);
        }

        if (NB_REQ_FOR_HASSENT > 0)
        {
          STACK_ARRAY(MPI_Status, NB_REQ_FOR_HASSENT, status);
          MPI_Waitall(NB_REQ_FOR_HASSENT, sreqs.get(), status.get());
        }
      }

      // 1. catch all sizes :
      int NB_SENT = 0;
      for (size_t r = 0; r < nranks; ++r) if (has_sent[r])++NB_SENT;

      int NB_REQ_FOR_SIZES = 7 * NB_SENT;
      STACK_ARRAY(MPI_Request, NB_REQ_FOR_SIZES, sreqs_sz);
      // 1.1. data
      STACK_ARRAY(int, NB_SENT, data_sz);
      // 1.2. datarange
      STACK_ARRAY(int, NB_SENT, datarange_sz);
      // 1.3. pgs
      STACK_ARRAY(int, NB_SENT, pgs_sz);
      // 1.4. szone
      STACK_ARRAY(int, NB_SENT, szone_sz);
      // 1.5. szonerange
      STACK_ARRAY(int, NB_SENT, szonerange_sz);
      // 1.4. szone
      STACK_ARRAY(int, NB_SENT, jzone_sz);
      // 1.5. szonerange
      STACK_ARRAY(int, NB_SENT, jzonerange_sz);

      int req_count = -1;
      for (size_t r = 0; r < nranks; ++r)
      {
        if (has_sent[r] == false) continue;

        int err = MPI_Irecv(&data_sz[r], 1, MPI_INT, r, int(TAG_DATA_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&datarange_sz[r], 1, MPI_INT, r, int(TAG_XRANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;
        
        err = MPI_Irecv(&pgs_sz[r], 1, MPI_INT, r, int(TAG_PGS_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szone_sz[r], 1, MPI_INT, r, int(TAG_SZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szonerange_sz[r], 1, MPI_INT, r, int(TAG_SZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;
        
        err = MPI_Irecv(&jzone_sz[r], 1, MPI_INT, r, int(TAG_JZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&jzonerange_sz[r], 1, MPI_INT, r, int(TAG_JZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;
      }
      
      assert((req_count + 1) == NB_REQ_FOR_SIZES);
      if (NB_REQ_FOR_SIZES > 0)
      {
        //std::cout << "req count vs nranks vs NB_REQ_FOR_SIZES : " << req_count << "/" << nranks << "/" << NB_REQ_FOR_SIZES << std::endl;
        STACK_ARRAY(MPI_Status, NB_REQ_FOR_SIZES, status);
        MPI_Waitall(NB_REQ_FOR_SIZES, sreqs_sz.get(), status.get());
      }


      for (size_t r = 0; r < nranks; ++r)
      {
        if (has_sent[r] == false) continue;
        //std::cout << "rank : " << rank << " has_sent[" << r << "] : " << has_sent[r] << std::endl;
        //std::cout << "rank : " << rank << " datarange_sz[" << r << "] : " << datarange_sz[r] << std::endl;
        //std::cout << "rank : " << rank << " data_sz[" << r << "] : " << data_sz[r] << std::endl;

        assert(data_sz[r] != 0);
        assert(datarange_sz[r] != 0);
        assert(pgs_sz[r] != 0);
        assert(szone_sz[r] != 0);
        assert(szonerange_sz[r] != 0);
        assert(jzone_sz[r] != 0);
        assert(jzonerange_sz[r] != 0);
      }

      // 2. get data
      std::map<int, join_msg_type> rank_to_data;
      {
        int NB_REQ_FOR_DATA = 7 * NB_SENT;
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_DATA, sreqs_data);

        req_count = -1;
        for (size_t r = 0; r < nranks; ++r)
        {
          if (has_sent[r] == false) continue;

          int err = NUGA::MPI::Irecv(data_sz[r], r, int(TAG_DATA_SZ) + 1, COM, rank_to_data[r].data, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(datarange_sz[r], r, int(TAG_XRANGE_SZ) + 1, COM, rank_to_data[r].datarange, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(pgs_sz[r], r, int(TAG_PGS_SZ) + 1, COM, rank_to_data[r].pgs, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szone_sz[r], r, int(TAG_SZONE_SZ) + 1, COM, rank_to_data[r].szone, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szonerange_sz[r], r, int(TAG_SZONERANGE_SZ) + 1, COM, rank_to_data[r].szonerange, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(jzone_sz[r], r, int(TAG_JZONE_SZ) + 1, COM, rank_to_data[r].jzone, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(jzonerange_sz[r], r, int(TAG_JZONERANGE_SZ) + 1, COM, rank_to_data[r].jzonerange, &sreqs_data[++req_count]);
          if (err) return 1;
        }

        ++req_count;
        if (req_count > 0)
        {
          //std::cout << "req count vs nranks vs NB_REQ_FOR_SIZES : " << req_count << "/" << nranks << "/" << NB_REQ_FOR_SIZES << std::endl;
          STACK_ARRAY(MPI_Status, req_count, status);
          MPI_Waitall(req_count, sreqs_data.get(), status.get());
        }
      }
      

      // 3. dispatch sensor data by zone
      for (auto & i : rank_to_data)
      {
        //int rankid = i.first;
        const join_msg_type& msg = i.second;

        const auto& szonerange = msg.szonerange;
        const auto& szone      = msg.szone;
        const auto& jzonerange = msg.jzonerange;
        const auto& jzone      = msg.jzone;
        const auto& datarange  = msg.datarange;
        const auto& data       = msg.data;
        const auto & pgs       = msg.pgs;

        assert(pgs.size() == datarange.size() - 1);
        assert(szone.size() == szonerange.size() - 1);
        assert(jzone.size() == jzonerange.size() - 1);

        for (size_t u = 0; u < szone.size(); ++u)
        {
          assert(u >= 0 && u < szone.size());
          int szid = szone[u];
          //if (rank == 1) std::cout << "szid : " << szid << std::endl;

          int jzonerange_beg = szonerange[u];
          int jzonerange_end = szonerange[u + 1];

          //if (rank == 1) std::cout << "jzonerange_beg : " << jzonerange_beg << std::endl;
          //if (rank == 1) std::cout << "jzonerange_end : " << jzonerange_end << std::endl;

          for (size_t v = jzonerange_beg; v < jzonerange_end; ++v)
          {
            assert(v >= 0 && v < jzone.size());
            int jzid = jzone[v];

            int n_beg = jzonerange[v];
            int n_end = jzonerange[v + 1];

            //if (rank == 1) std::cout << "n_beg : " << jzone_beg << std::endl;
            //if (rank == 1) std::cout << "n_end : " << jzone_end << std::endl;

            const auto it = zone_to_zone2jlists.find(jzid);
            assert(it != zone_to_zone2jlists.end());
            const auto & zone2jlists = it->second;
            const std::vector<int>* ptList{ nullptr };
            //std::cout << "zid : " << zid << std::endl;
            //std::cout << "nb of entries : " << zone2jlist.size() << std::endl;
            for (size_t i = 0; i < zone2jlists.size(); ++i)
              if (zone2jlists[i].first == szid)
              {
                ptList = &zone2jlists[i].second;
              }
            assert(ptList != nullptr);

            for (size_t n = n_beg; n < n_end; ++n)
            {
              assert(n >=0 && n < datarange.size());
              assert((n+1) < datarange.size());
              int data_beg = datarange[n];
              int data_end = datarange[n + 1];
              assert(n < pgs.size());
              int iloc = pgs[n];

              int sz = data_end - data_beg;
              int cols = sz / 4;
              //std::cout << "sz : " << sz << std::endl;
              //std::cout << "cols : " << cols << std::endl;

              assert(data_beg >= 0 && data_beg < data.size());
              const int* p = &data[data_beg];
              K_FLD::IntArray a;
              a.resize(4, cols);
              for (size_t c = 0; c < cols; ++c)
              {
                a(0, c) = *(p++);
                a(1, c) = *(p++);
                a(2, c) = *(p++);
                a(3, c) = *(p++);
              }
              assert(iloc >= 0 && iloc < ptList->size());
              int PGi = (*ptList)[iloc]-1;
              //std::cout << "PGi : " << PGi << std::endl;
              //std::cout << a << std::endl;

              zone_to_sensor_data[jzid][PGi] = a;
            }
          }

        }
      }
      
      return 0;
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
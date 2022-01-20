/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_MPI_MSG_TYPES_HXX
#define NUGA_MPI_MSG_TYPES_HXX

#include "Nuga/include/mpi_stl.hxx"


namespace NUGA
{
  enum eMPI_Tag { TAG_DATA_SZ = 2, TAG_DATA = 3, TAG_XRANGE_SZ = 4, TAG_XRANGE = 5, TAG_PGS_SZ = 6, TAG_PGS = 7, TAG_SZONE_SZ = 8, TAG_SZONERANGE_SZ = 10, TAG_JZONE_SZ = 12, TAG_JZONERANGE_SZ = 14, TAG_PTL_SZ = 16, TAG_MPI_EXCH_STATUS = 18, TAG_HASSENT = 99 };

  ///
  inline
  static void split_mpi_omp_joins
  (
      const std::map<int, std::map<int, std::vector<int>>>& zone_to_zone2jlists,
      int rank, const std::vector<int>& zonerank,
      std::map<int, std::map<int, std::vector<int>>>& zone_to_zone2jlists_omp,
      std::map<int, std::map<int, std::vector<int>>>& zone_to_zone2jlists_mpi

   )
  {
    zone_to_zone2jlists_omp.clear();
    zone_to_zone2jlists_mpi.clear();

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
          zone_to_zone2jlists_omp[zid][jzid] = ptlist;
        else // distant
          zone_to_zone2jlists_mpi[zid][jzid] = ptlist;
      }
    }
  }
  
  struct plan_msg_type
  {
    std::vector<int> data, datarange, pgs, szone, szonerange, jzone, jzonerange;

    void clear() { data.clear(); datarange.clear(); pgs.clear(); szone.clear(); szonerange.clear(); jzone.clear(); jzonerange.clear(); }

    ///
    static void convert_to_MPI_exchange_format
    (
      const std::map<int, std::map<int, std::map<int, K_FLD::IntArray>>> & sz_to_jz_to_PG_to_plan,
      const std::vector<int>& zonerank,
      std::map<int, plan_msg_type> & rank_to_data
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
      const std::map<int, plan_msg_type> & rank_to_data
    )
    {
      //if (rank==3) std::cout << "send_data 1" << std::endl;
      int nsranks = rank_to_data.size(); // nb of ranks to send to
      int NB_TOT_REQS = 7 * nsranks;    // for 7 vectors : data/datarange/pgs/szone/szonerange/jzone/jzonerange.  2 req per vec. => 14

      STACK_ARRAY(MPI_Request, NB_TOT_REQS, sreqs);

      STACK_ARRAY(bool, nrank/*all ranks*/, has_sent);
      for (size_t n = 0; n < nrank; ++n) has_sent[n] = false;

      int count_req{ -1 };
      for (auto& d : rank_to_data) // WARNING : reference is manadatory ! otherwise ISend might not finish before iteration d copy is deleted
      {
        int rankid = d.first;
        auto& data = d.second;

        has_sent[rankid] = true;

        NUGA::MPI::Isend(data.data, rankid, TAG_DATA_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.datarange, rankid, TAG_XRANGE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.pgs, rankid, TAG_PGS_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szone, rankid, TAG_SZONE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szonerange, rankid, TAG_SZONERANGE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.jzone, rankid, TAG_JZONE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.jzonerange, rankid, TAG_JZONERANGE_SZ, COM, &sreqs[++count_req]);
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
        if (r == rank) continue;
        //std::cout << "rank : " << rank  << " over (" << nrank << ") has sent to " << r << " ? : " << has_sent[r] << std::endl;
        MPI_Isend(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreq[++count_req]);
      }

      STACK_ARRAY(MPI_Status, nrank - 1, stats);
      MPI_Waitall(nrank - 1, sreq.get(), stats.get());
    }


    ///
    template <typename data_t>
    static int receive_data
    (
      int rank, int nranks, MPI_Comm COM,
      const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone2jlists,
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
      int rcount = -1;
      for (size_t r = 0; r < nranks; ++r)
      {
        if (has_sent[r] == false) continue;
        ++rcount;

        int err = MPI_Irecv(&data_sz[rcount], 1, MPI_INT, r, int(TAG_DATA_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&datarange_sz[rcount], 1, MPI_INT, r, int(TAG_XRANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&pgs_sz[rcount], 1, MPI_INT, r, int(TAG_PGS_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szone_sz[rcount], 1, MPI_INT, r, int(TAG_SZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szonerange_sz[rcount], 1, MPI_INT, r, int(TAG_SZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&jzone_sz[rcount], 1, MPI_INT, r, int(TAG_JZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&jzonerange_sz[rcount], 1, MPI_INT, r, int(TAG_JZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;
      }

      assert((req_count + 1) == NB_REQ_FOR_SIZES);
      if (NB_REQ_FOR_SIZES > 0)
      {
        //std::cout << "req count vs nranks vs NB_REQ_FOR_SIZES : " << req_count << "/" << nranks << "/" << NB_REQ_FOR_SIZES << std::endl;
        STACK_ARRAY(MPI_Status, NB_REQ_FOR_SIZES, status);
        MPI_Waitall(NB_REQ_FOR_SIZES, sreqs_sz.get(), status.get());
      }


      rcount=-1;
      for (size_t r = 0; r < nranks; ++r)
      {
        if (has_sent[r] == false) continue;
        ++rcount;
        //std::cout << "rank : " << rank << " has_sent[" << r << "] : " << has_sent[r] << std::endl;
        //std::cout << "rank : " << rank << " datarange_sz[" << r << "] : " << datarange_sz[r] << std::endl;
        //std::cout << "rank : " << rank << " data_sz[" << r << "] : " << data_sz[r] << std::endl;

        assert(data_sz[rcount] != 0);
        assert(datarange_sz[rcount] != 0);
        assert(pgs_sz[rcount] != 0);
        assert(szone_sz[rcount] != 0);
        assert(szonerange_sz[rcount] != 0);
        assert(jzone_sz[rcount] != 0);
        assert(jzonerange_sz[rcount] != 0);
      }

      // 2. get data
      std::map<int, plan_msg_type> rank_to_data;
      {
        int NB_REQ_FOR_DATA = 7 * NB_SENT;
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_DATA, sreqs_data);

        req_count = -1;
        rcount=-1;
        for (size_t r = 0; r < nranks; ++r)
        {
          if (has_sent[r] == false) continue;
          ++rcount;

          int err = NUGA::MPI::Irecv(data_sz[rcount], r, int(TAG_DATA_SZ) + 1, COM, rank_to_data[r].data, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(datarange_sz[rcount], r, int(TAG_XRANGE_SZ) + 1, COM, rank_to_data[r].datarange, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(pgs_sz[rcount], r, int(TAG_PGS_SZ) + 1, COM, rank_to_data[r].pgs, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szone_sz[rcount], r, int(TAG_SZONE_SZ) + 1, COM, rank_to_data[r].szone, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szonerange_sz[rcount], r, int(TAG_SZONERANGE_SZ) + 1, COM, rank_to_data[r].szonerange, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(jzone_sz[rcount], r, int(TAG_JZONE_SZ) + 1, COM, rank_to_data[r].jzone, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(jzonerange_sz[rcount], r, int(TAG_JZONERANGE_SZ) + 1, COM, rank_to_data[r].jzonerange, &sreqs_data[++req_count]);
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
        const plan_msg_type& msg = i.second;

        const auto& szonerange = msg.szonerange;
        const auto& szone = msg.szone;
        const auto& jzonerange = msg.jzonerange;
        const auto& jzone = msg.jzone;
        const auto& datarange = msg.datarange;
        const auto& data = msg.data;
        const auto & pgs = msg.pgs;

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

          // if (rank == 2) std::cout << "jzonerange_beg : " << jzonerange_beg << std::endl;
          // if (rank == 2) std::cout << "jzonerange_end : " << jzonerange_end << std::endl;
          // if (rank == 2) std::cout << "receive_data : 10 c 3" << std::endl;

          for (size_t v = jzonerange_beg; v < jzonerange_end; ++v)
          {
            assert(v >= 0 && v < jzone.size());
            int jzid = jzone[v];

            int n_beg = jzonerange[v];
            int n_end = jzonerange[v + 1];

            // if (rank == 2) std::cout << "n_beg : " << n_beg << std::endl;
            // if (rank == 2) std::cout << "n_end : " << n_end << std::endl;

            const auto it = zone_to_zone2jlists.find(jzid);
            assert(it != zone_to_zone2jlists.end());
            const auto & zone2jlists = it->second;

            const auto& itPtList = zone2jlists.find(szid);
            assert(itPtList != zone2jlists.end());
            const auto& ptlist = itPtList->second;

            for (size_t n = n_beg; n < n_end; ++n)
            {
              assert(n >= 0 && n < datarange.size());
              assert((n + 1) < datarange.size());
              int data_beg = datarange[n];
              int data_end = datarange[n + 1];
              assert(n < pgs.size());
              int iloc = pgs[n];

              int sz = data_end - data_beg;
              int cols = sz / 4;
              // if (rank == 2) std::cout << "sz : " << sz << std::endl;
              // if (rank == 2) std::cout << "cols : " << cols << std::endl;

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
              assert(iloc >= 0 && iloc < ptlist.size());
              int PGi = ptlist[iloc] - 1;
              //std::cout << "PGi : " << PGi << std::endl;
              //std::cout << a << std::endl;

              zone_to_sensor_data[jzid][PGi] = a;
            }
          }

        }
      }

      return 0;
    }
  };

  struct pointlist_msg_type
  {
    std::vector<int> ptList, szone, szonerange, jzone, jzonerange;

    void clear() { ptList.clear(); szone.clear(); szonerange.clear(); jzone.clear(); jzonerange.clear(); }

    ///
    static void exchange_pointlists
    (
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      int rank, int nranks,
      const std::map<int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone_to_list_owned,
      std::map<int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone_to_list_opp
    )
    {
      zone_to_zone_to_list_opp.clear();

      //separate MPI/OMP joins
      std::map<int, std::map<E_Int, std::vector<E_Int>>> zone_to_zone2jlists_mpi, zone_to_zone2jlists_omp;
      split_mpi_omp_joins(zone_to_zone_to_list_owned, rank, zonerank, zone_to_zone2jlists_omp, zone_to_zone2jlists_mpi);

      // OMP : just transpose
      if (!zone_to_zone2jlists_omp.empty())
      {
        for (const auto& it : zone_to_zone2jlists_omp)
        {
          int zid = it.first;
          const auto& zone2jlists = it.second;
          for (const auto & z2L : zone2jlists)
          {
            int jzid = z2L.first;
            const auto & ptL = z2L.second;

            zone_to_zone_to_list_opp[jzid][zid] = ptL;
          }
        }
      }

      // MPI
      if (zone_to_zone2jlists_mpi.empty()) return;

      // prepare data to send : rank/zone/ptL
      std::map<int, pointlist_msg_type> rank_to_mpi_data;
      convert_to_MPI_exchange_format(zone_to_zone2jlists_mpi, zonerank, rank_to_mpi_data);

      // Send MPI data   
      send_data(rank, nranks, COM, rank_to_mpi_data);

      // Receive MPI data and build sensor data by zone
      receive_data(rank, nranks, COM, zone_to_zone2jlists_mpi, zone_to_zone_to_list_opp);

      assert(zone_to_zone_to_list_owned.size() == zone_to_zone_to_list_opp.size());
    }

    
    ///
    static void convert_to_MPI_exchange_format
    (
      const std::map<int, std::map<E_Int, std::vector<E_Int>>> & zone_to_zone2jlists_mpi,
      const std::vector<int>& zonerank,
      std::map<int, pointlist_msg_type> & rank_to_data
    )
    {
      rank_to_data.clear();

      // Gather data by rank
      std::map<int, std::map<int, std::map<E_Int, std::vector<E_Int>>>> rank_to_sz_to_jz_to_PTL;
      for (auto & it : zone_to_zone2jlists_mpi)
      {
        int szid = it.first;
        const auto& jz_to_PTL = it.second;

        for (auto& it2 : jz_to_PTL)
        {
          int jzid = it2.first;
          auto& PTL = it2.second;

          assert(jzid > -1 && jzid < zonerank.size());
          int jrk = zonerank[jzid];

          rank_to_sz_to_jz_to_PTL[jrk][szid][jzid] = PTL;
        }
      }

      //
      for (auto& i : rank_to_sz_to_jz_to_PTL)
      {
        int rank = i.first;
        auto & sz_to_jz_to_PTL = i.second;

        auto & rankdata = rank_to_data[rank];

        for (auto & kk : sz_to_jz_to_PTL)
        {
          int szid = kk.first;
          auto & jz_to_PTL = kk.second;

          rankdata.szonerange.push_back(rankdata.jzone.size());
          rankdata.szone.push_back(szid);

          for (auto& j : jz_to_PTL)
          {
            int jzid = j.first;
            auto& PTL = j.second;

            rankdata.jzonerange.push_back(rankdata.ptList.size());
            rankdata.jzone.push_back(jzid);

            rankdata.ptList.insert(rankdata.ptList.end(), ALL(PTL));
          }
        }

        //one-pass-the-end
        rankdata.jzonerange.push_back(rankdata.ptList.size());
        rankdata.szonerange.push_back(rankdata.jzone.size());
      }
    }

    ///
    static void send_data
    (
      int rank, int nrank,
      MPI_Comm COM,
      const std::map<int, pointlist_msg_type> & rank_to_data
    )
    {
      int nsranks = rank_to_data.size(); // nb of ranks to send to
      int NB_TOT_REQS = 5 * nsranks;    // for 5 vectors : ptL/szone/szonerange/jzone/jzonerange. 

      STACK_ARRAY(MPI_Request, NB_TOT_REQS, sreqs);

      STACK_ARRAY(bool, nrank, has_sent);
      for (size_t n = 0; n < nrank; ++n) has_sent[n] = false;
      
      int count_req{ -1 };
      for (auto& d : rank_to_data) // WARNING : reference is manadatory ! otherwise ISend might not finish before iteration d copy is deleted
      {
        int rankid = d.first;
        auto& data = d.second;

        has_sent[rankid] = true;

        NUGA::MPI::Isend(data.ptList, rankid, TAG_PTL_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szone, rankid, TAG_SZONE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szonerange, rankid, TAG_SZONERANGE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.jzone, rankid, TAG_JZONE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.jzonerange, rankid, TAG_JZONERANGE_SZ, COM, &sreqs[++count_req]);
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

      STACK_ARRAY(MPI_Status, nrank - 1, stats);
      MPI_Waitall(nrank - 1, sreq.get(), stats.get());
    }

    ///
    static int receive_data
    (
      int rank, int nranks, MPI_Comm COM,
      const std::map<E_Int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone2jlists,
      std::map<int, std::map<E_Int, std::vector<E_Int>>>& zone_to_zone_to_list_opp
    )
    {
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

      int NB_REQ_FOR_SIZES = 5 * NB_SENT;
      STACK_ARRAY(MPI_Request, NB_REQ_FOR_SIZES, sreqs_sz);
      // 1.1. data
      STACK_ARRAY(int, NB_SENT, ptList_sz);
      // 1.4. szone
      STACK_ARRAY(int, NB_SENT, szone_sz);
      // 1.5. szonerange
      STACK_ARRAY(int, NB_SENT, szonerange_sz);
      // 1.4. szone
      STACK_ARRAY(int, NB_SENT, jzone_sz);
      // 1.5. szonerange
      STACK_ARRAY(int, NB_SENT, jzonerange_sz);

      int req_count = -1;
      int rcount = -1;
      for (size_t r = 0; r < nranks; ++r)
      {
        if (has_sent[r] == false) continue;
        ++rcount;

        int err = MPI_Irecv(&ptList_sz[rcount], 1, MPI_INT, r, int(TAG_PTL_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szone_sz[rcount], 1, MPI_INT, r, int(TAG_SZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szonerange_sz[rcount], 1, MPI_INT, r, int(TAG_SZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&jzone_sz[rcount], 1, MPI_INT, r, int(TAG_JZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&jzonerange_sz[rcount], 1, MPI_INT, r, int(TAG_JZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;
      }

      assert((req_count + 1) == NB_REQ_FOR_SIZES);
      if (NB_REQ_FOR_SIZES > 0)
      {
        //std::cout << "req count vs nranks vs NB_REQ_FOR_SIZES : " << req_count << "/" << nranks << "/" << NB_REQ_FOR_SIZES << std::endl;
        STACK_ARRAY(MPI_Status, NB_REQ_FOR_SIZES, status);
        MPI_Waitall(NB_REQ_FOR_SIZES, sreqs_sz.get(), status.get());
      }

      rcount = -1;
      for (size_t r = 0; r < nranks; ++r)
      {
        if (has_sent[r] == false) continue;
        ++rcount;
        //std::cout << "rank : " << rank << " has_sent[" << r << "] : " << has_sent[r] << std::endl;
        //std::cout << "rank : " << rank << " datarange_sz[" << r << "] : " << datarange_sz[r] << std::endl;
        //std::cout << "rank : " << rank << " data_sz[" << r << "] : " << data_sz[r] << std::endl;

        assert(ptList_sz[rcount] != 0);
        assert(szone_sz[rcount] != 0);
        assert(szonerange_sz[rcount] != 0);
        assert(jzone_sz[rcount] != 0);
        assert(jzonerange_sz[rcount] != 0);
      }

      // 2. get data
      std::map<int, pointlist_msg_type> rank_to_data;
      {
        int NB_REQ_FOR_DATA = 5 * NB_SENT;
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_DATA, sreqs_data);

        req_count = -1;
        rcount = -1;
        for (size_t r = 0; r < nranks; ++r)
        {
          if (has_sent[r] == false) continue;
          ++rcount;

          int err = NUGA::MPI::Irecv(ptList_sz[rcount], r, int(TAG_PTL_SZ) + 1, COM, rank_to_data[r].ptList, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szone_sz[rcount], r, int(TAG_SZONE_SZ) + 1, COM, rank_to_data[r].szone, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szonerange_sz[rcount], r, int(TAG_SZONERANGE_SZ) + 1, COM, rank_to_data[r].szonerange, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(jzone_sz[rcount], r, int(TAG_JZONE_SZ) + 1, COM, rank_to_data[r].jzone, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(jzonerange_sz[rcount], r, int(TAG_JZONERANGE_SZ) + 1, COM, rank_to_data[r].jzonerange, &sreqs_data[++req_count]);
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
        const pointlist_msg_type& msg = i.second;

        const auto& szonerange = msg.szonerange;
        const auto& szone = msg.szone;
        const auto& jzonerange = msg.jzonerange;
        const auto& jzone = msg.jzone;
        const auto& ptList = msg.ptList;

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

            zone_to_zone_to_list_opp[jzid][szid].insert(zone_to_zone_to_list_opp[jzid][szid].end(), &ptList[n_beg], &ptList[n_beg] + n_end);
          }

        }
      }
      return 0;
    }


  };
}

#endif

/*    
    Copyright 2013-2026 ONERA.

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

#ifndef NUGA_MPI_MSG_TYPES_HXX
#define NUGA_MPI_MSG_TYPES_HXX

#include "Nuga/include/mpi_stl.hxx"

namespace NUGA
{
  enum eMPI_Tag { TAG_DATA_SZ=2, TAG_XRANGE_SZ=4, TAG_XRANGE=5, TAG_PGS_SZ=6, TAG_PGS=7, TAG_SZONE_SZ=8, TAG_SZONERANGE_SZ=10, TAG_JZONE_SZ=12, TAG_JZONERANGE_SZ=14, TAG_PTL_SZ=16, TAG_MPI_EXCH_STATUS=18, TAG_HASSENT=99 };

  inline
  static int get_opp_zone(const std::map< int, std::pair<int,int> >& rid_to_zones, int rid, int zid)
  {
    const auto itjz = rid_to_zones.find(rid); assert(itjz != rid_to_zones.end());
    int jzid = (itjz->second.first == zid) ? itjz->second.second : itjz->second.first;

    return jzid;
  }
  ///
  inline
  static void split_mpi_omp_joins
  (
    const std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list,
    int rank, 
    const std::map<int, std::pair<int, int>>& rid_to_zones,
    const std::vector<int>& zonerank,
    std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_omp,
    std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_mpi
   )
  {
    zone_to_rid_to_list_omp.clear();
    zone_to_rid_to_list_mpi.clear();

    for (auto & it : zone_to_rid_to_list)
    {
      int zid = it.first;
      auto & rid_to_list = it.second;
      for (auto & it2 : rid_to_list)
      {
        int rid = it2.first;
        int jzid = get_opp_zone(rid_to_zones, rid, zid);
        auto& ptlist = it2.second;
        int rk = zonerank[jzid];
        if (rk == rank) // local join
          zone_to_rid_to_list_omp[zid][rid] = ptlist;
        else // distant
          zone_to_rid_to_list_mpi[zid][rid] = ptlist;
      }
    }
  }
  
  // data to serialize for exchange : one DynArray per face

  template <typename T>
  struct plan_msg_type
  {
    std::vector<T> data; // data contains all construction plans contiguouslsy. e.g a construction plan as a IntArray 4x3 is stored in data as 12 ints. get_data_stride() function in receiver side will allow to figure out the 4x3 structure
    std::vector<T> datarange; // i-th plan is stored between datarange[i] and datarange[i+1] in data
    std::vector<T> pgs;  // gives  the positions in pointlist : i-th plan is intended for pointlist[pgs[i]]
    std::vector<T> szone, szonerange; // "sending zone" and associated range : all plans concerned by szone[i] are in data between szonerange[i] and szonerange[i+1]
    std::vector<T> rac, racrange; // rid and associated range : all plans concerned by rid = rac[i] are in data between racrange[i] and racrange[i+1]

    void clear() { data.clear(); datarange.clear(); pgs.clear(); szone.clear(); szonerange.clear(); rac.clear(); racrange.clear(); }

    ///
    static void convert_to_MPI_exchange_format
    (
      const std::map<int, std::map<int, std::map<E_Int, K_FLD::DynArray<T>>>> & sz_to_rid_to_PG_to_plan,
      const std::map<int, std::pair<int,int>> & rid_to_zones,
      const std::vector<int>& zonerank,
      std::map<int, plan_msg_type> & rank_to_data
    )
    {
      rank_to_data.clear();

      // Gather data by rank
      std::map<int, std::map<int, std::map<int, std::map<E_Int, K_FLD::DynArray<T>>>>> rank_to_sz_to_rid_to_PG_to_plan;
      for (auto & it : sz_to_rid_to_PG_to_plan)
      {
        int szid = it.first;
        const auto& rid_to_PG_to_plan = it.second;

        for (auto& it2 : rid_to_PG_to_plan)
        {

          int rid = it2.first;
          int jzid = get_opp_zone(rid_to_zones, rid, szid);

          auto& PG_to_plan = it2.second;
          assert(jzid > -1 && (size_t)jzid < zonerank.size());
          int jrk = zonerank[jzid];

          rank_to_sz_to_rid_to_PG_to_plan[jrk][szid][rid] = PG_to_plan;
        }
      }

      //
      for (auto& i : rank_to_sz_to_rid_to_PG_to_plan)
      {
        int rank = i.first;
        auto & sz_to_rid_to_PG_to_plan = i.second;

        auto & rankdata = rank_to_data[rank]; //creation et referencement

        for (auto & kk : sz_to_rid_to_PG_to_plan)
        {
          int szid = kk.first;
          auto & rid_to_PG_to_plan = kk.second;

          rankdata.szonerange.push_back(rankdata.rac.size());
          rankdata.szone.push_back(szid);

          for (auto& r : rid_to_PG_to_plan)
          {
            int rid = r.first;
            //int jzid = get_opp_zone(rid_to_zones, rid, szid);

            auto& PG_to_plan = r.second;

            rankdata.racrange.push_back(rankdata.pgs.size());
            rankdata.rac.push_back(rid);

            for (auto& k : PG_to_plan)
            {
              E_Int iloc/*PGi*/ = k.first;
              auto p = k.second; //copy to ensure contiguous

              assert(p.cols() != 0);
              //E_Int sz0 = rankdata.data.size();
              rankdata.datarange.push_back(rankdata.data.size());
              rankdata.data.insert(rankdata.data.end(), p.begin(), p.end());
              //assert(rankdata.data.size() > sz0);
              rankdata.pgs.push_back(iloc);
            }
          }
        }

        //one-pass-the-end
        rankdata.datarange.push_back(rankdata.data.size());
        rankdata.racrange.push_back(rankdata.pgs.size());
        rankdata.szonerange.push_back(rankdata.rac.size());
      }
    }

    ///
    static void isend_data
    (
      int rank, int nrank,
      MPI_Comm COM,
      const std::map<int, plan_msg_type> & rank_to_data,
      bool* has_sent,
      std::vector<MPI_Request>& sreqs
    )
    {
      //if (rank==3) std::cout << "send_data 1" << std::endl;
      
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

        NUGA::MPI::Isend(data.rac, rankid, TAG_JZONE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.racrange, rankid, TAG_JZONERANGE_SZ, COM, &sreqs[++count_req]);
      }

      // tell if sent or not to all ranks (to avoid deadlock when receiving)

      for (int r = 0; r < nrank; ++r) 
      {
        if (r == rank) continue;
        //std::cout << "rank : " << rank  << " over (" << nrank << ") has sent to " << r << " ? : " << has_sent[r] << std::endl;
        MPI_Isend(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreqs[++count_req]);
      }
    }

    ///
    static int receive_data
    (
      int rank, int nranks, MPI_Comm COM,
      int data_stride,
      const std::map<int, std::pair<int,int>> & rid_to_zones,
      const std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list,
      std::vector<MPI_Request>& sender_reqs,
      std::map<int, std::map<E_Int, K_FLD::DynArray<T>>>& zone_to_sensor_data
    )
    {
      zone_to_sensor_data.clear();

      STACK_ARRAY(bool, nranks, has_sent);
      has_sent[rank] = false;

      int NB_REQ_FOR_HASSENT = nranks - 1;
      {
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_HASSENT, sreqs);
        int req_count = -1;
        for (int r = 0; r < nranks; ++r)
        {
          if (r == rank) continue;
          MPI_Irecv(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreqs[++req_count]);
        }

        if (NB_REQ_FOR_HASSENT > 0)
        {
          STACK_ARRAY(MPI_Status, NB_REQ_FOR_HASSENT, status);
          MPI_Waitall(NB_REQ_FOR_HASSENT, sreqs.get(), status.get());
        }
      }

      // 1. catch all sizes :
      int NB_SENT = 0;
      for (int r = 0; r < nranks; ++r) if (has_sent[r])++NB_SENT;

      int NB_REQ_FOR_SIZES = 7 * NB_SENT;
      STACK_ARRAY(MPI_Request, NB_REQ_FOR_SIZES, sreqs_sz);
      // 1.1. data
      STACK_ARRAY(int, NB_SENT, data_sz);
      // 1.2. datarange
      STACK_ARRAY(int, NB_SENT, datarange_sz);
      // 1.3. pgs
      STACK_ARRAY(E_Int, NB_SENT, pgs_sz);
      // 1.4. szone
      STACK_ARRAY(int, NB_SENT, szone_sz);
      // 1.5. szonerange
      STACK_ARRAY(int, NB_SENT, szonerange_sz);
      // 1.4. szone
      STACK_ARRAY(int, NB_SENT, rac_sz);
      // 1.5. szonerange
      STACK_ARRAY(int, NB_SENT, racrange_sz);

      int req_count = -1;
      int rcount = -1;
      for (int r = 0; r < nranks; ++r)
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

        err = MPI_Irecv(&rac_sz[rcount], 1, MPI_INT, r, int(TAG_JZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&racrange_sz[rcount], 1, MPI_INT, r, int(TAG_JZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;
      }

      assert((req_count + 1) == NB_REQ_FOR_SIZES);
      if (NB_REQ_FOR_SIZES > 0)
      {
        //std::cout << "req count vs nranks vs NB_REQ_FOR_SIZES : " << req_count << "/" << nranks << "/" << NB_REQ_FOR_SIZES << std::endl;
        MPI_Waitall(NB_REQ_FOR_SIZES, sreqs_sz.get(), MPI_STATUS_IGNORE);
      }


      rcount=-1;
      for (int r = 0; r < nranks; ++r)
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
        assert(rac_sz[rcount] != 0);
        assert(racrange_sz[rcount] != 0);
      }

      // 2. get data
      std::map<int, plan_msg_type> rank_to_data;
      {
        int NB_REQ_FOR_DATA = 7 * NB_SENT;
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_DATA, sreqs_data);

        req_count = -1;
        rcount=-1;
        for (int r = 0; r < nranks; ++r)
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

          err = NUGA::MPI::Irecv(rac_sz[rcount], r, int(TAG_JZONE_SZ) + 1, COM, rank_to_data[r].rac, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(racrange_sz[rcount], r, int(TAG_JZONERANGE_SZ) + 1, COM, rank_to_data[r].racrange, &sreqs_data[++req_count]);
          if (err) return 1;
        }

        MPI_Waitall(sender_reqs.size(), &sender_reqs[0], MPI_STATUS_IGNORE);

        ++req_count;
        if (req_count > 0)
        {
          //std::cout << "req count vs nranks vs NB_REQ_FOR_SIZES : " << req_count << "/" << nranks << "/" << NB_REQ_FOR_SIZES << std::endl;
          MPI_Waitall(req_count, sreqs_data.get(), MPI_STATUS_IGNORE);
        }
      }


      // 3. dispatch sensor data by zone
      for (auto & i : rank_to_data)
      {
        //int rankid = i.first;
        const plan_msg_type& msg = i.second;

        const auto& szonerange = msg.szonerange;
        const auto& szone = msg.szone;
        const auto& racrange = msg.racrange;
        const auto& rac = msg.rac;
        const auto& datarange = msg.datarange;
        const auto& data = msg.data;
        const auto & pgs = msg.pgs;

        assert(pgs.size() == datarange.size() - 1);
        assert(szone.size() == szonerange.size() - 1);
        assert(rac.size() == racrange.size() - 1);

        for (size_t u = 0; u < szone.size(); ++u)
        {
          assert(u >= 0 && u < szone.size());
          int szid = szone[u];
          //if (rank == 1) std::cout << "szid : " << szid << std::endl;

          int racrange_beg = szonerange[u];
          int racrange_end = szonerange[u + 1];

          // if (rank == 2) std::cout << "racrange_beg : " << racrange_beg << std::endl;
          // if (rank == 2) std::cout << "racrange_end : " << racrange_end << std::endl;
          // if (rank == 2) std::cout << "receive_data : 10 c 3" << std::endl;

          for (E_Int v = racrange_beg; v < racrange_end; ++v)
          {
            assert(v >= 0 && (size_t)v < rac.size());
            int rid = rac[v];
            int jzid = get_opp_zone(rid_to_zones, rid, szid);

            E_Int n_beg = racrange[v];
            E_Int n_end = racrange[v + 1];

            // if (rank == 2) std::cout << "n_beg : " << n_beg << std::endl;
            // if (rank == 2) std::cout << "n_end : " << n_end << std::endl;

            const auto it = zone_to_rid_to_list.find(jzid);
            assert(it != zone_to_rid_to_list.end());
            const auto & rid_to_list = it->second;

            const auto& itPtList = rid_to_list.find(rid);
            assert(itPtList != rid_to_list.end());
            const auto& ptlist = itPtList->second;

            for (E_Int n = n_beg; n < n_end; ++n)
            {
              assert(n >= 0 && (size_t)n < datarange.size());
              assert((size_t)(n + 1) < datarange.size());
              E_Int data_beg = datarange[n];
              E_Int data_end = datarange[n + 1];
              assert((size_t)n < pgs.size());
              E_Int iloc = pgs[n];

              E_Int sz = data_end - data_beg;
              

              E_Int cols = sz / data_stride;
             
              assert(data_beg >= 0 && (size_t)data_beg < data.size());
              const T* p = &data[data_beg];
              K_FLD::DynArray<T> a;
              a.resize(data_stride, cols);
              
              for (E_Int c = 0; c < cols; ++c)
              {
                for (E_Int k = 0; k < data_stride; ++k)
                  a(k, c) = *(p++);
              }

              assert(iloc >= 0 && (size_t)iloc < ptlist.size());
              E_Int PGi = ptlist[iloc] - 1;
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

  // data to serialize for exchange : one pointlist per rac

  struct pointlist_msg_type
  {
    std::vector<E_Int> ptList, szone, szonerange, rac, racrange;

    void clear() { ptList.clear(); szone.clear(); szonerange.clear(); rac.clear(); racrange.clear(); }

    ///
    static void transpose_pointlists
    (
      const std::map<int, std::pair<int,int>>& rid_to_zones,
      const std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
      std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_opp
    )
    {
      zone_to_rid_to_list_opp.clear();

      // OMP : just transpose
      if (!zone_to_rid_to_list_owned.empty())
      {
        for (const auto& it : zone_to_rid_to_list_owned)
        {
          int zid = it.first;
          const auto& rid_to_list = it.second;
          for (const auto & z2L : rid_to_list)
          {
            int rid = z2L.first;
            int jzid = get_opp_zone(rid_to_zones, rid, zid);
            const auto & ptL = z2L.second;

            zone_to_rid_to_list_opp[jzid][rid] = ptL;
          }
        }
      }

      assert(zone_to_rid_to_list_owned.size() == zone_to_rid_to_list_opp.size());
    }

    ///
    static void exchange_pointlists
    (
      const std::map<int, std::pair<int,int>>& rid_to_zones,
      const std::vector<int>& zonerank,
      MPI_Comm COM,
      int rank, int nranks,
      const std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_owned,
      std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_opp
    )
    {
      zone_to_rid_to_list_opp.clear();

      //separate MPI/OMP joins
      std::map<int, std::map<int, std::vector<E_Int>>> zone_to_rid_to_list_mpi, zone_to_rid_to_list_omp;
      split_mpi_omp_joins(zone_to_rid_to_list_owned, rank, rid_to_zones, zonerank, zone_to_rid_to_list_omp, zone_to_rid_to_list_mpi);

      // OMP : just transpose
      transpose_pointlists(rid_to_zones, zone_to_rid_to_list_omp, zone_to_rid_to_list_opp);

      // MPI
      if (zone_to_rid_to_list_mpi.empty()) return;

      // prepare data to send : rank/zone/ptL
      std::map<int, pointlist_msg_type> rank_to_mpi_data;
      convert_to_MPI_exchange_format(zone_to_rid_to_list_mpi, rid_to_zones, zonerank, rank_to_mpi_data);

      // Send MPI data   
      isend_data(rank, nranks, COM, rank_to_mpi_data);

      // Receive MPI data and build sensor data by zone
      receive_data(rank, nranks, COM, rid_to_zones, zone_to_rid_to_list_opp);

      assert(zone_to_rid_to_list_owned.size() == zone_to_rid_to_list_opp.size());
    }

    
    ///
    static void convert_to_MPI_exchange_format
    (
      const std::map<int, std::map<int, std::vector<E_Int>>> & zone_to_rid_to_list_mpi,
      const std::map<int, std::pair<int,int>>& rid_to_zones,
      const std::vector<int>& zonerank,
      std::map<int, pointlist_msg_type>& rank_to_data
    )
    {
      rank_to_data.clear();

      // Gather data by rank
      std::map<int, std::map<int, std::map<E_Int, std::vector<E_Int>>>> rank_to_sz_to_rid_to_PTL;
      for (auto & it : zone_to_rid_to_list_mpi)
      {
        int szid = it.first;
        const auto& rid_to_PTL = it.second;

        for (auto& it2 : rid_to_PTL)
        {
          int rid = it2.first;
          auto& PTL = it2.second;

          assert(rid > -1 && (size_t)rid < rid_to_zones.size());
          int jzid = get_opp_zone(rid_to_zones, rid, szid);
          int jrk = zonerank[jzid];

          rank_to_sz_to_rid_to_PTL[jrk][szid][rid] = PTL;
          //std::cout << "prepa send : ...[jrk][szid][rid] : [" << jrk << "][" << szid << "][" << rid << "] : sz PTL : " << PTL.size() << std::endl; 
        }
      }

      //
      for (auto& i : rank_to_sz_to_rid_to_PTL)
      {
        int rank = i.first;
        auto & sz_to_rid_to_PTL = i.second;

        auto & rankdata = rank_to_data[rank];

        for (auto & kk : sz_to_rid_to_PTL)
        {
          int szid = kk.first;
          auto & rid_to_PTL = kk.second;

          rankdata.szonerange.push_back(rankdata.rac.size());
          rankdata.szone.push_back(szid);

          for (auto& r : rid_to_PTL)
          {
            int rid = r.first;
            auto& PTL = r.second;

            rankdata.racrange.push_back(rankdata.ptList.size());
            rankdata.rac.push_back(rid);

            rankdata.ptList.insert(rankdata.ptList.end(), ALL(PTL));
          }
        }

        //one-pass-the-end
        rankdata.racrange.push_back(rankdata.ptList.size());
        rankdata.szonerange.push_back(rankdata.rac.size());
      }
    }

    ///
    static void isend_data
    (
      int rank, int nrank,
      MPI_Comm COM,
      const std::map<int, pointlist_msg_type> & rank_to_data
    )
    {
      int nsranks = rank_to_data.size(); // nb of ranks to send to
      int NB_TOT_REQS = 5 * nsranks;    // for 5 vectors : ptL/szone/szonerange/rac/racrange. 

      STACK_ARRAY(MPI_Request, NB_TOT_REQS, sreqs);

      STACK_ARRAY(bool, nrank, has_sent);
      for (E_Int n = 0; n < nrank; ++n) has_sent[n] = false;
      
      int count_req{ -1 };
      for (auto& d : rank_to_data) // WARNING: reference is mandatory ! otherwise ISend might not finish before iteration d copy is deleted
      {
        int rankid = d.first;
        auto& data = d.second;

        has_sent[rankid] = true;

        NUGA::MPI::Isend(data.ptList, rankid, TAG_PTL_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szone, rankid, TAG_SZONE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.szonerange, rankid, TAG_SZONERANGE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.rac, rankid, TAG_JZONE_SZ, COM, &sreqs[++count_req]);

        NUGA::MPI::Isend(data.racrange, rankid, TAG_JZONERANGE_SZ, COM, &sreqs[++count_req]);
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
      for (size_t r = 0; r < (size_t)nrank; ++r) {
        //std::cout << "rank : " << r << " has sent ? : " << has_sent[r] << std::endl;
        if (r == (size_t)rank) continue;
        MPI_Isend(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreq[++count_req]);
      }

      STACK_ARRAY(MPI_Status, nrank - 1, stats);
      MPI_Waitall(nrank - 1, sreq.get(), stats.get());
    }

    ///
    static int receive_data
    (
      int rank, int nranks, MPI_Comm COM,
      const std::map<int, std::pair<int,int>>& rid_to_zones,
      std::map<int, std::map<int, std::vector<E_Int>>>& zone_to_rid_to_list_opp
    )
    {
      STACK_ARRAY(bool, nranks, has_sent);
      has_sent[rank] = false;

      int NB_REQ_FOR_HASSENT = nranks - 1;
      {
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_HASSENT, sreqs);
        int req_count = -1;
        for (int r = 0; r < nranks; ++r)
        {
          if (r == rank) continue;
          /*int err = */MPI_Irecv(&has_sent[r], 1, MPI_C_BOOL, r, int(TAG_HASSENT), COM, &sreqs[++req_count]);
        }

        if (NB_REQ_FOR_HASSENT > 0)
        {
          STACK_ARRAY(MPI_Status, NB_REQ_FOR_HASSENT, status);
          MPI_Waitall(NB_REQ_FOR_HASSENT, sreqs.get(), status.get());
        }
      }
      
      // 1. catch all sizes :
      int NB_SENT = 0;
      for (int r = 0; r < nranks; ++r) if (has_sent[r])++NB_SENT;

      int NB_REQ_FOR_SIZES = 5 * NB_SENT;
      STACK_ARRAY(MPI_Request, NB_REQ_FOR_SIZES, sreqs_sz);
      // 1.1. data
      STACK_ARRAY(E_Int, NB_SENT, ptList_sz);
      // 1.4. szone
      STACK_ARRAY(E_Int, NB_SENT, szone_sz);
      // 1.5. szonerange
      STACK_ARRAY(E_Int, NB_SENT, szonerange_sz);
      // 1.4. szone
      STACK_ARRAY(E_Int, NB_SENT, rac_sz);
      // 1.5. szonerange
      STACK_ARRAY(E_Int, NB_SENT, racrange_sz);

      int req_count = -1;
      int rcount = -1;
      for (int r = 0; r < nranks; ++r)
      {
        if (has_sent[r] == false) continue;
        ++rcount;

        int err = MPI_Irecv(&ptList_sz[rcount], 1, MPI_INT, r, int(TAG_PTL_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szone_sz[rcount], 1, MPI_INT, r, int(TAG_SZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&szonerange_sz[rcount], 1, MPI_INT, r, int(TAG_SZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&rac_sz[rcount], 1, MPI_INT, r, int(TAG_JZONE_SZ), COM, &sreqs_sz[++req_count]);
        if (err) return err;

        err = MPI_Irecv(&racrange_sz[rcount], 1, MPI_INT, r, int(TAG_JZONERANGE_SZ), COM, &sreqs_sz[++req_count]);
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
      for (size_t r = 0; r < (size_t)nranks; ++r)
      {
        if (has_sent[r] == false) continue;
        ++rcount;
        //std::cout << "rank : " << rank << " has_sent[" << r << "] : " << has_sent[r] << std::endl;
        //std::cout << "rank : " << rank << " datarange_sz[" << r << "] : " << datarange_sz[r] << std::endl;
        //std::cout << "rank : " << rank << " data_sz[" << r << "] : " << data_sz[r] << std::endl;

        assert(ptList_sz[rcount] != 0);
        assert(szone_sz[rcount] != 0);
        assert(szonerange_sz[rcount] != 0);
        assert(rac_sz[rcount] != 0);
        assert(racrange_sz[rcount] != 0);
      }

      // 2. get data
      std::map<int, pointlist_msg_type> rank_to_data;
      {
        int NB_REQ_FOR_DATA = 5 * NB_SENT;
        STACK_ARRAY(MPI_Request, NB_REQ_FOR_DATA, sreqs_data);

        req_count = -1;
        rcount = -1;
        for (size_t r = 0; r < (size_t)nranks; ++r)
        {
          if (has_sent[r] == false) continue;
          ++rcount;

          int err = NUGA::MPI::Irecv(ptList_sz[rcount], r, int(TAG_PTL_SZ) + 1, COM, rank_to_data[r].ptList, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szone_sz[rcount], r, int(TAG_SZONE_SZ) + 1, COM, rank_to_data[r].szone, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(szonerange_sz[rcount], r, int(TAG_SZONERANGE_SZ) + 1, COM, rank_to_data[r].szonerange, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(rac_sz[rcount], r, int(TAG_JZONE_SZ) + 1, COM, rank_to_data[r].rac, &sreqs_data[++req_count]);
          if (err) return 1;

          err = NUGA::MPI::Irecv(racrange_sz[rcount], r, int(TAG_JZONERANGE_SZ) + 1, COM, rank_to_data[r].racrange, &sreqs_data[++req_count]);
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
        const auto& racrange = msg.racrange;
        const auto& rac = msg.rac;
        const auto& ptList = msg.ptList;

        assert(szone.size() == szonerange.size() - 1);
        assert(rac.size() == racrange.size() - 1);

        for (size_t u = 0; u < szone.size(); ++u)
        {
          assert(u >= 0 && u < szone.size());
          int szid = szone[u];
          //if (rank == 1) std::cout << "szid : " << szid << std::endl;

          int racrange_beg = szonerange[u];
          int racrange_end = szonerange[u + 1];

          //if (rank == 1) std::cout << "racrange_beg : " << racrange_beg << std::endl;
          //if (rank == 1) std::cout << "racrange_end : " << racrange_end << std::endl;

          for (size_t v = racrange_beg; v < (size_t)racrange_end; ++v)
          {
            assert(v >= 0 && v < rac.size());
            int rid = rac[v];
            int jzid = get_opp_zone(rid_to_zones, rid, szid);

            int n_beg = racrange[v];
            int n_end = racrange[v + 1];
            int sz = n_end - n_beg;

            zone_to_rid_to_list_opp[jzid][rid].insert(zone_to_rid_to_list_opp[jzid][rid].end(), &ptList[n_beg], &ptList[n_beg] + sz);
          }
        }
      }
      return 0;
    }


  };
}

#endif

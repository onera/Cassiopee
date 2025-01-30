/*    
    Copyright 2013-2025 Onera.

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
#ifndef _XCORE_XMPI_COMMUNICATOR_IMPL_STUB_HPP_
#define _XCORE_XMPI_COMMUNICATOR_IMPL_STUB_HPP_
#include <algorithm>
#include <iostream>

#include "xmpi/communicator.hpp"
#include "xmpi/constantes.hpp"
#include "xmpi/status.hpp"

namespace xcore {
struct communicator::Implementation {
  Implementation() : m_pt_sendbuffer(nullptr) {}
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  Implementation(const Implementation& impl, int color, int key)
      : m_pt_sendbuffer(nullptr) {}
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  Implementation(const Implementation& impl) : m_pt_sendbuffer(nullptr) {}
  // .............................................................
  Implementation(const Communicator_ext_t& com) : m_pt_sendbuffer(nullptr) {}
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ~Implementation() {}
  // .............................................................
  int getRank() const { return 0; }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  int getSize() const { return 1; }
  // .............................................................
  static int comm_world;
  const Communicator_ext_t& get_ext_comm() const { comm_world = 0; return comm_world; }
  Communicator_ext_t& get_ext_comm() { return comm_world; }
  // .................................................................
  template <typename K>
  error send(const K& obj, int dest, int tag) const {
    m_pt_sendbuffer = (void*)&obj;
    return error{};
  }
  // .............................................................
  template <typename K>
  error send(std::size_t nbItems, const K* sndbuff, int dest, int tag) const {
    if (dest != 0) return error::rank;
    if (sndbuff == nullptr) return error::buffer;
    m_nbItems = nbItems;
    m_pt_sendbuffer = sndbuff;
    return error::success;
  }
  // .................................................................
  template <typename K>
  request isend(const K& obj, int dest, int tag) const {
    (K*)m_pt_sendbuffer = &obj;
    return request{};
  }
  // .................................................................
  template <typename K>
  request issend(const K& obj, int dest, int tag) const {
    (K*)m_pt_sendbuffer = &obj;
    return request{};
  }
  // .............................................................
  template <typename K>
  request isend(std::size_t nbItems, const K* sndbuff, int dest,
                int tag) const {
    m_nbItems = nbItems;
    m_pt_sendbuffer = sndbuff;
    return request();
  }
  // .............................................................
  template <typename K>
  request issend(std::size_t nbItems, const K* sndbuff, int dest,
                 int tag) const {
    m_nbItems = nbItems;
    m_pt_sendbuffer = (void*)sndbuff;
    return request();
  }
  // .............................................................
  template <typename K>
  status recv(K& rcvobj, int sender, int tag) const {
    rcvobj = *static_cast<const K*>(m_pt_sendbuffer);
    return status{};
  }
  // .............................................................
  template <typename K>
  status recv(std::size_t nbItems, K* rcvbuff, int sender, int tag) const {
    if (m_pt_sendbuffer == nullptr) {
      return status{0, 0, error::buffer};
    }
    status stat;
    if (m_nbItems < nbItems) {
      stat = status{int(m_nbItems), 0, error::count};
      nbItems = m_nbItems;
    } else
      stat = status{int(nbItems), tag, error::success};
    if (m_pt_sendbuffer != rcvbuff)
      std::copy_n(static_cast<const K*>(m_pt_sendbuffer), nbItems, rcvbuff);
    m_pt_sendbuffer = nullptr;
    return stat;
  }
  // .............................................................
  template <typename K>
  request irecv(K& rcvobj, int sender, int tag) const {
    rcvobj = *static_cast<const K*>(m_pt_sendbuffer);
    return request{};
  }
  // .............................................................
  template <typename K>
  request irecv(std::size_t nbItems, K* rcvbuff, int sender, int tag) const {
    status stat;
    if (m_nbItems < nbItems) {
      nbItems = m_nbItems;
    }
    if (m_pt_sendbuffer != rcvbuff)
      std::copy_n(static_cast<const K*>(m_pt_sendbuffer), nbItems, rcvbuff);
    m_pt_sendbuffer = nullptr;
    return request();
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  status probe(int src, int tag) const {
    return status{int(m_nbItems), tag, error::success};
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  template <typename K>
  void broadcast(const K* obj_snd, K& obj_rcv, int root) const {
    if (&obj_rcv != obj_snd) obj_rcv = *obj_snd;
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  template <typename K>
  void broadcast(std::size_t nbItems, const K* bufsnd, K* bufrcv,
                 int root) const {
    if (bufsnd != bufrcv) std::copy_n(bufsnd, nbItems, bufrcv);
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  template <typename K>
  request ibroadcast(std::size_t nbItems, const K* bufsnd, K* bufrcv,
                     int root) const {
    if (bufsnd != bufrcv) std::copy_n(bufsnd, nbItems, bufrcv);
    return request();
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  void barrier() const {}
  // .............................................................
  template <typename K>
  void reduce(std::size_t nbItems, const K* objs, K* res, Operation op,
              int root) {}
  // .............................................................
  template <typename K>
  void reduce(const K& objs, K* res, Operation op, int root) {}
  // .............................................................
  template <typename K>
  void reduce(const K& objs, K* res, Operation op, bool commute, int root) {}
  // .........................................................................................
  template <typename K>
  void allreduce(const K& loc, K* glob, const Operation& op) const {
    if (&loc != glob) *glob = loc;
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  template <typename K>
  void allreduce(std::size_t nbItems, const K* objs, K* res, Operation op) {
    if (res != objs) std::copy(objs, objs + nbItems, res);
  }
  // .........................................................................................
  template <typename K, typename F>
  void allreduce(const K& loc, K* glob, const F& fct, bool is_commuting) const {
    if (&loc != glob) *glob = loc;
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  template <typename K>
  void gather(std::size_t nbItems, const K* objs, K* recv_objs, int root) {
    if (objs != recv_objs) std::copy(objs, objs + nbItems, recv_objs);
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  template <typename K>
  void gather(const K& objs, K* recv_objs, int root) {
    recv_objs[0] = objs;
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  void translateRanks(communicator::Implementation& o_impl, int nbRanks,
                      const int* ranks, int* tr_ranks) const {
    tr_ranks[0] = ranks[0];
  }
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  bool iprobe(int source, int tag, status& st) const { return true; }
  // .............................................................
 private:
  mutable std::size_t m_nbItems;
  mutable void* m_pt_sendbuffer;
};
// -----------------------------------------------------------------
}  // namespace xcore
#endif

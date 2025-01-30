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
#ifndef _XCORE_XMPI_COMMUNICATOR_MPI_IMPL_HPP_
#define _XCORE_XMPI_COMMUNICATOR_MPI_IMPL_HPP_

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>

// Communicator for MPI
#include "xmpi/constantes.hpp"
#include "xmpi/status.hpp"
#include "xmpi/xtype_traits.hpp"

namespace xcore {
namespace {
#if __cplusplus >= 201402L
template <typename K>
std::function<K(const K &, const K &)> reduce_functor;
// ay <= accumulator
template <typename K>
void reduce_user_function(void *x, void *y, int *length, MPI_Datatype *tp) {
  K val;
  K *ax = (K *)x;
  K *ay = (K *)y;
  for (int i = 0; i < *length; ++i) {
    val = reduce_functor<K>(ax[i], ay[i]);
    ay[i] = val;
  }
}
#else
template <typename K>
void reduce_user_function(void *x, void *y, int *length, MPI_Datatype *tp) {
  assert(false && "Only implemented for C++ 14 or greater version");
}
#endif
}  // namespace
// #################################################################################################
struct communicator::Implementation {
  Implementation() { MPI_Comm_dup(MPI_COMM_WORLD, &m_communicator); }
  // ...............................................................................................
  Implementation(const Implementation &impl, int color, int key) {
    MPI_Comm_split(impl.m_communicator, color, key, &m_communicator);
  }
  // ...............................................................................................
  Implementation(const Implementation &impl) {
    MPI_Comm_dup(impl.m_communicator, &m_communicator);
  }
  // ...............................................................................................
  Implementation(const MPI_Comm &excom) {
    MPI_Comm_dup(excom, &m_communicator);
  }
  // ...............................................................................................
  ~Implementation() { MPI_Comm_free(&m_communicator); }
  // ...............................................................................................
  Implementation &operator=(const Implementation &com) = delete;
  // -----------------------------------------------------------------------------------------------
  int getRank() const {
    int rank;
    MPI_Comm_rank(m_communicator, &rank);
    return rank;
  }
  // ...............................................................................................
  void translateRanks(communicator::Implementation &o_impl, int nbRanks,
                      const int *ranks, int *tr_ranks) const {
    MPI_Group group1, group2;
    MPI_Comm_group(m_communicator, &group1);
    MPI_Comm_group(o_impl.m_communicator, &group2);
    int error =
        MPI_Group_translate_ranks(group1, nbRanks, ranks, group2, tr_ranks);
    if (error != MPI_SUCCESS)
      std::cerr << "Warning : error translating ranks.(" << error << ")"
                << std::endl;
    // On utilise cette routine aussi avec des rangs invalides qu'on gère !
    // assert( tr_ranks[ 0 ] >= 0 );
  }
  // ...............................................................................................
  int getSize() const {
    int size;
    MPI_Comm_size(m_communicator, &size);
    return size;
  }
  // ...............................................................................................
  const Communicator_ext_t &get_ext_comm() const { return m_communicator; }
  Communicator_ext_t &get_ext_comm() { return m_communicator; }
  // ===============================================================================================
  // For communication, we handle container in special functions. The flag
  // here tells if the object
  // is a container or not. See below for the specialization of this
  // structure for containers
  // ( specialization of the functions mainly ).
  template <typename K, bool flag = false>
  struct Communication {
    static void send(const MPI_Comm &com, const K &snd_obj, int dest, int tag) {
      if (Type_MPI<K>::must_be_packed()) {
        MPI_Send(&snd_obj, sizeof(K), MPI_BYTE, dest, tag, com);
      } else {
        MPI_Send(&snd_obj, 1, Type_MPI<K>::mpi_type(), dest, tag, com);
      }
    }
    // .......................................................................................
    static request isend(const MPI_Comm &com, const K &snd_obj, int dest,
                         int tag) {
      request req;
      if (Type_MPI<K>::must_be_packed()) {
        MPI_Isend(&snd_obj, sizeof(K), MPI_BYTE, dest, tag, com,
                  req.m_req.get());
      } else {
        MPI_Isend(&snd_obj, 1, Type_MPI<K>::mpi_type(), dest, tag, com,
                  req.m_req.get());
      }
      return req;
    }
    // .......................................................................................
    static request issend(const MPI_Comm &com, const K &snd_obj, int dest,
                          int tag) {
      request req;
      if (Type_MPI<K>::must_be_packed()) {
        MPI_Issend(&snd_obj, sizeof(K), MPI_BYTE, dest, tag, com,
                   req.m_req.get());
      } else {
        MPI_Issend(&snd_obj, 1, Type_MPI<K>::mpi_type(), dest, tag, com,
                   req.m_req.get());
      }
      return req;
    }
    // .......................................................................................
    static status recv(const MPI_Comm &com, K &rcvobj, int sender, int tag) {
      status status;
      if (Type_MPI<K>::must_be_packed()) {
        MPI_Recv(&rcvobj, sizeof(K), MPI_BYTE, sender, tag, com,
                 &status.mpi_status);
      } else {
        MPI_Recv(&rcvobj, 1, Type_MPI<K>::mpi_type(), sender, tag, com,
                 &status.mpi_status);
      }
      return status;
    }
    // .......................................................................................
    static request irecv(const MPI_Comm &com, K &rcvobj, int sender, int tag) {
      request req;
      if (Type_MPI<K>::must_be_packed()) {
        MPI_Irecv(&rcvobj, sizeof(K), MPI_BYTE, sender, tag, com,
                  req.m_req.get());
      } else {
        MPI_Irecv(&rcvobj, 1, Type_MPI<K>::mpi_type(), sender, tag, com,
                  req.m_req.get());
      }
      return req;
    }
    // .......................................................................................
    static void broadcast(const MPI_Comm &com, const K *obj_snd, K &obj_rcv,
                          int root) {
      int rank;
      MPI_Comm_rank(com, &rank);
      if (root == rank) {
        assert(obj_snd != nullptr);
        obj_rcv = *obj_snd;
      }
      if (Type_MPI<K>::must_be_packed()) {
        MPI_Bcast(&obj_rcv, sizeof(K), MPI_BYTE, root, com);
      } else {
        MPI_Bcast(&obj_rcv, 1, Type_MPI<K>::mpi_type(), root, com);
      }
    }
    // .......................................................................................
    static void reduce(const MPI_Comm &com, const K &loc, K *glob,
                       const Operation &op, int root) {
      MPI_Reduce(&loc, glob, 1, Type_MPI<K>::mpi_type(), op, root, com);
    }
    // .......................................................................................
    static void allreduce(const MPI_Comm &com, const K &loc, K *glob,
                          const Operation &op) {
      MPI_Allreduce(&loc, glob, 1, Type_MPI<K>::mpi_type(), op, com);
    }
    // .......................................................................................
    static void gather(const MPI_Comm &com, const K &loc, K *gather_array,
                       int root) {
      MPI_Gather(&loc, 1, Type_MPI<K>::mpi_type(), gather_array, 1,
                 Type_MPI<K>::mpi_type(), root, com);
    }
  };
  // ---------------------------------------------------------------------------------------------
  // General send if above methods don't work :
  template <typename K>
  void send(std::size_t nbItems, const K *sndbuff, int dest, int tag) const {
    if (Type_MPI<K>::must_be_packed()) {
      MPI_Send(sndbuff, nbItems * sizeof(K), MPI_BYTE, dest, tag,
               m_communicator);
    } else {
      MPI_Send(sndbuff, nbItems, Type_MPI<K>::mpi_type(), dest, tag,
               m_communicator);
    }
  }
  // ...........................................................................................
  template <typename K>
  void send(const K &snd, int dest, int tag) const {
    Communication<K, is_container<K>::value>::send(m_communicator, snd, dest,
                                                   tag);
  }
  // -------------------------------------------------------------------------------------------
  template <typename K>
  request isend(std::size_t nbItems, const K *sndbuff, int dest,
                int tag) const {
    request req;
    if (Type_MPI<K>::must_be_packed()) {
      MPI_Isend(sndbuff, nbItems * sizeof(K), MPI_BYTE, dest, tag,
                m_communicator, req.m_req.get());
    } else {
      MPI_Isend(sndbuff, nbItems, Type_MPI<K>::mpi_type(), dest, tag,
                m_communicator, req.m_req.get());
    }
    return req;
  }
  // .........................................................................................
  template <typename K>
  request isend(const K &snd, int dest, int tag) const {
    request req = Communication<K, is_container<K>::value>::isend(
        m_communicator, snd, dest, tag);
    return req;
  }
  // -------------------------------------------------------------------------------------------
  template <typename K>
  request issend(std::size_t nbItems, const K *sndbuff, int dest,
                 int tag) const {
    int ierr;
    request req;
    if (Type_MPI<K>::must_be_packed()) {
      ierr = MPI_Issend(sndbuff, nbItems * sizeof(K), MPI_BYTE, dest, tag,
                        m_communicator, req.m_req.get());
    } else {
      ierr = MPI_Issend(sndbuff, nbItems, Type_MPI<K>::mpi_type(), dest, tag,
                        m_communicator, req.m_req.get());
    }
    if (ierr != MPI_SUCCESS) {
      char errMsg[1024];
      int lenStr;
      MPI_Error_string(ierr, errMsg, &lenStr);
      std::cerr << __PRETTY_FUNCTION__ << " : Erreur " << errMsg << std::endl;
      MPI_Abort(MPI_COMM_WORLD, ierr);
      exit(EXIT_FAILURE);
    }
    return req;
  }
  // .........................................................................................
  template <typename K>
  request issend(const K &snd, int dest, int tag) const {
    request req = Communication<K, is_container<K>::value>::issend(
        m_communicator, snd, dest, tag);
    return req;
  }
  // -------------------------------------------------------------------------------------------
  // Default receive :
  template <typename K>
  status recv(std::size_t nbItems, K *rcvbuff, int sender, int tag) const {
    status status;
    if (Type_MPI<K>::must_be_packed()) {
      MPI_Recv(rcvbuff, nbItems * sizeof(K), MPI_BYTE, sender, tag,
               m_communicator, &status.mpi_status);
    } else {
      MPI_Recv(rcvbuff, nbItems, Type_MPI<K>::mpi_type(), sender, tag,
               m_communicator, &status.mpi_status);
    }
    return status;
  }
  // .........................................................................................
  template <typename K>
  status recv(K &rcvobj, int sender, int tag) const {
    status stat = Communication<K, is_container<K>::value>::recv(
        m_communicator, rcvobj, sender, tag);
    return stat;
  }
  // -------------------------------------------------------------------------------------------
  template <typename K>
  request irecv(std::size_t nbItems, K *rcvbuff, int sender, int tag) const {
    request req;
    if (Type_MPI<K>::must_be_packed()) {
      MPI_Irecv(rcvbuff, nbItems * sizeof(K), MPI_BYTE, sender, tag,
                m_communicator, req.m_req.get());
    } else {
      MPI_Irecv(rcvbuff, nbItems, Type_MPI<K>::mpi_type(), sender, tag,
                m_communicator, req.m_req.get());
    }
    return req;
  }
  // .........................................................................................
  template <typename K>
  request irecv(K &rcvobj, int sender, int tag) const {
    request req = Communication<K, is_container<K>::value>::irecv(
        m_communicator, rcvobj, sender, tag);
    return req;
  }
  // -----------------------------------------------------------------------------------------
  // Basic Broadcast :
  template <typename K>
  void broadcast(std::size_t nbItems, const K *bufsnd, K *bufrcv,
                 int root) const {
    MPI_Barrier(m_communicator);
    assert(bufrcv != nullptr);
    if (root == getRank()) {
      assert(bufsnd != nullptr);
      if (bufsnd != bufrcv) {
        std::copy_n(bufsnd, nbItems, bufrcv);
      }
    }
    if (Type_MPI<K>::must_be_packed()) {
      MPI_Bcast(bufrcv, nbItems * sizeof(K), MPI_BYTE, root, m_communicator);
    } else {
      MPI_Bcast(bufrcv, nbItems, Type_MPI<K>::mpi_type(), root, m_communicator);
    }
    MPI_Barrier(m_communicator);
  }
  // .........................................................................................
  template <typename K>
  void broadcast(const K *obj_snd, K &obj_rcv, int root) const {
    Communication<K, is_container<K>::value>::broadcast(m_communicator, obj_snd,
                                                        obj_rcv, root);
  }
  // -----------------------------------------------------------------------------------------
  template <typename K>
  request ibroadcast(std::size_t nbItems, const K *bufsnd, K *bufrcv,
                     int root) const {
    assert(bufrcv != nullptr);
    if (root == getRank()) {
      assert(bufsnd != nullptr);
      if (bufsnd != bufrcv) std::copy_n(bufsnd, nbItems, bufrcv);
    }
    request req;
    MPI_Ibcast(bufrcv, nbItems, Type_MPI<K>::mpi_type(), root, m_communicator,
               &req);
    return req;
  }
  // -----------------------------------------------------------------------------------------
  void barrier() const { MPI_Barrier(m_communicator); }
  // ----------------------------------------------------------------------------------------------------
  status probe(int source, int tag) const {
    status status;
    int ierr = MPI_Probe(source, tag, m_communicator, &status.mpi_status);
    if (ierr != MPI_SUCCESS) {
      char errMsg[1024];
      int lenStr;
      MPI_Error_string(ierr, errMsg, &lenStr);
      std::cerr << __PRETTY_FUNCTION__ << " : Erreur " << errMsg << std::endl;
      MPI_Abort(MPI_COMM_WORLD, ierr);
      exit(EXIT_FAILURE);
    }
    return status;
  }
  // ----------------------------------------------------------------------------------------------------
  bool iprobe(int source, int tag, status &st) const {
    int flag, ierr;
    ierr = MPI_Iprobe(source, tag, m_communicator, &flag, &st.mpi_status);
    if (ierr != MPI_SUCCESS) {
      char errMsg[1024];
      int lenStr;
      MPI_Error_string(ierr, errMsg, &lenStr);
      std::cerr << __PRETTY_FUNCTION__ << " : Erreur " << errMsg << std::endl;
      MPI_Abort(MPI_COMM_WORLD, ierr);
      exit(EXIT_FAILURE);
    }
    return (flag != 0);
  }
  // -----------------------------------------------------------------------------------------
  template <typename K>
  void reduce(std::size_t nbItems, const K *objs, K *res, Operation op,
              int root) {
    assert(objs != nullptr);
    if (root == getRank()) {
      assert(res != nullptr);
      if (objs == res) {
        MPI_Reduce(MPI_IN_PLACE, res, nbItems, Type_MPI<K>::mpi_type(), op,
                   root, m_communicator);
      } else {
        MPI_Reduce(objs, res, nbItems, Type_MPI<K>::mpi_type(), op, root,
                   m_communicator);
      }
    } else
      MPI_Reduce(objs, res, nbItems, Type_MPI<K>::mpi_type(), op, root,
                 m_communicator);
  }
  // .........................................................................................
  template <typename K>
  void reduce(const K &loc, K *glob, const Operation &op, int root) const {
    Communication<K, is_container<K>::value>::reduce(m_communicator, loc, glob,
                                                     op, root);
  }
  // -----------------------------------------------------------------------------------------
  template <typename K, typename F>
  void reduce(std::size_t nbItems, const K *objs, K *res, const F &fct,
              bool commute, int root) {
    assert(objs != nullptr);
#if __cplusplus >= 201402L
    reduce_functor<K> = std::function<K(const K &, const K &)>(fct);

    MPI_Op op;
    MPI_Op_create(reduce_user_function<K>, (commute ? 1 : 0), &op);

    if (root == getRank()) {
      assert(res != nullptr);
      if (objs == res) {
        MPI_Reduce(MPI_IN_PLACE, res, nbItems, Type_MPI<K>::mpi_type(), op,
                   root, m_communicator);
      } else {
        MPI_Reduce(objs, res, nbItems, Type_MPI<K>::mpi_type(), op, root,
                   m_communicator);
      }
    } else
      MPI_Reduce(objs, res, nbItems, Type_MPI<K>::mpi_type(), op, root,
                 m_communicator);
    MPI_Op_free(&op);
#else
    assert(false && "Only implemented for C++ 14 or greater");
#endif
  }
  // .........................................................................................
  template <typename K, typename F>
  void reduce(const K &loc, K *glob, const F &fct, bool is_commuting,
              int root) const {
#if __cplusplus >= 201402L
    reduce_functor<K> = std::function<K(const K &, const K &)>(fct);

    MPI_Op op;
    MPI_Op_create(reduce_user_function<K>, (is_commuting ? 1 : 0), &op);

    Communication<K, is_container<K>::value>::reduce(m_communicator, loc, glob,
                                                     op, root);

    MPI_Op_free(&op);
#else
    assert(
        false &&
        "Only available with a c++ compiler compliant with C++ 14 or greater");
#endif
  }
  // ====================================================================================================
  template <typename K>
  void allreduce(std::size_t nbItems, const K *objs, K *res, Operation op) {
    assert(objs != nullptr);
    if (objs == res) {
      MPI_Allreduce(MPI_IN_PLACE, res, nbItems, Type_MPI<K>::mpi_type(), op,
                    m_communicator);
    } else {
      MPI_Allreduce(objs, res, nbItems, Type_MPI<K>::mpi_type(), op,
                    m_communicator);
    }
  }
  // .........................................................................................
  template <typename K>
  void allreduce(const K &loc, K *glob, const Operation &op) const {
    Communication<K, is_container<K>::value>::allreduce(m_communicator, loc,
                                                        glob, op);
  }
  // -----------------------------------------------------------------------------------------
  template <typename K, typename F>
  void allreduce(std::size_t nbItems, const K *objs, K *res, const F &fct,
                 bool commute) {
    assert(objs != nullptr);
#if __cplusplus >= 201402L
    reduce_functor<K> = std::function<K(const K &, const K &)>(fct);

    MPI_Op op;
    MPI_Op_create(reduce_user_function<K>, (commute ? 1 : 0), &op);

    assert(res != nullptr);
    if (objs == res) {
      MPI_Allreduce(MPI_IN_PLACE, res, nbItems, Type_MPI<K>::mpi_type(), op,
                    m_communicator);
    } else {
      MPI_Allreduce(objs, res, nbItems, Type_MPI<K>::mpi_type(), op,
                    m_communicator);
    }
    MPI_Op_free(&op);
#else
    assert(
        false &&
        "Only available with a c++ compiler compliant with C++ 14 or greater");
#endif
  }
  // .........................................................................................
  template <typename K, typename F>
  void allreduce(const K &loc, K *glob, const F &fct, bool is_commuting) const {
#if __cplusplus >= 201402L
    reduce_functor<K> = std::function<K(const K &, const K &)>(fct);

    MPI_Op op;
    MPI_Op_create(reduce_user_function<K>, (is_commuting ? 1 : 0), &op);

    Communication<K, is_container<K>::value>::allreduce(m_communicator, loc,
                                                        glob, op);

    MPI_Op_free(&op);
#else
    assert(
        false &&
        "Only available with a c++ compiler compliant with C++ 14 or greater");
#endif
  }
  // ====================================================================================================
  template <typename K>
  void gather(std::size_t nbItems, const K *objs, K *res, int root) {
    assert(objs != nullptr);
    MPI_Gather(objs, nbItems, Type_MPI<K>::mpi_type(), res, nbItems,
               Type_MPI<K>::mpi_type(), root, m_communicator);
  }
  // .........................................................................................
  template <typename K>
  void gather(const K &loc, K *glob, int root) const {
    Communication<K, is_container<K>::value>::gather(m_communicator, loc, glob,
                                                     root);
  }

 private:
  MPI_Comm m_communicator;
};
// ###############################################################################################
// # Specialization of communication functions for containers :
template <typename K>
struct communicator::Implementation::Communication<K, true> {
  static void send(const MPI_Comm &com, const K &snd_arr, int dest, int tag) {
    std::vector<typename K::value_type, typename K::allocator_type> *snd;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      snd = (std::vector<typename K::value_type, typename K::allocator_type>
                 *)&snd_arr;
    } else {
      snd = new std::vector<typename K::value_type, typename K::allocator_type>(
          snd_arr.size());
      std::copy(snd_arr.begin(), snd_arr.end(), snd->begin());
    }

    if (Type_MPI<typename K::value_type>::must_be_packed()) {
      MPI_Send(snd->data(), snd->size() * sizeof(typename K::value_type),
               MPI_BYTE, dest, tag, com);
    } else {
      MPI_Send(snd->data(), snd->size(),
               Type_MPI<typename K::value_type>::mpi_type(), dest, tag, com);
    }
    if (!std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      delete snd;
    }
  }
  // .......................................................................................
  static request isend(const MPI_Comm &com, const K &snd_obj, int dest,
                       int tag) {
    std::vector<typename K::value_type, typename K::allocator_type> *snd;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      snd = (std::vector<typename K::value_type, typename K::allocator_type>
                 *)&snd_obj;
    } else {
      snd = new std::vector<typename K::value_type, typename K::allocator_type>(
          snd_obj.size());
      std::copy(snd_obj.begin(), snd_obj.end(), snd->begin());
    }
    request req;
    if (Type_MPI<typename K::value_type>::must_be_packed()) {
      MPI_Isend(snd->data(), snd->size() * sizeof(typename K::value_type),
                MPI_BYTE, dest, tag, com, req.m_req.get());
    } else {
      MPI_Isend(snd->data(), snd->size(),
                Type_MPI<typename K::value_type>::mpi_type(), dest, tag, com,
                req.m_req.get());
    }
    if (!std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      delete snd;
    }
    return req;
  }
  // .......................................................................................
  static request issend(const MPI_Comm &com, const K &snd_obj, int dest,
                        int tag) {
    std::vector<typename K::value_type, typename K::allocator_type> *snd;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      snd = (std::vector<typename K::value_type, typename K::allocator_type>
                 *)&snd_obj;
    } else {
      snd = new std::vector<typename K::value_type, typename K::allocator_type>(
          snd_obj.size());
      std::copy(snd_obj.begin(), snd_obj.end(), snd->begin());
    }

    request req;
    if (Type_MPI<typename K::value_type>::must_be_packed()) {
      MPI_Issend(snd->data(), snd->size() * sizeof(typename K::value_type),
                 MPI_BYTE, dest, tag, com, req.m_req.get());
    } else {
      MPI_Issend(snd->data(), snd->size(),
                 Type_MPI<typename K::value_type>::mpi_type(), dest, tag, com,
                 req.m_req.get());
    }
    if (!std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      delete snd;
    }
    return req;
  }
  // .......................................................................................
  static status recv(const MPI_Comm &com, K &rcvobj, int sender, int tag) {
    MPI_Status st;
    MPI_Probe(sender, tag, com, &st);
    int szMsg;

    MPI_Get_count(&st, Type_MPI<typename K::value_type>::mpi_type(), &szMsg);
    std::vector<typename K::value_type, typename K::allocator_type> *rcv;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      rcv = (std::vector<typename K::value_type, typename K::allocator_type>
                 *)&rcvobj;
      if (rcv->size() < (std::size_t)szMsg) {
        std::vector<typename K::value_type, typename K::allocator_type>(szMsg)
            .swap(*rcv);
      }
    } else {
      rcv = new std::vector<typename K::value_type, typename K::allocator_type>(
          szMsg);
    }
    if (Type_MPI<typename K::value_type>::must_be_packed()) {
      MPI_Recv(rcv->data(), rcv->size() * sizeof(typename K::value_type),
               MPI_BYTE, sender, tag, com, &st);
    } else {
      MPI_Recv(rcv->data(), rcv->size(),
               Type_MPI<typename K::value_type>::mpi_type(), sender, tag, com,
               &st);
    }
    if (!std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      rcvobj = K(rcv->begin(), rcv->end());
      delete rcv;
    }
    return status(st);
  }
  // .......................................................................................
  static request irecv(const MPI_Comm &com, K &rcvobj, int sender, int tag) {
    std::vector<typename K::value_type, typename K::allocator_type> *rcv;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      rcv = (std::vector<typename K::value_type, typename K::allocator_type>
                 *)&rcvobj;
    } else {
      rcv = new std::vector<typename K::value_type, typename K::allocator_type>(
          rcvobj.size());
    }

    request req;
    if (Type_MPI<typename K::value_type>::must_be_packed()) {
      MPI_Irecv(rcv->data(), rcv->size() * sizeof(typename K::value_type),
                MPI_BYTE, sender, tag, com, req.m_req.get());
    } else {
      MPI_Irecv(rcv->data(), rcv->size(),
                Type_MPI<typename K::value_type>::mpi_type(), sender, tag, com,
                req.m_req.get());
    }

    if (!std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      std::copy(rcvobj.begin(), rcvobj.end(), rcv->begin());
      delete rcv;
    }

    return req;
  }
  // .......................................................................................
  static void broadcast(const MPI_Comm &com, const K *obj_snd, K &obj_rcv,
                        int root) {
    std::size_t szMsg = (obj_snd != nullptr ? obj_snd->size() : obj_rcv.size());
    std::vector<typename K::value_type, typename K::allocator_type> *rcv;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      rcv = (std::vector<typename K::value_type, typename K::allocator_type>
                 *)&obj_rcv;
      if (szMsg > rcv->size()) {
        std::vector<typename K::value_type, typename K::allocator_type>(
            obj_snd->size())
            .swap(*rcv);
      }
    } else {
      rcv = new std::vector<typename K::value_type, typename K::allocator_type>(
          szMsg);
    }

    int rank;
    MPI_Comm_rank(com, &rank);
    if ((root == rank) && (&obj_rcv != obj_snd)) {
      assert(obj_snd != nullptr);
      std::copy(obj_snd->begin(), obj_snd->end(), obj_rcv.begin());
    }
    if (Type_MPI<typename K::value_type>::must_be_packed()) {
      MPI_Bcast(rcv->data(), rcv->size() * sizeof(typename K::value_type),
                MPI_BYTE, root, com);
    } else {
      MPI_Bcast(rcv->data(), rcv->size(),
                Type_MPI<typename K::value_type>::mpi_type(), root, com);
    }
    if (!std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      std::copy(obj_rcv.begin(), obj_rcv.end(), rcv->begin());
      delete rcv;
    }
  }
  // .......................................................................................
  static void reduce(const MPI_Comm &com, const K &loc, K *glob,
                     const Operation &op, int root) {
    std::size_t szMsg = loc.size();
    int rank;
    MPI_Comm_rank(com, &rank);
    std::vector<typename K::value_type, typename K::allocator_type> *glb;
    std::vector<typename K::value_type, typename K::allocator_type> *lc;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      glb = (std::vector<typename K::value_type, typename K::allocator_type> *)
          glob;
      lc = (std::vector<typename K::value_type, typename K::allocator_type>
                *)&loc;
      if (glb != nullptr) {
        if (szMsg > glb->size()) {
          std::vector<typename K::value_type, typename K::allocator_type>(szMsg)
              .swap(*glb);
        }
      }
    } else {
      if (rank == root)
        glb =
            new std::vector<typename K::value_type, typename K::allocator_type>(
                szMsg);
      else
        glb = nullptr;
      lc = new std::vector<typename K::value_type, typename K::allocator_type>(
          szMsg);
      std::copy(loc.begin(), loc.end(), lc->begin());
    }
    if (glb != nullptr)
      MPI_Reduce(lc->data(), glb->data(), szMsg,
                 Type_MPI<typename K::value_type>::mpi_type(), op, root, com);
    else
      MPI_Reduce(lc->data(), nullptr, szMsg,
                 Type_MPI<typename K::value_type>::mpi_type(), op, root, com);
    if (rank == root)
      if (!std::is_base_of<
              std::vector<typename K::value_type, typename K::allocator_type>,
              K>::value) {
        std::copy(glb->begin(), glb->end(), glob->begin());
        delete glb;
        delete lc;
      }
  }
  // ----------------------------------------------------------------------------------------------------
  static void allreduce(const MPI_Comm &com, const K &loc, K *glob,
                        const Operation &op) {
    std::size_t szMsg = loc.size();
    int rank;
    MPI_Comm_rank(com, &rank);
    std::vector<typename K::value_type, typename K::allocator_type> *glb;
    std::vector<typename K::value_type, typename K::allocator_type> *lc;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      glb = (std::vector<typename K::value_type, typename K::allocator_type> *)
          glob;
      lc = (std::vector<typename K::value_type, typename K::allocator_type>
                *)(&loc);
      if (glb != nullptr) {
        if (szMsg > glb->size()) {
          std::vector<typename K::value_type, typename K::allocator_type>(szMsg)
              .swap(*glb);
        }
      }
    } else {
      glb = new std::vector<typename K::value_type, typename K::allocator_type>(
          szMsg);
      lc = new std::vector<typename K::value_type, typename K::allocator_type>(
          szMsg);
      std::copy(loc.begin(), loc.end(), lc->begin());
    }
    MPI_Allreduce(lc->data(), glb->data(), szMsg,
                  Type_MPI<typename K::value_type>::mpi_type(), op, com);
    if (!std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      std::copy(glb->begin(), glb->end(), glob->begin());
      delete glb;
      delete lc;
    }
  }
  // ----------------------------------------------------------------------------------------------------
  template <typename C>
  static void gather(const MPI_Comm &com, const K &loc, K *glob, int root) {
    std::size_t szMsg = loc.size();
    int rank;
    MPI_Comm_rank(com, &rank);
    int sz;
    MPI_Comm_size(com, &sz);
    std::vector<typename K::value_type, typename K::allocator_type> *glb;
    std::vector<typename K::value_type, typename K::allocator_type> *lc;
    if (std::is_base_of<
            std::vector<typename K::value_type, typename K::allocator_type>,
            K>::value) {
      glb = (std::vector<typename K::value_type, typename K::allocator_type> *)
          glob;
      lc = (std::vector<typename K::value_type, typename K::allocator_type>
                *)&loc;
      if (glb != nullptr) {
        if (szMsg > glb->size()) {
          std::vector<typename K::value_type, typename K::allocator_type>(
              szMsg * sz)
              .swap(*glb);
        }
      }
    } else {
      if (rank == root)
        glb =
            new std::vector<typename K::value_type, typename K::allocator_type>(
                szMsg * sz);
      else
        glb = nullptr;
      lc = new std::vector<typename K::value_type, typename K::allocator_type>(
          szMsg);
      std::copy(loc.begin(), loc.end(), lc->begin());
    }
    if (glb != nullptr)
      MPI_Gather(lc->data(), szMsg,
                 Type_MPI<typename K::value_type>::mpi_type(), glb->data(),
                 szMsg, Type_MPI<typename K::value_type>::mpi_type(), root,
                 com);
    else
      MPI_Gather(lc->data(), szMsg,
                 Type_MPI<typename K::value_type>::mpi_type(), nullptr, szMsg,
                 Type_MPI<typename K::value_type>::mpi_type(), root, com);
    if (rank == root)
      if (!std::is_base_of<
              std::vector<typename K::value_type, typename K::allocator_type>,
              K>::value) {
        std::copy(glb->begin(), glb->end(), glob->begin());
        delete glb;
        delete lc;
      }
  }
};
}  // namespace xcore

#endif

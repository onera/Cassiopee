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

#ifndef _XCORE_XMPI_COMMUNICATOR_HPP_
#define _XCORE_XMPI_COMMUNICATOR_HPP_

#include <cstdlib>
#include <functional>
#include <memory>
#include <vector>

#include "xmpi/request.hpp"
#include "xmpi/status.hpp"

namespace xcore {
/*!   \class communicator
 *    \brief This class manages the data message exchanges ( point to point or
 * collective )
 *           inside a communication group.
 *
 *    The communicator class create instances which group some processes
 *    together. Two ways to create an instance of communicator class :
 *    1. With the default constructor : a global communicator instance is
 *       created which contains all processes executed in the parallel session
 *    2. With a color and a key : to create a partition of the processes in
 *       several communicators.
 *
 *    Probably than future versions of the library will provide other
 *    services to create new groups.
 *
 */
class communicator {
 public:
  /*!
   *   \brief Default constructor : instance a global communicator
   *
   *   The default constructor build an instance which contains all
   *   processes executed for the parallel session.
   */
  communicator();
  /*!
   *   \brief Split a communicator into a groupe of sub-communicators
   *
   *   This constructor instances a new communicator by splitting
   *   the communicator \ref com into a group of sub--communicators
   *   based on the input values \ref color and \ref key. It's important
   *   to note here that the original communicator doesn't go away,
   *   but a new communicator is created on each process.
   *
   *   \param com   It's the communicator that will be used as the basis
   *                for the new communicators
   *   \param color Determines to which new communicator each processes
   *                will belong. All processes which pass in the same
   *                value for \ref color are assigned to the same
   *                communicator. If the color is the constant undefined,
   *                that process won't be included in any of the new
   *                communicators.
   *   \param key   Determines the ordering ( rank ) within each
   *                new communicator. The process which passes in the
   *                smallest value for \ref key will be rank 0, the
   *                next smallest will be rank 1, and so on. If there is
   *                a tie, the process that had the lower rank in the original
   *                communicator will be first.
   *
   */
  communicator(const communicator &com, int color, int key);
  /*!*
   *  \brief Convert a communicator coming from external library used
   *         for Parallel library in Parallel communicator.
   *
   *  \param com  The external communicator coming from library used
   *              for implementation
   */
  communicator(const Communicator_ext_t &com);
  /*!
   *   \brief Duplicate a communicator in a new instance.
   *
   *   This constructor is used to create a new communicator that has a new
   *   communication context but contains the same group of processes
   *   as the input communicator.
   *
   *   \param com communicator to be duplicated
   */
  communicator(const communicator &com);

  communicator(communicator &&com) = delete;
  /*!
   *    Destructor. Destroy the communicator in the parallel context.
   */
  ~communicator();

  communicator &operator=(const communicator &com) = delete;
  communicator &operator=(communicator &&com) = delete;

  Communicator_ext_t &get_implementation();
  const Communicator_ext_t &get_implementation() const;

  // ===============================================================================================
  //                               Context of the communicator
  int rank; /*!< Rank of the current process inside the communicator instance */
  int size; /*!< Size of the communicator instance ( a.k.a number of processes
included in the communicator ) */

  /**
   * @brief      Gets the translated rank in other_com.
   *
   * @param[in]  other_com  An other communicator
   *
   * @return     The rank in the other_com communicator
   */
  int translateRank(const communicator &other_com) const;

  /**
   * @brief      Translate the rank rk for othercom communicator
   *
   * @param[in]  othercom  An other communicator
   * @param[in]  rk        The rank to translate
   *
   * @return     The translated rank
   */
  int translateRank(const communicator &othercom, int rk) const;

  /**
   * @brief      Translate an array of ranks
   *
   * @param[in]  othercom          An other communicator
   * @param[in]  ranksToTranslate  The ranks to translate
   *
   * @return     The translated ranks
   */
  std::vector<int> translateRanks(const communicator &othercom,
                                  const std::vector<int> &ranksToTranslate);
  // ===============================================================================================
  //                               Point to point communication
  /*!
   *    \brief Perform a blocking send to send an object to another process
   *
   *    This method send an object \ref obj to a process \ref dest with
   *    an indentifier \ref tag ( default value 0 ). This method performs
   *    a blocking send for large objects and a non blocking send for
   *    small objects.
   *
   *    NB : For a container, the send method send the data contained in the
   *         container... Serialiaze if not a vector.
   *
   *    \param obj  The object to send
   *    \param dest The rank of the destination
   *    \param tag  The message tag ( value default is zero )
   */
  template <typename K>
  void send(const K &obj, int dest, int tag = 0) const;
  /*!
   *    \brief Perform a blocking send for a buffer of objects to another
   * process
   *
   *    This method send a buffer \ref buff to a process \ref dest with
   *    an indentifier \ref tag ( default value 0 ). This method performs
   *    a blocking send for large array and a non blocking send for
   *    small array.
   *
   *    \param nbObjs Number of objects contained in the buffer
   *    \param buff   The buffer to send
   *    \param dest   The rank of the destination
   *    \param tag    The message tag ( value default is zero )
   */
  template <typename K>
  void send(std::size_t nbObjs, const K *buff, int dest, int tag = 0) const;
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /*!
   *    \brief Perform an asynchonous send on an object to another process
   *
   *    The method isend starts a standard-mode asynchronous send and return
   *    a request object which can be used later to query the status of the
   *    communication or wait for the completion of the send.
   *
   *    \param obj  The object to send
   *    \param dest The destination rank
   *    \param tag  The message tag
   *    \return     The request object associated at the send call
   */
  template <typename K>
  request isend(const K &obj, int dest, int tag = 0) const;
  /*!
   *    \brief Perform an asynchonous send on a buffer of objects to another
   * process
   *
   *    The method isend starts a standard-mode asynchronous send and return
   *    a request object which can be used later to query the status of the
   *    communication or wait for the completion of the send.
   *
   *    \param nbItems The number of items stored in the buffer
   *    \param obj     The buffer to send
   *    \param dest    The destination rank
   *    \param tag     The message tag
   *    \return        The request object associated at the send call
   */
  template <typename K>
  request isend(size_t nbItems, const K *obj, int dest, int tag = 0) const;
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /*!
   *    \brief Perform an asynchonous send on an object to another process
   *
   *    The method isend starts a standard-mode asynchronous send and return
   *    a request object which can be used later to query the status of the
   *    communication or wait for the completion of the send.
   *
   *    \param obj  The object to send
   *    \param dest The destination rank
   *    \param tag  The message tag
   *    \return     The request object associated at the send call
   */
  template <typename K>
  request issend(const K &obj, int dest, int tag = 0) const;
  /*!
   *    \brief Perform an asynchonous send on a buffer of objects to another
   * process
   *
   *    The method isend starts a standard-mode asynchronous send and return
   *    a request object which can be used later to query the status of the
   *    communication or wait for the completion of the send.
   *
   *    \param nbItems The number of items stored in the buffer
   *    \param obj     The buffer to send
   *    \param dest    The destination rank
   *    \param tag     The message tag
   *    \return        The request object associated at the send call
   */
  template <typename K>
  request issend(size_t nbItems, const K *obj, int dest, int tag = 0) const;
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /*!
   *    \brief Perform a blocking receive to receive an object sended by another
   * process
   *
   *    This method receive an object \ref obj from a process \ref sender
   *    with an indentifier \ref tag ( default value any_tag ). This method
   *    performs a blocking receive.
   *
   *    \param obj    The target object where one receive data of the sended
   * object
   *    \param sender The rank of the sender
   *    \param tag    The excepted message tag ( default value any_tag )
   *    \return Status The status of the received message.
   */
  template <typename K>
  status recv(K &obj, int sender = any_source, int tag = any_tag) const;
  /*!
   *    \brief Perform a blocking receive for a buffer of objects.
   *
   *    This method performs a receive Operation to receive a buffer of objects
   *    from \ref sender with ref \ref tag. This method performs a blocking
   * receive.
   *
   *    \param nbObjs  The number of objects to receive into the buffer.
   *    \param buff    The target buffer where store the receiving data.
   *    \param sender  The rank of the sender
   *    \param tag     The excepted message tag ( default value any_tag )
   *    \return Status The status of the received message.
   */
  template <typename K>
  status recv(std::size_t nbObjs, K *buff, int sender = any_source,
              int tag = any_tag) const;
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /*!
   *    \brief Perform a asynchronous receive for an object
   *
   *    This method performs a asynchronous receive operation. This method
   *    return before receiving the message and return a request. The request
   *    object allow to test if the message is received or not. The system
   *    begin to start writing data in the object after return the request
   *    object.
   *
   *    \param obj    The receive object
   *    \param sender Rank of the source
   *    \param tag    Message tag.
   *    \return       The request object associated at the receive message.
   */
  template <typename K>
  request irecv(K &obj, int sender = any_source, int tag = any_tag) const;
  /*!
   *    \brief Perform a asynchronous receive for a buffer of objects.
   *
   *    This method performs a asynchronous receive operation. This method
   *    return before receiving the message and return a request. The request
   *    object allow to test if the message is received or not. The system
   *    begin to start writing data in the buffer after return the request
   *    object.
   *
   *    \param nbItems Number of item to receive into the buffer
   *    \param obj     The receive vector
   *    \param sender  Rank of the source
   *    \param tag     Message tag.
   *    \return        The request object associated at the receive message.
   */
  template <typename K>
  request irecv(std::size_t nbItems, K *obj, int sender = any_source,
                int tag = any_tag) const;
  // ===============================================================================================
  //                                     Collective communication
  /*!
   *    \brief Perform a broadcast from a process to other processes.
   *
   *    This method performs a broadcast from the root process to
   *    other processes.
   *
   *    \param o_snd The object to broadcast.
   *    \param o_rcv The object where receive the broadcasted object.
   *    \param root  The rank of the root process
   */
  template <typename K>
  void bcast(const K &o_snd, K &o_rcv, int root = 0) const;
  /*!
   *    \brief Perform a broadcast from a process to other processes.
   *
   *    This method performs a broadcast from the root process to
   *    other processes. Don't call this method with the root process !
   *
   *    \param o_rcv The object where receive the broadcasted object.
   *    \param root  The rank of the root process
   */
  template <typename K>
  void bcast(K &o_rcv, int root = 0) const;
  /*!
   *    \brief Perform a broadcast from a process to other processes.
   *
   *    This method performs a broadcast from the root process to
   *    other processes.
   *
   *    \param Number of items to broadcast.
   *    \param b_snd The buffer of objects to broadcast.
   *    \param b_rcv The buffer of objects where receive broadcasted objects (
   * must be allocated
   *                 before the call of this method )
   *    \param root  The rank of the root process
   */
  template <typename K>
  void bcast(std::size_t nbObjs, const K *b_snd, K *b_rcv, int root = 0) const;
  /*!
   *    \brief Perform a broadcast from a process to other processes.
   *
   *    This method performs a broadcast from the root process to
   *    other processes. Don't call this method with the root process !
   *
   *    \param Number of items to broadcast.
   *    \param b_rcv The buffer of objects where receive broadcasted objects (
   * must be allocated
   *                 before the call of this method )
   *    \param root  The rank of the root process
   */
  template <typename K>
  void bcast(std::size_t nbObjs, K *b_rcv, int root = 0) const;
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /*!
   *    \brief Blocks until all processor inside the communicator have reached
   * this routine.
   *
   *    The processus wait until all processes contained in the current
   * communicator reached
   *    this routine.
   *
   */
  void barrier() const;
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members
   * of the
   *   current communicator. The reduction Operation must be here one of
   * predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process.
   *
   *   \param obj   An object used in the reduction operation
   *   \param res   The result object computed in the root process ( significant
   * only on root process )
   *   \param op    The pre-defined operation to do in the reduction operation
   *   \param root  The rank of the process where store the result of the
   * reduction operation
   */
  template <typename K>
  void reduce(const K &obj, K &res, const Operation &op, int root = 0) const;
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members
   * of the
   *   current communicator. The reduction Operation must be here one of
   * predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process. This version of this method must be
   * used only on non root process !
   *
   *   \param obj   An object used in the reduction operation
   *   \param op    The pre-defined operation to do in the reduction operation
   *   \param root  The rank of the process where store the result of the
   * reduction operation
   */
  template <typename K>
  void reduce(const K &obj, const Operation &op, int root = 0) const;
  // ----------------------------------------------------------------------------------
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members
   * of the
   *   current communicator. The reduction Operation must be here one of
   * predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process.
   *
   *   \param nbItems The number of items stored in the local buffer.
   *   \param obj     A buffer used in the reduction operation
   *   \param res     The result buffer computed in the root process (
   * significant only on root process )
   *   \param op      The pre-defined operation to do in the reduction operation
   *   \param root    The rank of the process where store the result of the
   * reduction operation
   */
  template <typename K>
  void reduce(std::size_t nbObjs, const K *b_objs, K *b_res, Operation op,
              int root = 0) const;
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members
   * of the
   *   current communicator. The reduction Operation must be here one of
   * predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process. This version of this method must be
   * used only on non root process !
   *
   *   \param nbItems The number of items stored in the local buffer.
   *   \param obj     A buffer used in the reduction operation
   *   \param op      The pre-defined operation to do in the reduction operation
   *   \param root    The rank of the process where store the result of the
   * reduction operation
   */
  template <typename K>
  void reduce(std::size_t nbObjs, const K *b_objs, Operation op,
              int root = 0) const;
  // ...............................................................................................
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members
   * of the
   *   current communicator. The reduction Operation must be here one of
   * predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process.
   *
   *   \param obj           An object used in the reduction operation
   *   \param res           The result object computed in the root process (
   * significant only on root process )
   *   \param op            An functor on two objects
   *   \param is_commutable A boolean to say if the function as commutable
   * parameters ( a.k.a op(x,y) = op(y,x) )
   *   \param root          The rank of the process where store the result of
   * the reduction operation
   */
  template <typename K, typename Func>
  void reduce(const K &obj, K &res, const Func &op, bool is_commutable = false,
              int root = 0) const;
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members
   * of the
   *   current communicator. The reduction Operation must be here one of
   * predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process. This version of this method must be
   * used only on non root process !
   *
   *   \param obj           An object used in the reduction operation
   *   \param op            An functor on two objects
   *   \param is_commutable A boolean to say if the function as commutable
   * parameters ( a.k.a op(x,y) = op(y,x) )
   *   \param root          The rank of the process where store the result of
   * the reduction operation
   */
  template <typename K, typename Func>
  void reduce(const K &obj, const Func &op, bool commute = false,
              int root = 0) const;
  // ...............................................................................................
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members of the current communicator. The
   * reduction Operation must be here one of predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process.
   *
   *   \param nbItems       The number of the objects in the buffer
   *   \param obj           A buffer used in the reduction operation
   *   \param res           The result buffer computed in the root process (
   * significant only on root process )
   *   \param op            An functor on two objects
   *   \param is_commutable A boolean to say if the function as commutable
   * parameters ( a.k.a op(x,y) = op(y,x) )
   *   \param root          The rank of the process where store the result of
   * the reduction operation
   */
  template <typename K, typename Func>
  void reduce(std::size_t nbObjs, const K *b_objs, K *b_res, const Func &op,
              bool is_commutable, int root = 0) const;
  /*!
   *   \brief Reduce values on all processes within current communicator
   *
   *   This method performs a global reduce Operation ( such as sum, max,
   * logical AND, etc. ) across all members
   * of the
   *   current communicator. The reduction Operation must be here one of
   * predefined list of Operations.
   *
   *   This method combine the element provided by each process in the
   * communicator, using the Operation op, and
   * returns the
   *   combined value for the root process. This version of this method must be
   * used only on non root process !
   *
   *   \param nbItems       The number of the objects in the buffer
   *   \param obj           A buffer used in the reduction operation
   *   \param op            An functor on two objects
   *   \param is_commutable A boolean to say if the function as commutable
   * parameters ( a.k.a op(x,y) = op(y,x) )
   *   \param root          The rank of the process where store the result of
   * the reduction operation
   */
  template <typename K, typename Func>
  void reduce(std::size_t nbObjs, const K *b_objs, const Func &op,
              bool is_commutable, int root = 0) const;
  // ====================================================================================================
  template <typename K>
  void allreduce(const K &obj, K &res, const Operation &op) const;
  template <typename K>
  void allreduce(std::size_t nbObjs, const K *b_objs, K *b_res,
                 Operation op) const;
  template <typename K, typename Func>
  void allreduce(const K &obj, K &res, const Func &op,
                 bool is_commutable = false) const;
  template <typename K, typename Func>
  void allreduce(std::size_t nbObjs, const K *b_objs, K *b_res, const Func &op,
                 bool is_commutable = false) const;
  // ===================================================================
  template <typename K>
  void gather(const K &obj, K *arr_of_objs = nullptr, int root = 0) const;
  template <typename K>
  void gather(std::size_t nb_objs_to_send, const K *objs,
              K *arr_of_revc_objs = nullptr, int root = 0);
  // ===================================================================
  status probe(int source = any_source, int tag = any_tag) const;
  // Return status with  if none message with specified source and tag is
  // available
  // else return the matching status.
  bool iprobe(status &status, int source = any_source, int tag = any_tag) const;

 private:
  struct Implementation;
  Implementation *m_impl;
};
}  // namespace xcore
#include "xmpi/communicator_impl.hpp"
#endif

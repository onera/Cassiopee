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

#ifndef _XCORE_XMPI_COMMUNICATOR_IMPL_HPP_
#define _XCORE_XMPI_COMMUNICATOR_IMPL_HPP_

#include <cassert>

#include "xmpi/communicator.hpp"
#if defined( _MPI )
#  include "xmpi/communicator_mpi_impl.hpp"
#else
#  include "xmpi/communicator_stub_impl.hpp"
#endif

namespace xcore {
    template <typename K>
    void communicator::send( const K& obj, int dest, int tag ) const {
        m_impl->send( obj, dest, tag );
    }
    // .................................................................
    template <typename K>
    void communicator::send( std::size_t nbObjs, const K* buff, int dest, int tag ) const {
        m_impl->send( nbObjs, buff, dest, tag );
    }
    // .................................................................
    template <typename K>
    request communicator::isend( const K& obj, int dest, int tag ) const {
        return m_impl->isend( obj, dest, tag );
    }
    // .................................................................
    template <typename K>
    request communicator::isend( std::size_t nbItems, const K* obj, int dest, int tag ) const {
        return m_impl->isend( nbItems, obj, dest, tag );
    }
    // .................................................................
    template <typename K>
    request communicator::issend( const K& obj, int dest, int tag ) const {
        return m_impl->issend( obj, dest, tag );
    }
    // .................................................................
    template <typename K>
    request communicator::issend( std::size_t nbItems, const K* obj, int dest, int tag ) const {
        return m_impl->issend( nbItems, obj, dest, tag );
    }
    // .................................................................
    template <typename K>
    status communicator::recv( K& obj, int sender, int tag ) const {
        return m_impl->recv( obj, sender, tag );
    }
    // .................................................................
    template <typename K>
    status communicator::recv( std::size_t nbObjs, K* buff, int sender, int tag ) const {
        return m_impl->recv( nbObjs, buff, sender, tag );
    }
    // .................................................................
    template <typename K>
    request communicator::irecv( K& obj, int sender, int tag ) const {
        return m_impl->irecv( obj, sender, tag );
    }
    // .................................................................
    template <typename K>
    request communicator::irecv( std::size_t nbObjs, K* buff, int sender, int tag ) const {
        return m_impl->irecv( nbObjs, buff, sender, tag );
    }
    // =================================================================
    // Opérations collectives :
    template <typename K>
    void communicator::bcast( const K& objsnd, K& objrcv, int root ) const {
        m_impl->broadcast( &objsnd, objrcv, root );
    }
    // .................................................................
    template <typename K>
    void communicator::bcast( K& objrcv, int root ) const {
        m_impl->broadcast( static_cast<K*>( nullptr ), objrcv, root );
    }
    // .................................................................
    template <typename K>
    void communicator::bcast( std::size_t nbObjs, const K* b_snd, K* b_rcv, int root ) const {
        m_impl->broadcast( nbObjs, b_snd, b_rcv, root );
    }
    // .................................................................
    template <typename K>
    void communicator::bcast( std::size_t nbObjs, K* b_rcv, int root ) const {
        m_impl->broadcast( nbObjs, (const K*)nullptr, b_rcv, root );
    }
    // =================================================================
    template <typename K>
    void communicator::reduce( const K& obj, K& res, const Operation& op, int root ) const {
        m_impl->reduce( obj, &res, op, root );
    }
    // .................................................................
    template <typename K>
    void communicator::reduce( const K& obj, const Operation& op, int root ) const {
        assert( root != rank );
        m_impl->reduce( 1, &obj, nullptr, op, root );
    }
    // _________________________________________________________________
    template <typename K, typename Func>
    void communicator::reduce( const K& obj, K& res, const Func& op, bool commute, int root ) const {
        m_impl->reduce( obj, &res, op, commute, root );
    }
    // .................................................................
    template <typename K, typename Func>
    void communicator::reduce( const K& obj, const Func& op, bool commute, int root ) const {
        assert( root != rank );
        m_impl->reduce( obj, nullptr, op, commute, root );
    }
    // _________________________________________________________________
    template <typename K>
    void communicator::reduce( std::size_t nbItems, const K* obj, K* res, Operation op, int root ) const {
        m_impl->reduce( nbItems, obj, res, op, root );
    }
    // .................................................................
    template <typename K>
    void communicator::reduce( std::size_t nbItems, const K* obj, Operation op, int root ) const {
        assert( rank != root );
        m_impl->reduce( nbItems, obj, nullptr, op, root );
    }
    // _________________________________________________________________
    template <typename K, typename Func>
    void communicator::reduce( std::size_t nbItems, const K* obj, K* res, const Func& op, bool commute,
                               int root ) const {
        m_impl->reduce( nbItems, obj, res, op, commute, root );
    }
    // .................................................................
    template <typename K, typename Func>
    void communicator::reduce( std::size_t nbItems, const K* obj, const Func& op, bool commute, int root ) const {
        assert( rank != root );
        m_impl->reduce( nbItems, obj, nullptr, op, commute, root );
    }
    // =================================================================
    template <typename K>
    void communicator::allreduce( const K& obj, K& res, const Operation& op ) const {
        m_impl->allreduce( obj, &res, op );
    }
    // _________________________________________________________________
    template <typename K, typename Func>
    void communicator::allreduce( const K& obj, K& res, const Func& op, bool commute ) const {
        m_impl->allreduce( obj, &res, op, commute );
    }
    // _________________________________________________________________
    template <typename K>
    void communicator::allreduce( std::size_t nbItems, const K* obj, K* res, Operation op ) const {
        m_impl->allreduce( nbItems, obj, res, op );
    }
    // _________________________________________________________________
    template <typename K, typename Func>
    void communicator::allreduce( std::size_t nbItems, const K* obj, K* res, const Func& op, bool commute ) const {
        m_impl->allreduce( nbItems, obj, res, op, commute );
    }
    // .................................................................
    template<typename K> 
    void communicator::gather( const K &obj, K* arr_of_objs, int root ) const
    {
        m_impl->gather( obj, arr_of_objs, root );
    }
    // _________________________________________________________________
    template<typename K> 
    void communicator::gather( std::size_t nb_objs_to_send, const K* objs, K* arr_of_revc_objs, int root )
    {
        m_impl->gather( nb_objs_to_send, objs, arr_of_revc_objs, root );
    }
}
#endif

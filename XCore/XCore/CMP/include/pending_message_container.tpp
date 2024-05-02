/*    
    Copyright 2013-2020 Onera.

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
// pending_message_container.tpp
///////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "pending_message_container.hpp"

namespace CMP {
    template<typename MsgBuffer>
    template<typename... Args> void 
    PendingMsgContainer<MsgBuffer>:: emplace_back(Args const & ...args)
    {
        m_pending_message.emplace_back(std::make_shared<MsgBuffer>(args...) );
    }

    template <typename MsgBuffer>
    void PendingMsgContainer<MsgBuffer>::push_back( pointer_message_buffer& pt_buff ) {
        m_pending_message.push_back( PendingMsgContainer::pending_message( pt_buff ) );
    }
    // -----------------------------------------------------------------------------------------------------
    template <typename MsgBuffer>
    void PendingMsgContainer<MsgBuffer>::pop( const iterator& it ) {
        ( *it ) = back( );
        m_pending_message.pop_back( );
    }
    // -----------------------------------------------------------------------------------------------------
    template <typename MsgBuffer>
    void PendingMsgContainer<MsgBuffer>::waitAll( ) {
        for ( iterator it = begin( ); it != end( ); ++it ) ( *it ).message_buffer->wait( );
    }
    // -----------------------------------------------------------------------------------------------------
    // If iterator is equal to end(), not complete message found
    template <typename MsgBuffer>
    typename PendingMsgContainer<MsgBuffer>::iterator
    PendingMsgContainer<MsgBuffer>::get_first_complete_message( ) {
        for ( iterator it = begin( ); it != end( ); ++it )
            if ( ( *it ).message_buffer->test( ) ) return it;
        return end( );
    }
    // -----------------------------------------------------------------------------------------------------
    template <typename MsgBuffer>
    typename PendingMsgContainer<MsgBuffer>::pending_message& PendingMsgContainer<MsgBuffer>::
    operator[]( std::size_t i ) {
        assert( i < size( ) );
        return m_pending_message[i];
    }
    // -----------------------------------------------------------------------------------------------------
    template <typename MsgBuffer>
    const typename PendingMsgContainer<MsgBuffer>::pending_message&
        PendingMsgContainer<MsgBuffer>::operator[]( std::size_t i ) const {
        assert( i < size( ) );
        return m_pending_message[i];
    }
}

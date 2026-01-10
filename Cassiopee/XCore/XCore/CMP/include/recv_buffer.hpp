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
// RecvBuffer.hpp
#ifndef _CMP_RECVBUFFER_HPP_
#define _CMP_RECVBUFFER_HPP_
#include <cassert>
#include <string>
#include <memory>
using std::shared_ptr;

#include "xmpi/xmpi.hpp"
#include "vector_view.hpp"

namespace CMP {
    class RecvBuffer {
    public:
        using value_t = unsigned char;
        RecvBuffer( int source = xcore::any_source, int id_tag = xcore::any_tag );
        RecvBuffer( int source, int id_tag, const xcore::communicator& comm );
        RecvBuffer( const RecvBuffer& r_buf );
        ~RecvBuffer( );

        template <typename K>
        RecvBuffer& operator>>( K& val ) {
            std::size_t sz;
            const value_t* pt_data = unpack( sz );
            assert( sz == sizeof( K ) );
            val = *(const K*)( pt_data );
            return *this;
        }
        template <typename K>
        RecvBuffer& operator>>( vector_view<K>& array ) {
            std::size_t sz;
            const value_t* pt_data = unpack( sz );
            std::size_t sz_arr  = sz / sizeof( K );
            assert( sz_arr * sizeof( K ) == sz );
            array = vector_view<K>( (const K*)( pt_data ), sz_arr );
            return *this;
        }
        RecvBuffer& operator>>( std::string& str ) {
            std::size_t sz;
            const value_t* pt_data = unpack( sz );
            str                 = std::string( (char*)pt_data, sz );
            return *this;
        }

        int  irecv( );
        bool test( xcore::status* pt_status = NULL );
        void wait( xcore::status* pt_status = NULL );
        int         source( ) const;
        int         tag( ) const;
        int         size( ) const;
        const value_t* data( ) const;
        const value_t* current_data( ) const;

    private:
        const value_t* unpack( std::size_t& sz );
        class Implementation;
        shared_ptr<Implementation> m_pt_implementation;
    };
}
#endif

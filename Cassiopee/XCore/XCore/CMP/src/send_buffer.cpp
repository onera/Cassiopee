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
// SendBuffer.cpp
#include <iostream>
#include <vector>
#include <cassert>

#include "CMP/include/send_buffer.hpp"
#include "xmpi/context.hpp"
//#define CMP_DEBUG

namespace CMP {
    // ===========================================================================================
    // Déclaration de la mise en oeuvre du buffer d'envoi
    class SendBuffer::Implementation {
    public:
        typedef std::vector<SendBuffer::PackedData> Datas;
        typedef std::vector<value_t>          Buffer;

        Implementation( int dest, int tag, const xcore::communicator& comm );
        Implementation( const Implementation& impl ) = delete;
        ~Implementation( ) {}

        Implementation& operator = ( const Implementation& ) = delete;

        void finalize();
        int  isend( );
        bool test( xcore::status* pt_status );
        int  wait( xcore::status* pt_status );

        int         receiver( ) const { return m_recv_rank; }
        int         tag( ) const { return m_id_tag; }
        std::size_t size( ) const { return m_cur_size; }
        
              value_t* data( )       { return m_arr_buffer.data(); }
        const value_t* data( ) const { return m_arr_buffer.data(); }

         SendBuffer::PackedData& pack( const SendBuffer::PackedData& data ) {
#           if defined(CMP_DEBUG) 
            auto& logg = xcore::context::logger();
            logg << LogTrace << std::endl;
#           endif
            m_cur_size += sizeof( std::size_t ) + data.size( );
            m_arr_pkg_data.push_back( data );
            return m_arr_pkg_data.back();
        }

      void clear() { m_cur_size = 0; m_arr_pkg_data.clear(); }

    private:
        void           copyPackagedDataInBuffer( );
        int            m_recv_rank;
        int            m_id_tag;
        const xcore::communicator& m_ref_comm;
        Datas          m_arr_pkg_data;
        Buffer         m_arr_buffer;
        std::size_t    m_cur_size;
        xcore::request m_request;
        bool           m_is_data_copied;
    };
    // ===========================================================================================
    // Définition de la mise en oeuvre du buffer d'envoi
    SendBuffer::Implementation::Implementation( int dest, int tag, const xcore::communicator& comm )
        : m_recv_rank( dest ),
          m_id_tag( tag ),
          m_ref_comm( comm ),
          m_arr_pkg_data( ),
          m_arr_buffer( ),
          m_cur_size( 0 ),
          m_request(),
          m_is_data_copied(false) {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        m_arr_pkg_data.reserve( 16384);//1024 );
    }
    // -------------------------------------------------------------------------------------------
    void SendBuffer::Implementation::finalize()
    {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        if (not m_is_data_copied) {
            copyPackagedDataInBuffer( );
            m_is_data_copied = true;
        }
    }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::Implementation::isend( ) {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        if (not m_is_data_copied) {
            copyPackagedDataInBuffer( );
            m_is_data_copied = true;
        }
        int         length  = int( m_arr_buffer.size( ) );
#       if defined(CMP_DEBUG) 
        logg << LogInformation << "Longueur message : " << length << std::endl;
#       endif
        m_request = m_ref_comm.issend(length, data(), m_recv_rank, m_id_tag);
        m_is_data_copied = false;
        return xcore::success;
    }
    // -------------------------------------------------------------------------------------------
    bool SendBuffer::Implementation::test( xcore::status* pt_status ) {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
       if (m_arr_buffer.size() == 0) return true;

        bool flag = m_request.test();
        if ( pt_status != NULL ) *pt_status = m_request.get_status();
        return flag;
    }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::Implementation::wait( xcore::status* pt_status ) {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        int ierr = xcore::success;
        if (m_arr_buffer.size() == 0) return ierr;
        m_request.wait();
#       if defined(CMP_DEBUG) 
        logg << LogInformation << m_ref_comm.rank <<  " : wait finish" << std::endl;
#       endif
        if ( pt_status != NULL ) *pt_status = m_request.get_status();
#       if defined(CMP_DEBUG) 
        logg << LogInformation << m_ref_comm.rank << " : Returning the error" << std::endl;
#       endif
        if ( pt_status != NULL ) return pt_status->error();
#       if defined(CMP_DEBUG) 
        logg << LogInformation << m_ref_comm.rank << " : Returning the error..." << std::endl;
#       endif
        return ierr;
        //return m_request.get_status().error();
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void SendBuffer::Implementation::copyPackagedDataInBuffer( ) {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        const std::size_t data_chunk                 = 2048;            // Copie par paquet de 2048 octets maximum
        const std::size_t min_size_for_parallel_copy = 8 * data_chunk;  // Taille minimal pour faire une copie parallèle
        if (m_cur_size > m_arr_buffer.size())
            std::vector<value_t>( m_cur_size ).swap( m_arr_buffer );
        std::size_t count = 0;
        //m_arr_buffer.resize( m_cur_size );
        Buffer::iterator itB = m_arr_buffer.begin( );

        for ( Datas::iterator it_d = m_arr_pkg_data.begin( ); it_d != m_arr_pkg_data.end( ); ++it_d ) {

            assert( itB < m_arr_buffer.end( ) );
            value_t* pt           = &( *itB );
            *(std::size_t*)( pt ) = ( *it_d ).size( );
            itB += sizeof( std::size_t );
            if (not (*it_d).must_copy()) {
                (*it_d).set_data(&(*itB));// A implementer.
            } else
            if ( ( *it_d ).size( ) < min_size_for_parallel_copy ) {
                std::copy( ( *it_d ).begin( ), ( *it_d ).end( ), itB );
		} else {
                std::size_t nb_chunks = ( *it_d ).size( ) / data_chunk;

#pragma omp                 parallel for
                for ( std::size_t i = 0; i < nb_chunks; ++i ) {


                    std::copy( ( *it_d ).begin( ) + i * data_chunk, ( *it_d ).begin( ) + ( i + 1 ) * data_chunk,
                               itB + i * data_chunk );
		}
                std::size_t size_copied = data_chunk * nb_chunks;

                std::copy( ( *it_d ).begin( ) + size_copied, ( *it_d ).end( ), itB + size_copied );
		}
            itB += ( *it_d ).size( );
            assert( itB <= m_arr_buffer.end( ) );
            count += ( *it_d ).size( );
        }
    }
    // ===========================================================================================
    // Définition du buffer d'envoi
    SendBuffer::SendBuffer( int recv_rank, int id_tag )
        : m_pt_implementation( ) 
    {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        xcore::communicator& comm = xcore::context::globalCommunicator();
        m_pt_implementation = std::make_shared<SendBuffer::Implementation>( recv_rank, id_tag, comm );
    }
    SendBuffer::SendBuffer( int recv_rank, int id_tag, const xcore::communicator& comm )
        : m_pt_implementation( new SendBuffer::Implementation( recv_rank, id_tag, comm ) ) 
    {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif        
    }
    // -------------------------------------------------------------------------------------------
    SendBuffer::SendBuffer( const SendBuffer& s_buf ) : m_pt_implementation( s_buf.m_pt_implementation ) {}
    // -------------------------------------------------------------------------------------------
    SendBuffer::~SendBuffer( ) {}
    // -------------------------------------------------------------------------------------------
    SendBuffer& SendBuffer::operator=( const SendBuffer& s_buf ) {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        if ( this != &s_buf ) { m_pt_implementation = s_buf.m_pt_implementation; }
        return *this;
    }
    // -------------------------------------------------------------------------------------------
    SendBuffer& SendBuffer::operator<<( const PackedData& data ) {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        m_pt_implementation->pack( data );
        return *this;
    }
    // -------------------------------------------------------------------------------------------
    SendBuffer::PackedData&
    SendBuffer::push_inplace_array( size_t size )
    {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        return m_pt_implementation->pack( SendBuffer::PackedData(size) );
    }
    // -------------------------------------------------------------------------------------------
    void SendBuffer::finalize_and_copy()
    {
#       if defined(CMP_DEBUG) 
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        return m_pt_implementation->finalize();
    }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::isend( ) { return m_pt_implementation->isend( ); }
    // -------------------------------------------------------------------------------------------
    bool SendBuffer::test( xcore::status* pt_status ) { return m_pt_implementation->test( pt_status ); }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::wait( xcore::status* pt_status ) { return m_pt_implementation->wait( pt_status ); }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::receiver( ) const { return m_pt_implementation->receiver( ); }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::tag( ) const { return m_pt_implementation->tag( ); }
    // -------------------------------------------------------------------------------------------
    std::size_t SendBuffer::size( ) const { return m_pt_implementation->size( ); }
    // -------------------------------------------------------------------------------------------
    auto SendBuffer::data( ) const -> const value_t* { return m_pt_implementation->data( ); }
    // -------------------------------------------------------------------------------------------
    void SendBuffer::clear() { m_pt_implementation->clear(); }
}

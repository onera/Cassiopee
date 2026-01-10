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
// recv_buffer.cpp
#include <iostream>
#include <cassert>

//#define CMP_DEBUG
#include "CMP/include/recv_buffer.hpp"
#include "xmpi/context.hpp"
namespace CMP {
    // Déclaration de la realisation du buffer de reception
    // ================================================================================================
    class RecvBuffer::Implementation {
    public:
        using value_t = typename RecvBuffer::value_t;
        typedef std::vector<value_t> buffer;
        typedef const value_t*       pointer_data;

        Implementation( int sender, int tag, const xcore::communicator& comm )
            : m_sender( sender ),
              m_id_tag( tag ),
              m_current_size( 0 ),
              m_rcv_buffer( ),
              m_cur_iterator(),
              m_request( ),
              m_ref_comm( comm )
            {
                m_cur_iterator = m_rcv_buffer.begin();
            }
        ~Implementation( ) {}

        int  irecv( );
        bool test( xcore::status* pt_status );
        void wait( xcore::status* pt_status );
        int         sender( ) const { return m_sender; }
        int         tag( ) const { return m_id_tag; }
        int         size( ) const { return m_rcv_buffer.size( ); }
        int current_size( ) const { return m_current_size; }
              value_t* data( )       { return &m_rcv_buffer[0]; }
        const value_t* data( ) const { return &m_rcv_buffer[0]; }
        const value_t* current_data( ) const { return &( *m_cur_iterator ); }

        const value_t* unpack( std::size_t& sz );

    private:
        void                              testAndGetLengthOfMessage( );
        void                              waitAndGetLengthOfMessage( );
        int                               m_sender;
        int                               m_id_tag;
        int                               m_current_size;
        std::vector<value_t>        m_rcv_buffer;
        std::vector<value_t>::const_iterator m_cur_iterator;

        xcore::request                    m_request ;
        const xcore::communicator&        m_ref_comm;
    };
    // Definition de la realisation du buffer de reception
    // ================================================================================================
    int RecvBuffer::Implementation::irecv( ) {
#       if defined(CMP_DEBUG)
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#endif
        int ierr = xcore::success;
        if (m_current_size == 0) testAndGetLengthOfMessage( );
#       if defined(CMP_DEBUG)    
        logg << LogInformation << "size to receive : " << m_current_size << std::endl;
#       endif
        if ( m_current_size > 0 ) {
            m_request = m_ref_comm.irecv(size(), data(), m_sender, m_id_tag);
        }
        return ierr;
    }
    // ------------------------------------------------------------------------------------------------
    bool RecvBuffer::Implementation::test( xcore::status* pt_status ) {
#       if defined(CMP_DEBUG)        
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#endif
        if ( m_current_size == 0 ) {
            testAndGetLengthOfMessage( );
#           if defined(CMP_DEBUG)        
            logg << LogInformation << "size to receive : " << size() << std::endl;
#           endif
            if ( m_current_size > 0 ) {
                m_request = m_ref_comm.irecv(size(), data(), m_sender, m_id_tag);
            } else
                return false;
        }
        bool flag = m_request.test();
        if ( pt_status != NULL )
            *pt_status = m_request.get_status();
        return flag;
    }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::Implementation::wait( xcore::status* pt_status ) {
#       if defined(CMP_DEBUG)        
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#endif
        if ( m_current_size == 0 ) {
            waitAndGetLengthOfMessage( );
            xcore::status loc_status = m_ref_comm.recv(size(), data(), m_sender, m_id_tag);
            if ( pt_status != NULL ) *pt_status = loc_status;
        } else {
#           if defined(CMP_DEBUG)        
            logg << LogInformation << "size : " << size() << std::endl;
#           endif            
            m_request.wait();            
            if ( pt_status != NULL ) *pt_status = m_request.get_status();
        }
        m_current_size = 0;
    }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::Implementation::testAndGetLengthOfMessage( ) 
    {
#       if defined(CMP_DEBUG)        
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#endif
        xcore::status status;
        bool is_available = m_ref_comm.iprobe(status, m_sender, m_id_tag);
#       if defined(CMP_DEBUG)        
        logg << LogInformation << "is_available : " << is_available << std::endl;
#       endif        
        if ( is_available ) 
        {
            m_current_size = status.count<unsigned char>();
#           if defined(CMP_DEBUG)
            logg << LogInformation << "Longueur message a recevoir : " 
                 << m_current_size << std::endl;
#           endif
            int m_rcv_buffer_size = m_rcv_buffer.size();
            if (m_current_size > m_rcv_buffer_size)
                std::vector<value_t>( m_current_size ).swap( m_rcv_buffer );
            m_cur_iterator = m_rcv_buffer.begin( );
        }
    }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::Implementation::waitAndGetLengthOfMessage( ) {
#       if defined(CMP_DEBUG)        
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#endif
        xcore::status status = m_ref_comm.probe(m_sender, m_id_tag);
        m_current_size = status.count<value_t>();
#       if defined(CMP_DEBUG)        
        logg << LogInformation << "Longueur message a recevoir... : " << m_current_size << std::endl;
#       endif        
        int m_rcv_buffer_size = m_rcv_buffer.size();
        if (m_current_size > m_rcv_buffer_size)
            std::vector<value_t>( m_current_size ).swap( m_rcv_buffer );
        m_cur_iterator = m_rcv_buffer.begin( );
    }
    // ------------------------------------------------------------------------------------------------
    const unsigned char* RecvBuffer::Implementation::unpack( std::size_t& sz ) {
#       if defined(CMP_DEBUG)        
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
        logg << LogInformation 
             << "m_cur_iterator distance : " << m_cur_iterator - m_rcv_buffer.begin()
             << std::endl;
#endif
        m_current_size = 0;
        assert(size()>0);
        assert(m_cur_iterator >= m_rcv_buffer.begin() );
        assert( m_cur_iterator + sizeof( sz ) < m_rcv_buffer.end( ) );
        const value_t* pt_cur = current_data( );
        sz                 = *( (const std::size_t*)( pt_cur ) );
        m_cur_iterator += sizeof( std::size_t );
        assert( m_cur_iterator + sz <= m_rcv_buffer.end( ) );
        const value_t* data = current_data( );
        m_cur_iterator += sz;
#       if defined(CMP_DEBUG)        
        logg << LogInformation << "Existing unpack function"  << std::endl;
#       endif       
        return data;
    }
    // Définition du buffer de réception
    // ================================================================================================
    RecvBuffer::RecvBuffer( int source, int id_tag )
        : m_pt_implementation( ) 
    {
#       if defined(CMP_DEBUG)
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
        xcore::communicator& comm = xcore::context::globalCommunicator();
        m_pt_implementation = std::make_shared<RecvBuffer::Implementation>( source, id_tag, comm );
    }
    RecvBuffer::RecvBuffer( int source, int id_tag, const xcore::communicator& comm )
        : m_pt_implementation( std::make_shared<RecvBuffer::Implementation>( source, id_tag, comm ) ) 
    {
#       if defined(CMP_DEBUG)
        auto& logg = xcore::context::logger();
        logg << LogTrace << std::endl;
#       endif
    }
    // ------------------------------------------------------------------------------------------------
    RecvBuffer::RecvBuffer( const RecvBuffer& r_buf ) : m_pt_implementation( r_buf.m_pt_implementation ) {}
    // ------------------------------------------------------------------------------------------------
    RecvBuffer::~RecvBuffer( ) {}
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::irecv( ) { return m_pt_implementation->irecv( ); }
    // ------------------------------------------------------------------------------------------------
    bool RecvBuffer::test( xcore::status* pt_status ) { return m_pt_implementation->test( pt_status ); }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::wait( xcore::status* pt_status ) { m_pt_implementation->wait( pt_status ); }
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::source( ) const { return m_pt_implementation->sender( ); }
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::tag( ) const { return m_pt_implementation->tag( ); }
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::size( ) const { return m_pt_implementation->size( ); }
    // ------------------------------------------------------------------------------------------------
    auto RecvBuffer::data( ) const -> const value_t* { return m_pt_implementation->data( ); }
    // ------------------------------------------------------------------------------------------------
    auto RecvBuffer::current_data( ) const -> const value_t*  { return m_pt_implementation->current_data( ); }
    // ------------------------------------------------------------------------------------------------
    auto RecvBuffer::unpack( std::size_t& sz ) -> const value_t*  { return m_pt_implementation->unpack( sz ); }
    // ------------------------------------------------------------------------------------------------
    }

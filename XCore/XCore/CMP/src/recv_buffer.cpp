// recv_buffer.cpp
#include "recv_buffer.hpp"
#include <iostream>
#include <cassert>

namespace CMP {
#if defined(_MPI)
    // Déclaration de la réalisation du buffer de réception
    // ================================================================================================
    class RecvBuffer::Implementation {
    public:
        typedef std::vector<char> buffer;
        typedef const char*       pointer_data;

        Implementation( int sender, int tag, const MPI_Comm& comm )
            : m_sender( sender ),
              m_id_tag( tag ),
              m_rcv_buffer( ),
              m_cur_iterator( ),
              m_request( ),
              m_ref_comm( comm ) {}
        ~Implementation( ) {}

        int  irecv( );
        bool test( MPI_Status* pt_status );
        void wait( MPI_Status* pt_status );
        int         sender( ) const { return m_sender; }
        int         tag( ) const { return m_id_tag; }
        int         size( ) const { return m_rcv_buffer.size( ); }
        const char* data( ) const { return &m_rcv_buffer[0]; }
        const char* current_data( ) const { return &( *m_cur_iterator ); }

        const char* unpack( std::size_t& sz );

    private:
        void                              testAndGetLengthOfMessage( );
        void                              waitAndGetLengthOfMessage( );
        int                               m_sender;
        int                               m_id_tag;
        std::vector<char>                 m_rcv_buffer;
        std::vector<char>::const_iterator m_cur_iterator;
        MPI_Request                       m_request;
        MPI_Comm                          m_ref_comm;
    };
    // Définition de la réalisation du buffer de réception
    // ================================================================================================
    int RecvBuffer::Implementation::irecv( ) {
      //std::cout << __PRETTY_FUNCTION__ << std::endl;
        int ierr = MPI_SUCCESS;
        if ( size( ) == 0 ) testAndGetLengthOfMessage( );
        if ( size( ) > 0 ) {
            m_cur_iterator = m_rcv_buffer.begin( );
            ierr = MPI_Irecv( &m_rcv_buffer[0], m_rcv_buffer.size( ), MPI_BYTE, m_sender, m_id_tag, m_ref_comm,
                              &m_request );
        }
        return ierr;
    }
    // ------------------------------------------------------------------------------------------------
    bool RecvBuffer::Implementation::test( MPI_Status* pt_status ) {
      //std::cout << __PRETTY_FUNCTION__ << std::endl;
        if ( size( ) == 0 ) {
            testAndGetLengthOfMessage( );
            if ( size( ) > 0 ) {
                MPI_Irecv( &m_rcv_buffer[0], m_rcv_buffer.size( ), MPI_BYTE, m_sender, m_id_tag, m_ref_comm,
                           &m_request );
            } else
                return false;
        }
        int flag;
        if ( pt_status == NULL )
            MPI_Test( &m_request, &flag, MPI_STATUS_IGNORE );
        else
            MPI_Test( &m_request, &flag, pt_status );
        return flag != 0;
    }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::Implementation::wait( MPI_Status* pt_status ) {
      //std::cout << __PRETTY_FUNCTION__ << std::endl;
        if ( size( ) == 0 ) {
            waitAndGetLengthOfMessage( );
            MPI_Status loc_status;
            if ( pt_status == NULL ) pt_status = &loc_status;
            MPI_Recv( &m_rcv_buffer[0], m_rcv_buffer.size( ), MPI_BYTE, m_sender, m_id_tag, m_ref_comm, pt_status );
        } else {
            if ( pt_status == NULL )
                MPI_Wait( &m_request, MPI_STATUS_IGNORE );
            else
                MPI_Wait( &m_request, pt_status );
        }
    }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::Implementation::testAndGetLengthOfMessage( ) {
      //std::cout << __PRETTY_FUNCTION__ << std::endl;
        MPI_Status status;
        int        is_available;
        MPI_Iprobe( m_sender, m_id_tag, m_ref_comm, &is_available, &status );
        if ( is_available ) {
            int length;
            MPI_Get_count( &status, MPI_BYTE, &length );
            //std::cout << "Recv buffer : " << length << std::endl;
            if (length > m_rcv_buffer.size() )
                std::vector<char>( length ).swap( m_rcv_buffer );
            m_cur_iterator = m_rcv_buffer.begin( );
        }
    }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::Implementation::waitAndGetLengthOfMessage( ) {
      //std::cout << __PRETTY_FUNCTION__ << std::endl;
        MPI_Status status;
        MPI_Probe( m_sender, m_id_tag, m_ref_comm, &status );
        int length;
        MPI_Get_count( &status, MPI_BYTE, &length );
        //std::cout << "Recv buffer : " << length << std::endl;
        if (length > m_rcv_buffer.size())
            std::vector<char>( length ).swap( m_rcv_buffer );
        m_cur_iterator = m_rcv_buffer.begin( );
    }
    // ------------------------------------------------------------------------------------------------
    const char* RecvBuffer::Implementation::unpack( std::size_t& sz ) {
      /*std::cout << __PRETTY_FUNCTION__ << std::endl;
      std::cout << "distance fin buffer et pointeur courant : " << m_rcv_buffer.end() - m_cur_iterator
                << std::endl;
                std::cout << "Taille donne a rajoutee : " << sizeof(sz) << " + ";*/
        assert( m_cur_iterator + sizeof( sz ) < m_rcv_buffer.end( ) );
        const char* pt_cur = current_data( );
        sz                 = *( (const std::size_t*)( pt_cur ) );
        //std::cout << sz << std::endl;
        m_cur_iterator += sizeof( std::size_t );
        assert( m_cur_iterator + sz <= m_rcv_buffer.end( ) );
        const char* data = current_data( );
        m_cur_iterator += sz;
        return data;
    }
#else
#endif
    // Définition du buffer de réception
    // ================================================================================================
    RecvBuffer::RecvBuffer( int source, int id_tag, const MPI_Comm& comm )
        : m_pt_implementation( new RecvBuffer::Implementation( source, id_tag, comm ) ) {}
    // ------------------------------------------------------------------------------------------------
    RecvBuffer::RecvBuffer( const RecvBuffer& r_buf ) : m_pt_implementation( r_buf.m_pt_implementation ) {}
    // ------------------------------------------------------------------------------------------------
    RecvBuffer::~RecvBuffer( ) {}
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::irecv( ) { return m_pt_implementation->irecv( ); }
    // ------------------------------------------------------------------------------------------------
    bool RecvBuffer::test( MPI_Status* pt_status ) { return m_pt_implementation->test( pt_status ); }
    // ------------------------------------------------------------------------------------------------
    void RecvBuffer::wait( MPI_Status* pt_status ) { m_pt_implementation->wait( pt_status ); }
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::source( ) const { return m_pt_implementation->sender( ); }
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::tag( ) const { return m_pt_implementation->tag( ); }
    // ------------------------------------------------------------------------------------------------
    int RecvBuffer::size( ) const { return m_pt_implementation->size( ); }
    // ------------------------------------------------------------------------------------------------
    const char* RecvBuffer::data( ) const { return m_pt_implementation->data( ); }
    // ------------------------------------------------------------------------------------------------
    const char* RecvBuffer::current_data( ) const { return m_pt_implementation->current_data( ); }
    // ------------------------------------------------------------------------------------------------
    const char* RecvBuffer::unpack( std::size_t& sz ) { return m_pt_implementation->unpack( sz ); }
    // ------------------------------------------------------------------------------------------------
    }

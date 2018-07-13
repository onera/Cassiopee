// SendBuffer.cpp
#include "send_buffer.hpp"
#include <iostream>
#include <vector>
#include <cassert>

namespace CMP {
    // ===========================================================================================
    // Déclaration de la mise en oeuvre du buffer d'envoi
    class SendBuffer::Implementation {
    public:
        typedef std::vector<SendBuffer::PackedData> Datas;
        typedef std::vector<char>                   Buffer;

        Implementation( int dest, int tag, const MPI_Comm& comm );
        Implementation( const Implementation& impl ) : m_ref_comm( impl.m_ref_comm ) {
            exit( -1 );
        }
        ~Implementation( ) {}

        void finalize();
        int  isend( );
        bool test( MPI_Status* pt_status );
        int wait( MPI_Status* pt_status );

        int         receiver( ) const { return m_recv_rank; }
        int         tag( ) const { return m_id_tag; }
        std::size_t size( ) const { return m_cur_size; }
        const void* data( ) const { return static_cast<const void*>( &m_arr_buffer[0] ); }

         SendBuffer::PackedData& pack( const SendBuffer::PackedData& data ) {
            m_cur_size += sizeof( std::size_t ) + data.size( );
            m_arr_pkg_data.push_back( data );
            return m_arr_pkg_data.back();
        }

    private:
        void           copyPackagedDataInBuffer( );
        int            m_recv_rank;
        int            m_id_tag;
        const MPI_Comm m_ref_comm;
        Datas          m_arr_pkg_data;
        Buffer         m_arr_buffer;
        std::size_t    m_cur_size;
        MPI_Request    m_request;
        bool           m_is_data_copied;
    };
    // ===========================================================================================
    // Définition de la mise en oeuvre du buffer d'envoi
    SendBuffer::Implementation::Implementation( int dest, int tag, const MPI_Comm& comm )
        : m_recv_rank( dest ),
          m_id_tag( tag ),
          m_ref_comm( comm ),
          m_arr_pkg_data( ),
          m_arr_buffer( ),
          m_cur_size( 0 ),
          m_is_data_copied(false) {
        m_arr_pkg_data.reserve( 16384);//1024 );
    }
    // -------------------------------------------------------------------------------------------
    void SendBuffer::Implementation::finalize()
    {
        if (not m_is_data_copied) {
            copyPackagedDataInBuffer( );
            m_is_data_copied = true;
        }
    }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::Implementation::isend( ) {
        if (not m_is_data_copied) {
            copyPackagedDataInBuffer( );
            m_is_data_copied = true;
        }
        const void* pt_data = &m_arr_buffer[0];
        int         length  = int( m_arr_buffer.size( ) );
        int ret = MPI_Issend( (void*)pt_data, length, MPI_BYTE, m_recv_rank, m_id_tag, m_ref_comm, &m_request );
        m_is_data_copied = false;
        return ret;
    }
    // -------------------------------------------------------------------------------------------
    bool SendBuffer::Implementation::test( MPI_Status* pt_status ) {
        int flag, ierr;
        if (m_arr_buffer.size() == 0) return true;

        if ( pt_status == NULL )
            ierr = MPI_Test( &m_request, &flag, MPI_STATUS_IGNORE );
        else
            ierr = MPI_Test( &m_request, &flag, pt_status );
#if defined( PCM_DEBUG )
        if ( ierr != MPI_SUCCESS ) {
            char errMsg[1024];
            int  lenStr;
            MPI_Error_string( ierr, errMsg, &lenStr );
            std::cerr << __PRETTY_FUNCTION__ << " : Erreur " << errMsg << std::endl;
            MPI_Abort( MPI_COMM_WORLD, ierr );
            exit( EXIT_FAILURE );
        }
#endif
        assert( ierr == MPI_SUCCESS );
        return ( flag != 0 );
    }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::Implementation::wait( MPI_Status* pt_status ) {
        int ierr = MPI_SUCCESS;
        if (m_arr_buffer.size() == 0) return ierr;
        if ( pt_status == NULL )
            ierr = MPI_Wait( &m_request, MPI_STATUS_IGNORE );
        else
            ierr = MPI_Wait( &m_request, pt_status );
#if defined( PCM_DEBUG )
        if ( ierr != MPI_SUCCESS ) {
            char errMsg[1024];
            int  lenStr;
            MPI_Error_string( ierr, errMsg, &lenStr );
            std::cerr << __PRETTY_FUNCTION__ << " : Erreur " << errMsg << std::endl;
            MPI_Abort( MPI_COMM_WORLD, ierr );
            exit( EXIT_FAILURE );
        }
#endif
        return ierr;
    }
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void SendBuffer::Implementation::copyPackagedDataInBuffer( ) {
        const std::size_t data_chunk                 = 2048;            // Copie par paquet de 2048 octets maximum
        const std::size_t min_size_for_parallel_copy = 8 * data_chunk;  // Taille minimal pour faire une copie parallèle
        if (m_cur_size != m_arr_buffer.size())
            std::vector<char>( m_cur_size ).swap( m_arr_buffer );
        std::size_t count = 0;
        // m_arr_buffer.resize( m_cur_size );
        Buffer::iterator itB = m_arr_buffer.begin( );

        for ( Datas::iterator it_d = m_arr_pkg_data.begin( ); it_d != m_arr_pkg_data.end( ); ++it_d ) {

            assert( itB < m_arr_buffer.end( ) );
            char* pt              = &( *itB );
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
    SendBuffer::SendBuffer( int recv_rank, int id_tag, const MPI_Comm& comm )
        : m_pt_implementation( new SendBuffer::Implementation( recv_rank, id_tag, comm ) ) {}
    // -------------------------------------------------------------------------------------------
    SendBuffer::SendBuffer( const SendBuffer& s_buf ) : m_pt_implementation( s_buf.m_pt_implementation ) {}
    // -------------------------------------------------------------------------------------------
    SendBuffer::~SendBuffer( ) {}
    // -------------------------------------------------------------------------------------------
    SendBuffer& SendBuffer::operator=( const SendBuffer& s_buf ) {
        if ( this != &s_buf ) { m_pt_implementation = s_buf.m_pt_implementation; }
        return *this;
    }
    // -------------------------------------------------------------------------------------------
    SendBuffer& SendBuffer::operator<<( const PackedData& data ) {
        m_pt_implementation->pack( data );
        return *this;
    }
    // -------------------------------------------------------------------------------------------
    SendBuffer::PackedData&
    SendBuffer::push_inplace_array( size_t size )
    {
        return m_pt_implementation->pack( SendBuffer::PackedData(size) );
    }
    // -------------------------------------------------------------------------------------------
    void SendBuffer::finalize_and_copy()
    {
        return m_pt_implementation->finalize();
    }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::isend( ) { return m_pt_implementation->isend( ); }
    // -------------------------------------------------------------------------------------------
    bool SendBuffer::test( MPI_Status* pt_status ) { return m_pt_implementation->test( pt_status ); }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::wait( MPI_Status* pt_status ) { return m_pt_implementation->wait( pt_status ); }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::receiver( ) const { return m_pt_implementation->receiver( ); }
    // -------------------------------------------------------------------------------------------
    int SendBuffer::tag( ) const { return m_pt_implementation->tag( ); }
    // -------------------------------------------------------------------------------------------
    std::size_t SendBuffer::size( ) const { return m_pt_implementation->size( ); }
    // -------------------------------------------------------------------------------------------
    const void* SendBuffer::data( ) const { return m_pt_implementation->data( ); }
}


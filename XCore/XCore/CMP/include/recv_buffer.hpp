// RecvBuffer.hpp
#ifndef _CMP_RECVBUFFER_HPP_
#define _CMP_RECVBUFFER_HPP_
#if defined(_WIN64)
# define __int64 long long
#endif
#include <mpi.h>
#include <cassert>
#include <string>
#include <memory>
using std::shared_ptr;
#include "vector_view.hpp"

namespace CMP {
    class RecvBuffer {
    public:
        RecvBuffer( int source = MPI_ANY_SOURCE, int id_tag = MPI_ANY_TAG, const MPI_Comm& comm = MPI_COMM_WORLD );
        RecvBuffer( const RecvBuffer& r_buf );
        ~RecvBuffer( );

        template <typename K>
        RecvBuffer& operator>>( K& val ) {
            std::size_t sz;
            const char* pt_data = unpack( sz );
            assert( sz == sizeof( K ) );
            val = *(const K*)( pt_data );
            return *this;
        }
        template <typename K>
        RecvBuffer& operator>>( vector_view<K>& array ) {
            std::size_t sz;
            const char* pt_data = unpack( sz );
            std::size_t sz_arr  = sz / sizeof( K );
            assert( sz_arr * sizeof( K ) == sz );
            array = vector_view<K>( (const K*)( pt_data ), sz_arr );
            return *this;
        }
        RecvBuffer& operator>>( std::string& str ) {
            std::size_t sz;
            const char* pt_data = unpack( sz );
            str                 = std::string( pt_data, sz );
            return *this;
        }

        int  irecv( );
        bool test( MPI_Status* pt_status = NULL );
        void wait( MPI_Status* pt_status = NULL );
        int         source( ) const;
        int         tag( ) const;
        int         size( ) const;
        const char* data( ) const;
        const char* current_data( ) const;

    private:
        const char* unpack( std::size_t& sz );
        class Implementation;
        shared_ptr<Implementation> m_pt_implementation;
    };
}
#endif

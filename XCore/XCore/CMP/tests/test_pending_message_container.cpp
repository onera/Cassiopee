#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include "pending_message_container.h"
#include "recv_buffer.hpp"
#include "send_buffer.hpp"

class DataTest {
public:
    static int nbDataTest;
    DataTest( const int* pt_arr, const std::string& label ) : m_pt_data( pt_arr ), m_label( label ) { nbDataTest += 1; }
    ~DataTest( ) { nbDataTest -= 1; };

    const int*         get_data( ) const { return m_pt_data; }
    const std::string& get_label( ) const { return m_label; }

private:
    const int*  m_pt_data;
    std::string m_label;
};
std::ostream& operator<<( std::ostream& out, const DataTest& d ) {
    out << "Data " << d.get_label( ) << " with data ";
    for ( int i = 0; i < 3; ++i ) out << d.get_data( )[i] << " ";
    return out;
}
int DataTest::nbDataTest = 0;
//
typedef CMP::PendingMsgContainer<CMP::RecvBuffer, DataTest> RecvQueue;
typedef CMP::PendingMsgContainer<CMP::SendBuffer> SendQueue;
// ---------------------------------------------------------------------------------------
int main( int nargs, char* argv[] ) {
    MPI_Init( &nargs, &argv );
    int szCom, rank;
    MPI_Comm_size( MPI_COMM_WORLD, &szCom );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    const int        nbData = 3;
    std::vector<int> array( nbData );
    int              ind = nbData * rank + 1;
    for ( int j = 0; j < nbData; ++j ) array[j] = ind++;
    int*      pt_array                          = &array[0];

    SendQueue squeue( szCom );
    RecvQueue rqueue( szCom );

    int                 ia = rank * 7;
    int                 ib;
    const std::size_t   nbElts = 1000000;
    std::vector<double> sdata( nbElts );
    for ( std::size_t        i = 0; i < nbElts; ++i ) sdata[i] = rank + i + 0.5;
    CMP::vector_view<double> v_data;
    CMP::vector_view<double> inpl_data;
    for ( int i = 0; i < szCom; ++i ) {
        if ( i != rank ) {
            squeue.push_back( new CMP::SendBuffer( i, 404 ) );
            rqueue.push_back( new CMP::RecvBuffer( i, 404 ), new DataTest( pt_array, std::string( "...data..." ) ) );
            CMP::SendBuffer& sbuf = squeue.back_message_buffer( );
            CMP::RecvBuffer& rbuf = rqueue.back_message_buffer( );
            rbuf.irecv( );
            sbuf << ia << sdata;
            CMP::SendBuffer::PackedData& pck_arr = sbuf.push_inplace_array(10*sizeof(double));
            sbuf.finalize_and_copy();
            double* pt_data = pck_arr.data<double>();
            for ( int j = 0; j < pck_arr.size()/sizeof(double); ++j )
                pt_data[j] = i+j+3.14;
            sbuf.isend( );
        }
    }

    while ( not rqueue.empty( ) ) {
        RecvQueue::iterator it = rqueue.get_first_complete_message( );
        if ( it != rqueue.end( ) ) {
            CMP::RecvBuffer& rbuf = ( *it ).get_message_buffer( );
            rbuf >> ib >> v_data >> inpl_data;
            std::cout << "local data : " << ( *it ).get_local_data( ) << std::endl;
            for ( std::size_t i = 0; i < std::min( v_data.size( ), 10UL ); ++i )
                std::cout << std::setprecision( 10 ) << v_data[i] << " ";
            if ( v_data.size( ) > 10 ) {
                std::cout << " ... ";
                for ( std::size_t i = v_data.size( ) - 10; i < v_data.size( ); ++i )
                    std::cout << std::setprecision( 10 ) << v_data[i] << " ";
            }
            std::cout << "\nreceived inplace data :\n";
            for ( int i = 0; i < 10; ++i )
                std::cout << inpl_data[i] << " "; 
            std::cout << std::endl;
            rqueue.pop( it );
        }
    }

    squeue.waitAll( );
    squeue.clear( );
    MPI_Finalize( );
    return EXIT_SUCCESS;
}

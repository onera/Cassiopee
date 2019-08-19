#include <cassert>
namespace CMP {
#define Queue_T Queue<SBuffer, RBuffer, Data>

    template <typename SBuffer, typename RBuffer, typename Data>
    Queue_T::Queue( int maxSize ) : m_queue( ) {
        if ( maxSize > 0 ) m_queue.reserve( maxSize );
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    void Queue_T::push_back( SBuffer* sbuffer, RBuffer* rbuffer, Data* dat ) {
        m_queue.resize( m_queue.size( ) + 1 );
        Queue_T::value& item = m_queue.back( );
        item.send_buff       = sbuffer;
        item.recv_buff       = rbuffer;
        item.data            = dat;
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    typename Queue_T::value& Queue_T::push_back( ) {
        m_queue.resize( m_queue.size( ) + 1 );
        return m_queue.back( );
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    void Queue_T::pop( const Queue_T::iterator& iter ) {
        ( *iter ).swap( m_queue.back( ) );
        m_queue.pop_back( );
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    void Queue_T::wait_complete_reception( ) {
        for ( Queue_T::iterator it = begin( ); it != end( ); ++it ) ( *it ).getRecvBuffer( ).wait( );
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    void Queue_T::wait_complete_sending( ) {
        for ( Queue_T::iterator it = begin( ); it != end( ); ++it ) ( *it ).getSendBuffer( ).wait( );
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    bool Queue_T::test_all( ) {
        bool flag = true;
        for ( Queue_T::iterator it = begin( ); ( it != end( ) ) && ( flag == true ); ++it ) {
            flag &= ( *it ).getRecvBuffer( ).test( ) && ( *it ).getSendBuffer( ).test( );
        }
        return flag;
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    bool Queue_T::test_complete_sending( ) {
        bool flag = true;
        for ( Queue_T::iterator it = begin( ); ( it != end( ) ) && ( flag == true ); ++it ) {
            flag &= ( *it ).getSendBuffer( ).test( );
        }
        return flag;
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    bool Queue_T::test_complete_reception( ) {
        bool flag = true;
        for ( Queue_T::iterator it = begin( ); ( it != end( ) ) && ( flag == true ); ++it ) {
            flag &= ( *it ).getRecvBuffer( ).test( );
        }
        return flag;
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    typename Queue_T::iterator Queue_T::test_any_sending( ) {
        bool              flag = false;
        Queue_T::iterator it;
        for ( it = begin( ); it != end( ); ++it ) {
            flag = ( *it ).getSendBuffer( ).test( );
            if ( flag ) break;
        }
        return it;
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    typename Queue_T::iterator Queue_T::test_any_reception( ) {
        bool              flag = false;
        Queue_T::iterator it;
        for ( it = begin( ); it != end( ); ++it ) {
            flag = ( *it ).getRecvBuffer( ).test( );
            if ( flag ) break;
        }
        return it;
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    typename Queue_T::value& Queue_T::operator[]( int i ) {
        assert( i >= 0 );
        assert( i < size( ) );
        return m_queue[i];
    }
    // --------------------------------------------------------------
    template <typename SBuffer, typename RBuffer, typename Data>
    const typename Queue_T::value& Queue_T::operator[]( int i ) const {
        assert( i >= 0 );
        assert( i < size( ) );
        return m_queue[i];
    }
#undef Queue_T
}

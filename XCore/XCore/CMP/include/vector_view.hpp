// vector_view.hpp
#ifndef _CMP_VECTOR_VIEW_HPP_
#define _CMP_VECTOR_VIEW_HPP_
#include <stdexcept>
#include <vector>

namespace CMP {
    template <typename K>
    class vector_view {
    public:
        typedef K                                     value_type;
        typedef const K*                              pointer;
        typedef pointer                               const_pointer;
        typedef const K&                              reference;
        typedef reference                             const_reference;
        typedef pointer                               iterator;
        typedef const_pointer                         const_iterator;
        typedef std::reverse_iterator<iterator>       reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        typedef std::size_t                           size_type;
        typedef std::ptrdiff_t                        difference_type;

        vector_view( ) : m_pt_data( NULL ), m_size( 0 ) {}
        vector_view( const K* data, size_type sz ) : m_pt_data( data ), m_size( sz ) {}
        vector_view( const std::vector<K>& arr ) : m_pt_data( &arr[0] ), m_size( arr.size( ) ) {}
        template <typename InpRandIter>
        vector_view( InpRandIter beg, InpRandIter end ) : m_pt_data( &( *beg ) ), m_size( end - beg ) {}
        vector_view( const vector_view& other ) : m_pt_data( other.m_pt_data ), m_size( other.m_size ) {}
        ~vector_view( ) {}

        vector_view& operator=( const vector_view& other ) {
            m_pt_data = other.m_pt_data;
            m_size    = other.m_size;
            return *this;
        }

        iterator               begin( ) { return m_pt_data; }
        iterator               end( ) { return m_pt_data + m_size; }
        const_iterator         begin( ) const { return cbegin( ); }
        const_iterator         end( ) const { return cend( ); }
        const_iterator         cbegin( ) const { return m_pt_data; }
        const_iterator         cend( ) const { return m_pt_data + m_size; }
        reverse_iterator       rbegin( ) { return end( ) - 1; }
        reverse_iterator       rend( ) { return begin( ) - 1; }
        const_reverse_iterator crbegin( ) { return cend( ) - 1; }
        const_reverse_iterator crend( ) { return cbegin( ) - 1; }
        const_reverse_iterator rbegin( ) const { return crbegin( ); }
        const_reverse_iterator rend( ) const { return crend( ); }

        size_type size( ) const { return m_size; }
        size_type length( ) const { return m_size; }
        size_type max_size( ) const { return m_size; }
        bool      empty( ) const { return m_size == 0; }

        const_reference operator[]( const size_type pos ) { return m_pt_data[pos]; }
        const K& at( size_type pos ) const {
            if ( pos >= size( ) ) throw std::out_of_range( "Wrong index" );
            return m_pt_data[pos];
        }
        const K&      front( ) const { return m_pt_data[0]; }
        const K&      back( ) const { return m_pt_data[size( ) - 1]; }
        const_pointer data( ) const { return m_pt_data; }

        void clear( ) { *this = vector_view( ); }
        void swap( vector_view& s ) {
            std::swap( m_pt_data, s.m_pt_data );
            std::swap( m_size, s.m_size );
        }

        operator std::vector<K>( ) const { return std::vector<K>( begin( ), end( ) ); }

    private:
        const_pointer m_pt_data;
        size_type     m_size;
    };
}
#endif

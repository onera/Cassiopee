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
#ifndef _KCORE_VECTOR_VIEW_HPP_
#define _KCORE_VECTOR_VIEW_HPP_
#include <stdexcept>
#include <vector>
#include <cassert>

namespace K_MEMORY {
    /**
     * @brief      Same interface as vector but not own the pointed memory.
     *
     * @tparam     K     Kind of value contained by the vector_view
     */
    template <typename K>
    class vector_view {
    public:
        typedef K                                     value_type;
        typedef K*                                    pointer;
        typedef const K*                              const_pointer;
        typedef K&                                    reference;
        typedef const K&                              const_reference;
        typedef pointer                               iterator;
        typedef const_pointer                         const_iterator;
        typedef std::reverse_iterator<iterator>       reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        typedef std::size_t                           size_type;
        typedef std::ptrdiff_t                        difference_type;

        vector_view( ) : m_pt_data( NULL ), m_size( 0 ) {}
        vector_view( K* data, size_type sz ) : m_pt_data( data ), m_size( sz ) {}
        vector_view( std::vector<K>& arr ) : m_pt_data( &arr[0] ), m_size( arr.size( ) ) {}
        template <typename InpRandIter>
        vector_view( InpRandIter beg, InpRandIter end ) : m_pt_data( &( *beg ) ), m_size( end - beg ) {}
        vector_view( const vector_view& other ) : m_pt_data( other.m_pt_data ), m_size( other.m_size ) {}
        ~vector_view( ) {}

        vector_view& operator=( const vector_view& other ) {
            m_pt_data = const_cast<pointer>(other.m_pt_data);
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

        reference operator[] ( const size_type pos ) {
            assert(pos < size());
            return m_pt_data[pos]; }
        const_reference operator[] ( const size_type pos ) const { 
            assert(pos < size());
            return m_pt_data[pos]; }
        reference at( size_type pos ) {
            if ( pos >= size( ) ) throw std::out_of_range( "Wrong index" );
            return m_pt_data[pos];
        }
        const_reference at( size_type pos ) const {
            if ( pos >= size( ) ) throw std::out_of_range( "Wrong index" );
            return m_pt_data[pos];
        }
        reference       front( ) { return m_pt_data[0]; }
        reference       back( ) { return m_pt_data[size( ) - 1]; }
        const_reference front( ) const { return m_pt_data[0]; }
        const_reference back( ) const { return m_pt_data[size( ) - 1]; }
        pointer data( ) { return m_pt_data; }
        const_pointer data( ) const { return m_pt_data; }

        void clear( ) { *this = vector_view( ); }
        void swap( vector_view& s ) {
            std::swap( m_pt_data, s.m_pt_data );
            std::swap( m_size, s.m_size );
        }

        operator std::vector<K>( ) const { return std::vector<K>( begin( ), end( ) ); }

    private:
        pointer m_pt_data;
        size_type     m_size;
    };
}
#endif

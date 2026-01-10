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
#ifndef _CMP_SHARED_PTR_HPP_
#define _CMP_SHARED_PTR_HPP_
#include <cassert>
#include <cstdlib>
#include <utility>

namespace CMP {
    template <typename K>
    class shared_ptr {
    public:
        typedef K        element_type;
        typedef K*       pointer;
        typedef const K* const_pointer;
        typedef std::pair<std::size_t, pointer> value;

        // Constructors and destructor
        // ===========================
        shared_ptr( ) : m_nref_pointer( new value ) {
            set_pointer( NULL );
            set_count( 0 );
        }
        // .....................................
        template <typename L>
        shared_ptr( L* ptr ) : m_nref_pointer( new value ) {
            set_pointer( static_cast<pointer>( ptr ) );
            set_count( 1 );
        }
        // .....................................
        shared_ptr( shared_ptr& spt ) : m_nref_pointer( spt.m_nref_pointer ) { increference( ); }
        // .....................................
        shared_ptr( const shared_ptr& spt ) : m_nref_pointer( spt.m_nref_pointer ) { increference( ); }
        // .....................................
        ~shared_ptr( ) { dereference( ); }

        // Opérateurs
        // ================================================================================
        shared_ptr& operator=( const shared_ptr& pt ) {
            if ( this != &pt ) {
                dereference( );
                m_nref_pointer = pt.m_nref_pointer;
                increference( );
            }
            return *this;
        }

        void reset( ) {
            dereference( );
            m_nref_pointer = new value;
            set_pointer( NULL );
            set_count( 0 );
        }
        template <typename Y>
        void reset( Y* ptr ) {
            dereference( );
            m_nref_pointer = new value;
            set_pointer( const_cast<pointer>( ptr ) );
            set_count( 1 );
        }
        void swap( shared_ptr& r ) {
            value* tmp       = m_nref_pointer;
            m_nref_pointer   = r.m_nref_pointer;
            r.m_nref_pointer = tmp;
        }

        pointer       get( ) { return get_pointer( ); }
        const_pointer get( ) const { return get_pointer( ); }

        element_type& operator*( ) { return *get_pointer( ); }
        const element_type& operator*( ) const { return *get_pointer( ); }
        pointer operator->( ) { return get_pointer( ); }
        const_pointer operator->( ) const { return get_pointer( ); }

        std::size_t use_count( ) const { return get_count( ); }

        operator bool( ) { return get_pointer( ) == NULL; }

    private:
        pointer       get_pointer( ) { return m_nref_pointer == NULL ? NULL : m_nref_pointer->second; }
        const_pointer get_pointer( ) const { return m_nref_pointer == NULL ? NULL : m_nref_pointer->second; }
        void set_pointer( K* ptK ) {
            assert( m_nref_pointer != NULL );
            m_nref_pointer->second = ptK;
        }

        void set_count( std::size_t nbref ) {
            assert( m_nref_pointer != NULL );
            m_nref_pointer->first = nbref;
        }
        std::size_t get_count( ) const { return m_nref_pointer == NULL ? 0 : m_nref_pointer->first; }
        void        increference( ) {
            assert( m_nref_pointer != NULL );
            m_nref_pointer->first += 1;
        }
        void dereference( ) {
            if ( m_nref_pointer == NULL ) return;
            if ( m_nref_pointer->first > 0 ) {
                m_nref_pointer->first -= 1;
                if ( m_nref_pointer->first == 0 ) {
                    delete m_nref_pointer->second;
                    delete m_nref_pointer;
                    m_nref_pointer = NULL;
                }
            }
        }

        mutable value* m_nref_pointer;
    };
    template <typename T, typename U>
    bool operator==( const shared_ptr<T>& lhs, const shared_ptr<T>& rhs ) {
        return lhs.get( ) == rhs.get( );
    }
    template <typename T, typename U>
    bool operator!=( const shared_ptr<T>& lhs, const shared_ptr<T>& rhs ) {
        return lhs.get( ) != rhs.get( );
    }
    template <typename T, typename U>
    bool operator<( const shared_ptr<T>& lhs, const shared_ptr<T>& rhs ) {
        return lhs.get( ) < rhs.get( );
    }
    template <typename T, typename U>
    bool operator>( const shared_ptr<T>& lhs, const shared_ptr<T>& rhs ) {
        return lhs.get( ) > rhs.get( );
    }
    template <typename T, typename U>
    bool operator<=( const shared_ptr<T>& lhs, const shared_ptr<T>& rhs ) {
        return lhs.get( ) <= rhs.get( );
    }
    template <typename T, typename U>
    bool operator>=( const shared_ptr<T>& lhs, const shared_ptr<T>& rhs ) {
        return lhs.get( ) >= rhs.get( );
    }

    template <typename T, typename U>
    bool operator==( const shared_ptr<T>& lhs, const U* rhs ) {
        return lhs.get( ) == rhs;
    }
    template <typename T, typename U>
    bool operator==( const T* lhs, const shared_ptr<U>& rhs ) {
        return lhs == rhs.get( );
    }
    template <typename T, typename U>
    bool operator!=( const shared_ptr<T>& lhs, const U* rhs ) {
        return lhs.get( ) != rhs;
    }
    template <typename T, typename U>
    bool operator!=( const T* lhs, const shared_ptr<U>& rhs ) {
        return lhs != rhs.get( );
    }
    template <typename T, typename U>
    bool operator>( const shared_ptr<T>& lhs, const U* rhs ) {
        return lhs.get( ) > rhs;
    }
    template <typename T, typename U>
    bool operator>( const T* lhs, const shared_ptr<U>& rhs ) {
        return lhs > rhs.get( );
    }
    template <typename T, typename U>
    bool operator<( const shared_ptr<T>& lhs, const U* rhs ) {
        return lhs.get( ) < rhs;
    }
    template <typename T, typename U>
    bool operator<( const T* lhs, const shared_ptr<U>& rhs ) {
        return lhs < rhs.get( );
    }
    template <typename T, typename U>
    bool operator>=( const shared_ptr<T>& lhs, const U* rhs ) {
        return lhs.get( ) >= rhs;
    }
    template <typename T, typename U>
    bool operator>=( const T* lhs, const shared_ptr<U>& rhs ) {
        return lhs >= rhs.get( );
    }
    template <typename T, typename U>
    bool operator<=( const shared_ptr<T>& lhs, const U* rhs ) {
        return lhs.get( ) <= rhs;
    }
    template <typename T, typename U>
    bool operator<=( const T* lhs, const shared_ptr<U>& rhs ) {
        return lhs <= rhs.get( );
    }
    template <class T, class U, class V>
    std::basic_ostream<U, V>& operator<<( std::basic_ostream<U, V>& os, const shared_ptr<T>& ptr ) {
        os << (void*)ptr.get( );
        return os;
    }
}
#endif

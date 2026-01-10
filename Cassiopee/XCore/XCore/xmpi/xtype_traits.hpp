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
#ifndef _XCORE_XMPI_TYPE_TRAITS_HPP_
#define _XCORE_XMPI_TYPE_TRAITS_HPP_
#include <type_traits>

namespace xcore
{
    template <bool B, class T = void>
    using enable_if = std::enable_if<B, T>;

    template <typename C>
    struct is_allocatable_container {
        static const bool value = ( std::is_constructible<C, std::size_t>::value ) && ( not std::is_arithmetic<C>::value );
    };

    template <typename K>
    struct is_direct_access_container {
        template <typename C>
        static constexpr decltype( std::declval<C>().operator[]( 0 ), bool() ) test( int )
        {
            return true;
        }
        template <typename C>
        static constexpr bool test( ... )
        {
            return false;
        }

        static constexpr bool value = test<K>( int{} );
    };

    template <typename T>
    struct has_const_iterator {
    private:
        typedef char yes;
        typedef struct {
            char array[ 2 ];
        } no;

        template <typename C>
        static yes test( typename C::const_iterator * );
        template <typename C>
        static no test( ... );

    public:
        static const bool value = sizeof( test<T>( 0 ) ) == sizeof( yes );
        typedef T type;
    };

    /**
   Detect if the type T has a begin and end methods
*/
    template <typename T>
    struct has_begin_end {
        template <typename C>
        static char (
            &f( typename std::enable_if<
                std::is_same<decltype( static_cast<typename C::const_iterator (
                                           C::* )() const>( &C::begin ) ),
                             typename C::const_iterator ( C::* )() const>::value,
                void>::type * ) )[ 1 ];

        template <typename C>
        static char ( &f( ... ) )[ 2 ];

        template <typename C>
        static char (
            &g( typename std::enable_if<
                std::is_same<decltype( static_cast<typename C::const_iterator (
                                           C::* )() const>( &C::end ) ),
                             typename C::const_iterator ( C::* )() const>::value,
                void>::type * ) )[ 1 ];

        template <typename C>
        static char ( &g( ... ) )[ 2 ];

        static bool const beg_value = sizeof( f<T>( 0 ) ) == 1;
        static bool const end_value = sizeof( g<T>( 0 ) ) == 1;
    };

    /**
   Detect if T is a container : if T has const_iterator AND begin and end
   methods
*/
    template <typename T>
    struct is_container
        : std::integral_constant<bool, has_const_iterator<T>::value && has_begin_end<T>::beg_value && has_begin_end<T>::end_value> {
    };

#define DEFINE_HAS_SIGNATURE(traitsName, funcName, signature)               \
    template <typename U>                                                   \
    class traitsName                                                        \
    {                                                                       \
    private:                                                                \
        template<typename T, T> struct helper;                              \
        template<typename T>                                                \
        static std::uint8_t check(helper<signature, &funcName>*);           \
        template<typename T> static std::uint16_t check(...);               \
    public:                                                                 \
        static                                                              \
        constexpr bool value = sizeof(check<U>(0)) == sizeof(std::uint8_t); \
    }

}  // namespace xcore
#endif

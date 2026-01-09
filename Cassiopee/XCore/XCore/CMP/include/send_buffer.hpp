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
// SendBuffer.hpp
#ifndef _CMP_SENDBUFFER_HPP_
#define _CMP_SENDBUFFER_HPP_
#include <iostream>
#include <string>
#include <memory>
#include <cassert>
using std::shared_ptr;

#include "xmpi/xmpi.hpp"
#include "vector_view.hpp"

namespace CMP {
    class SendBuffer {
    public:
        using value_t = unsigned char;
        class PackedData {
        public:
            using value_t = typename SendBuffer::value_t;
            typedef value_t*       iterator;
            typedef const value_t* const_iterator;

            template <typename K>
            PackedData( const K& val )
                : m_nbBytes( sizeof( K ) ), m_pt_data( (value_t*)new K( val ) ), m_is_owner( true ), m_must_copy(true) {}

            template <typename K>
            PackedData( const std::size_t nbElts, const K* pt_data )
                : m_nbBytes( nbElts * sizeof( K ) ), m_pt_data( (value_t*)pt_data ), m_is_owner( pt_data == NULL ), m_must_copy(true) {
                if ( m_is_owner ) m_pt_data = (unsigned char*)new K[nbElts];
            }

            PackedData( const std::size_t nbElts )
                : m_nbBytes( nbElts ), m_pt_data( NULL ), m_is_owner( false ), m_must_copy(false) 
            {
            }

            template <typename InpIterator>
            PackedData( InpIterator beg, InpIterator end ) : m_nbBytes( 0 ), m_pt_data( NULL ), m_is_owner( false ), m_must_copy(true) {
                value_t* pt_beg = (value_t*)&( *beg );
                value_t* pt_end = (value_t*)&( *end );
                m_nbBytes    = pt_end - pt_beg;
                m_pt_data    = pt_beg;
            }

            PackedData( const std::string& str );

            /** Copie plus dans l'esprit du déplacement */
            PackedData( const PackedData& data )
                : m_nbBytes( data.m_nbBytes ), m_pt_data( data.m_pt_data ), m_is_owner( data.m_is_owner ) {
                data.m_is_owner = false; m_must_copy = data.m_must_copy;
            }
            /** Destructeur
            */
            ~PackedData( ) {
                if ( m_is_owner ) delete m_pt_data;
            }

            PackedData& operator=( const PackedData& data ) {
                if ( this != &data ) {
                    if ( m_is_owner ) delete m_pt_data;
                    m_nbBytes       = data.m_nbBytes;
                    m_pt_data       = data.m_pt_data;
                    m_is_owner      = data.m_is_owner;
                    m_must_copy     = data.m_must_copy;
                    data.m_is_owner = false;
                }
                return *this;
            }

            /** Adresse des données à sérialiser */
            template <typename K>
            K* data( ) {
                return (K*)( m_pt_data );
            }
            /** Adresse des données à sérialiser
            */
            template <typename K>
            const K* data( ) const {
                return (const K*)( m_pt_data );
            }

            template<typename K> void set_data( K* addr ) { assert(m_pt_data == NULL); m_pt_data = (value_t*)addr; }

            std::size_t    size( ) const { return m_nbBytes; }  // Taille en octet des données sérialisés
            iterator       begin( ) { return static_cast<iterator>( m_pt_data ); }
            const_iterator begin( ) const { return static_cast<const_iterator>( m_pt_data ); }
            iterator       end( ) { return static_cast<iterator>( m_pt_data ) + size( ); }
            const_iterator end( ) const { return static_cast<const_iterator>( m_pt_data ) + size( ); }
            bool must_copy() const { return m_must_copy; }

        private:
            std::size_t  m_nbBytes;   // Taille en octet des données
            value_t*        m_pt_data;   // Adresse
            mutable bool m_is_owner;  // Données à détruire par la suite ?
            bool         m_must_copy; // Uniquement un espace à réserver ( faux ) ou données qu'on doit copier ( true ) ?
        };
        // -----------------------------------------------------------------------------
        SendBuffer( int recv_rank, int id_tag = xcore::any_tag );
        SendBuffer( int recv_rank, int id_tag, const xcore::communicator& comm );
        SendBuffer( const SendBuffer& s_buf );
        ~SendBuffer( );

        SendBuffer& operator=( const SendBuffer& s_buf );

        SendBuffer& operator<<( const PackedData& data );
        template <typename K>
        SendBuffer& operator<<( const K& val ) {
            *this << PackedData( val );
            return *this;
        }
        template <typename K>
        SendBuffer& operator<<( const std::vector<K>& tab ) {
            *this << PackedData( tab.begin( ), tab.end( ) );
            return *this;
        }
        template <typename K>
        SendBuffer& operator<<( const vector_view<K>& tab ) {
            *this << PackedData( tab.begin( ), tab.end( ) );
            return *this;
        }

        PackedData& push_inplace_array( size_t size );

        void finalize_and_copy();
        int  isend( );  // Envoie asynchrone du buffer, retour un entier pour signaler une erreur éventuelle
        bool test( xcore::status* pt_status = NULL );  // Renvoie vrai si l'envoi est fini.
        int wait(  xcore::status* pt_status = NULL );   // Attend que l'envoi se finisse
        int         receiver( ) const;              // Retourne le rang du receveur
        int         tag( ) const;                   // Retourne l'identifiant du message envoyé par le buffer
        std::size_t size( ) const;                  // Retourne la taille en octet pris par le buffer d'envoi
        const value_t* data( ) const;                  // Renvoie l'adresse du buffer d'envoi ( pour déboguage )
        void clear();
    private:  // ######################### PRIVATE PART BELOW #############################################
        // PIMPL template
        class Implementation;
        shared_ptr<Implementation> m_pt_implementation;
    };
    inline std::ostream& operator<<( std::ostream& out, const SendBuffer& s_buf ) {
        out << "SendBuffer{ receiver : " << s_buf.receiver( ) << ", id tag : " << s_buf.tag( )
            << ", size of buffer : " << s_buf.size( ) << ", addr buffer : " << (void*)s_buf.data( ) << " }";
        return out;
    }
}

#endif

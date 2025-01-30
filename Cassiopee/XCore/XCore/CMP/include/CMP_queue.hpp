/*    
    Copyright 2013-2025 Onera.

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
#ifndef _CMP_CMP_QUEUE_HPP_
#define _CMP_CMP_QUEUE_HPP_
#include <cassert>
#include <vector>

namespace CMP {
    /** \brief
         Queue gérant un buffer d'envoie, un de réception et un de données associées à ces deux échanges
     */
    template <typename SBuffer, typename RBuffer, typename Data>
    class Queue {
    public:
        // Structure permettant d'associer un buffer d'envoie avec un buffer de réception et des données internes
        // ======================================================================================================
        struct item_queue {
            SBuffer* send_buff;
            RBuffer* recv_buff;
            Data*    data;

            item_queue( ) : send_buff( NULL ), recv_buff( NULL ), data( NULL ) {}
            item_queue( SBuffer* sbuff, RBuffer* rbuff, Data* d ) : send_buff( sbuff ), recv_buff( rbuff ), data( d ) {}
            ~item_queue( ) {
                delete send_buff;
                delete recv_buff;
                delete data;
            }

            SBuffer& getSendBuffer( ) { return *send_buff; }
            RBuffer& getRecvBuffer( ) { return *recv_buff; }
            Data&    getData( ) { return *data; }

            void setSendBuffer( SBuffer* buf ) { send_buff = buf; }
            void setRecvBuffer( RBuffer* buf ) { recv_buff = buf; }
            void setData( Data* dat ) { data = dat; }

            void swap( item_queue& it ) {
                SBuffer* stmp = send_buff;
                RBuffer* rtmp = recv_buff;
                Data*    dtmp = data;
                send_buff     = it.send_buff;
                recv_buff     = it.recv_buff;
                data          = it.data;
                it.send_buff  = stmp;
                it.recv_buff  = rtmp;
                it.data       = dtmp;
            }
            item_queue( const item_queue& q ) : send_buff( NULL ), recv_buff( NULL ), data( NULL ) {
                assert( q.send_buff == NULL );
                assert( q.recv_buff == NULL );
                assert( q.data == NULL );
            }
            item_queue& operator=( const item_queue& q ) {
                assert( q.send_buff == NULL );
                assert( q.recv_buff == NULL );
                assert( q.data == NULL );
                delete send_buff;
                delete recv_buff;
                delete data;
                send_buff = NULL;
                recv_buff = NULL;
                data      = NULL;
                return *this;
            }

#if __cplusplus > 199711L
            template <typename... K>
            SBuffer& make_sendBuffer( K&... args ) {
                send_buff = new SBuffer( std::forward( args... ) );
                return *send_buff;
            }
            template <typename... K>
            RBuffer& make_recvBuffer( K&... args ) {
                recv_buff = new RBuffer( std::forward( args... ) );
                return *recv_buff;
            }
            template <typename... K>
            Data& make_data( K&... args ) {
                data = new Data( std::forward( args... ) );
                return data;
            }
            item_queue( item_queue&& iq ) {
                send_buff    = iq.send_buff;
                iq.send_buff = nullptr;
                recv_buff    = iq.recv_buff;
                iq.recv_buff = nullptr;
                data         = iq.data;
                iq.data      = nullptr;
            }
            item_queue& operator=( item_queue&& iq ) {
                if ( this != &iq ) {
                    delete send_buff;
                    delete recv_buff;
                    delete data;
                    send_buff    = iq.send_buff;
                    iq.send_buff = nullptr;
                    recv_buff    = iq.recv_buff;
                    iq.recv_buff = nullptr;
                    data         = iq.data;
                    iq.data      = nullptr;
                }
                return *this;
            }
#endif
        };
        // =======================================================================
        // Début de la déclaration de la classe Queue :
        typedef typename Queue<SBuffer, RBuffer, Data>::item_queue value;
        typedef typename std::vector<value>::iterator       iterator;
        typedef typename std::vector<value>::const_iterator const_iterator;

        Queue( ) : m_queue( ) {}
        Queue( int maxSize );
#if __cplusplus > 199711L
        Queue( const Queue& queue ) = delete;
        Queue& operator=( const Queue& queue ) = delete;

        Queue( Queue&& queue ) = delete;
        Queue& operator=( Queue&& queue ) = delete;
#endif
        ~Queue( ) {}

        void push_back( SBuffer* sbuffer, RBuffer* rbuffer, Data* dat );
        value& push_back( );
        void pop( const iterator& it );

        bool empty( ) { return m_queue.empty( ); }
        int  size( ) { return m_queue.size( ); }
        void clear( ) { m_queue.clear( ); }

        void wait_all( ) {
            wait_complete_reception( );
            wait_complete_sending( );
        }
        void wait_complete_sending( );
        void wait_complete_reception( );
        bool test_all( );
        bool test_complete_sending( );
        bool test_complete_reception( );

        // Renvoie end() si aucun buffer de prêt
        // sinon renvoie un itérateur sur le premier buffer prêt de trouvé
        iterator test_any_reception( );
        iterator test_any_sending( );

        value& operator[]( int i );
        const value& operator[]( int i ) const;
        value&       back( ) { return m_queue.back( ); }
        const value& back( ) const { return m_queue.back( ); }
        value&       front( ) { return m_queue.front( ); }
        const value& front( ) const { return m_queue.front( ); }

        iterator begin( ) { return m_queue.begin( ); }
#if __cplusplus > 199711L
        const_iterator cbegin( ) { return m_queue.cbegin( ); }
#endif
        const_iterator begin( ) const { return m_queue.begin( ); }

        iterator end( ) { return m_queue.end( ); }
#if __cplusplus > 199711L
        const_iterator cend( ) { return m_queue.cend( ); }
#endif
        const_iterator end( ) const { return m_queue.end( ); }

    private:
#if __cplusplus < 201103L
        Queue( const Queue& queue );
        Queue& operator=( const Queue& queue );
#endif
        std::vector<value> m_queue;
    };
}

#endif

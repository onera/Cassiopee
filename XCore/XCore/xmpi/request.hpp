#ifndef _XCORE_XMPI_REQUEST_HPP_
#define _XCORE_XMPI_REQUEST_HPP_

#include "xmpi/status.hpp"

namespace xcore
{
#       if defined(_MPI)        
        class request
        {
        public:
            request() = default;
            request( const MPI_Request& req ) : m_req( req )
            {}
            request(const request& ) = default;
            ~request() = default;

            request& operator = ( const request& ) = default;
            bool test() {
                int flag; MPI_Test(&m_req, &flag, &m_status);
                return (flag != 0);
            }
            void wait() {MPI_Wait( &m_req, &m_status);}
            void cancel() { MPI_Cancel( &m_req); }
            status get_status() const { return status(m_status); }
        private:
            MPI_Request m_req;
            MPI_Status  m_status;
        };
#       else
        class request
        {
        public:
            request() = default;
            request( const request& ) = default;
            ~request()= default;

            request& operator = ( const request& ) = default;

            bool test() { return true; }
            void wait() {}
            void cancel() {}
            status get_status() const { return {0,0,0}; }
        };
#       endif
}

#endif

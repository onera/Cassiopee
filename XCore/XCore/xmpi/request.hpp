#ifndef _XCORE_XMPI_REQUEST_HPP_
#define _XCORE_XMPI_REQUEST_HPP_
#include <memory>
using std::shared_ptr;
#include "xmpi/status.hpp"

namespace xcore
{
#       if defined(_MPI)        
        class request
        {
        public:
            request() : m_req(std::make_shared<MPI_Request>()) {}
            request(const request& ) = default;
            ~request() = default;

            request& operator = ( const request& ) = default;
            bool test() {
                int flag; MPI_Test(m_req.get(), &flag, &m_status);
                return (flag != 0);
            }
            void wait() {MPI_Wait( m_req.get(), &m_status);}
            void cancel() { MPI_Cancel( m_req.get()); }
            status get_status() const { return status(m_status); }
            std::shared_ptr<MPI_Request> m_req;
        private:
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

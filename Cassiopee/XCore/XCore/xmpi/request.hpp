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

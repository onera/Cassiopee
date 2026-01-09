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
#ifndef _MPI

#include "Logger/log_from_distributed_file.hpp"
#include "xmpi/context.hpp"
#include "Logger/log_from_distributed_file.hpp"

namespace xcore
{
    communicator *context::pt_global_com = nullptr;
    K_LOGGER::logger*       pt_xmpi_logg = nullptr;

    context::context() : m_provided(context::thread_support::Single)
    {}
    // ...............................................................................................
    context::context( int &nargc, char *argv[], bool isMultithreaded )
        : context::context( nargc, argv,
                            ( isMultithreaded
                                  ? context::thread_support::Multiple
                                  : context::thread_support::Single ) ) {}
    // ...............................................................................................
    context::context( int &nargc, char *argv[],
                      context::thread_support thread_level_support )
        : m_provided( thread_level_support )
    {
    }
    // ...............................................................................................
    context::~context()
    {
    }
    // ...............................................................................................
    communicator &context::globalCommunicator()
    {
        if ( pt_global_com == nullptr ) pt_global_com = new communicator;
        return *pt_global_com;
    }
    // ...............................................................................................
    K_LOGGER::logger& context::logger()
    {
        if ( pt_xmpi_logg == nullptr ) 
        {
            pt_xmpi_logg = new K_LOGGER::logger;
            pt_xmpi_logg->subscribe(new K_LOGGER::log_from_distributed_file(K_LOGGER::logger::listener::listen_for_all, 
                                                                            "XMPI_Output", 0));
        }
        return *pt_xmpi_logg;
    }
}
#endif

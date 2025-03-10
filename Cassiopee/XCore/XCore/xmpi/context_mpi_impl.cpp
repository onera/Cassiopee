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
#if defined (_MPI)
#include "xmpi/x_mpi.h"
#include "xmpi/context.hpp"
#include "Logger/log_from_distributed_file.hpp"

namespace xcore
{
    communicator *context::pt_global_com = nullptr;
    K_LOGGER::logger*  pt_xmpi_logg = nullptr;

    context::context()
    {
        int provided;
        int initialized;
        has_initialized_parallel = true;
        MPI_Initialized(&initialized);
        if ( not initialized )
        {
            MPI_Init_thread(NULL,NULL, MPI_THREAD_MULTIPLE, &provided);
            has_initialized_parallel = false;
        }
        else
            MPI_Query_thread( &provided );
        switch ( provided ) {
            case MPI_THREAD_SINGLE:
                m_provided = context::thread_support::Single;
                break;
            case MPI_THREAD_FUNNELED:
                m_provided = context::thread_support::Funneled;
                break;
            case MPI_THREAD_SERIALIZED:
                m_provided = context::thread_support::Serialized;
                break;
            default:
                m_provided = context::thread_support::Multiple;
        }
    }
    // ...............................................................................................    
    context::context( int &nargc, char *argv[], bool isMultithreaded )
        : context::context(
              nargc, argv, ( isMultithreaded ? context::thread_support::Multiple : context::thread_support::Single ) ) {
    }
    // ...............................................................................................
    context::context( int &nargc, char *argv[], context::thread_support thread_level_support )
        : m_provided( thread_level_support ) {
        int initialized;
        has_initialized_parallel = true;
        MPI_Initialized(&initialized);
        if ( initialized ) // MPI a deja ete initialise auparavant
        {
            has_initialized_parallel = false;
            int provided;
            MPI_Query_thread( &provided );
            if ( provided < thread_level_support )
                throw std::runtime_error("MPI is already initialized, but with a thread support lower than the required !");
            switch ( provided ) {
                case MPI_THREAD_SINGLE:
                    m_provided = context::thread_support::Single;
                    break;
                case MPI_THREAD_FUNNELED:
                    m_provided = context::thread_support::Funneled;
                    break;
                case MPI_THREAD_SERIALIZED:
                    m_provided = context::thread_support::Serialized;
                   break;
                default:
                    m_provided = context::thread_support::Multiple;
            }
            return;
        }
        if ( thread_level_support == context::thread_support::Single )
            if ( (nargc == 0) and (argv == nullptr) )
                MPI_Init( NULL, NULL );
            else
                MPI_Init( &nargc, &argv );
        else {
            int level_support;
            switch ( thread_level_support ) {
            case thread_support::Funneled:
                level_support = MPI_THREAD_FUNNELED;
                break;
            case thread_support::Serialized:
                level_support = MPI_THREAD_SERIALIZED;
                break;
            default:
                level_support = MPI_THREAD_MULTIPLE;
            }
            int provided;
            if ( (nargc == 0) and (argv == nullptr) )
                MPI_Init_thread( NULL, NULL, level_support, &provided );
            else
                MPI_Init_thread( &nargc, &argv, level_support, &provided );
            switch ( provided ) {
                case MPI_THREAD_SINGLE:
                    m_provided = context::thread_support::Single;
                    break;
                case MPI_THREAD_FUNNELED:
                    m_provided = context::thread_support::Funneled;
                    break;
             case MPI_THREAD_SERIALIZED:
                    m_provided = context::thread_support::Serialized;
                    break;
                default:
                    m_provided = context::thread_support::Multiple;
            }
        }
    }
    // ...............................................................................................
    context::~context( ) { if (has_initialized_parallel)  MPI_Finalize( ); }
    // ...............................................................................................
    communicator &context::globalCommunicator( ) {
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
                                                                            "XMPI_Output", globalCommunicator( ).rank));
        }
        return *pt_xmpi_logg;
    }
}
#endif

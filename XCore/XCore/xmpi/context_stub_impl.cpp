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

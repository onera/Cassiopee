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
#ifndef _XCORE_XMPI_CONTEXT_HPP_
#define _XCORE_XMPI_CONTEXT_HPP_

#include "Logger/logger.hpp"
#include "xmpi/communicator.hpp"
/*!  \namespace xcore
 *
 *   Namespace gathering the objects managing the parallel services
 */
namespace xcore
{
        /*!   \class context
     *    \brief Class which manages parallel context according to the
     *           choosen implementation
     *
     *    One instance of this class must be created. This instance
     *    initialize the underlying library choosen for the
     *    implementation ( MPI, PVM and so .). At the destruction of
     *    this unique object, the processus are synchronized and the
     *    parallel context is destroyed.
     *
     */
        class context
        {
        public:
            /*!
         *    \enum    thread_support
         *    \brief   Enumerate support levels for multithreads.
         *
         *    This enumeration describes the support level of the parallel
         *    library for multithreads. This level is defined by the
         *    choice of the library's user and by the support of the
         *    underlying library for multithreads.
         */
            enum thread_support {
                Single,     /*!< Only the main thread can make parallel calls */
                Funneled,   /*!< Only the thread that creates context will make parallel calls. */
                Serialized, /*!< Only one thread will make Parallel library calls at one time. */
                Multiple    /*!< Multiple threads may call MPI at once with no restrictions. */
            };
            /*!
            Default constructor. Only called if the MPI context initialization is done in third
            library ( mpi4py by example ).
            */
            context();
            /*!
         *    Construction initializing the parallel context
         *
         *    The constructor must be call only one time in the program
         *    ( in the beginning of the program by example ).
         *    This constructor provides two multithread supports :
         *    1. By default : The multiple thread support
         *    2. The Single multithread support
         *
         *     \param   nargc Number of arguments ( included the executable
         *                    name )
         *     \param   argv  Argument vector
         *     \param   isMultithreaded Boolean to choose multithread level support :
         *         True ( by default ) : Multiple thread level support
         *         False               : Single thread level support
         */
            context( int &nargc, char *argv[], bool isMultithreaded = true );
            /*!
         *     context constructor with fine control of multithread level support
         *
         *     \param   nargc  Number of arguments
         *     \param   argv   Argument vector
         *     \param   thread_level_support \link thread_support multithread level support
         */
            context( int &nargc, char *argv[], thread_support thread_level_support );

            /*!
         *    Destructor. Synchronize all processes and destroy the parallel context.
         */
            ~context();
            /*!
         *     Return the actual multithread level support
         */
            thread_support levelOfThreadSupport() const { return m_provided; }

            // The global communicator
            static communicator &globalCommunicator();
            static K_LOGGER::logger& logger();
        private:
            bool has_initialized_parallel; /*!< Has context initialize the parallel context ? */
            thread_support m_provided; /*!< Actual multithread level support */
            static communicator *pt_global_com;
        };
}  // namespace xcore
#endif

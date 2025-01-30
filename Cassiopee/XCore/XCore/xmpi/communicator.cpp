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

#include <cassert>

#include "xmpi/communicator.hpp"

#if defined( _MPI )
#include "xmpi/communicator_mpi_impl.hpp"
#else
#include "xmpi/communicator_stub_impl.hpp"
int xcore::communicator::Implementation::comm_world = 1;
#endif

namespace xcore
{
    communicator::communicator()
        : m_impl( new communicator::Implementation )
    {
        rank = m_impl->getRank();
        size = m_impl->getSize();
    }
    // .............................................................................
    communicator::communicator( const communicator &com, int color, int key )
        : m_impl( new communicator::Implementation( *com.m_impl, color, key ) )
    {
        rank = m_impl->getRank();
        size = m_impl->getSize();
    }
    // .............................................................................
    communicator::communicator( const communicator &com )
        : m_impl( new communicator::Implementation( *com.m_impl ) )
    {
        rank = m_impl->getRank();
        size = m_impl->getSize();
    }
    // .............................................................................
    communicator::communicator( const Communicator_ext_t& excom )
        : m_impl( new communicator::Implementation( excom ) )
    {
        rank = m_impl->getRank();
        size = m_impl->getSize();
    }
    // .............................................................................
    communicator::~communicator() { delete m_impl; }
    // .............................................................................
    Communicator_ext_t& communicator::get_implementation()
    {
        return m_impl->get_ext_comm();
    }
    // .............................................................................
    const Communicator_ext_t& communicator::get_implementation() const
    {
        return m_impl->get_ext_comm();
    }
    // =============================================================================
    int communicator::translateRank( const communicator &other_com ) const
    {
        int tr_rank;
        m_impl->translateRanks( *other_com.m_impl, 1, &rank, &tr_rank );
        return tr_rank;
    }
    // .............................................................................
    int communicator::translateRank( const communicator &other_com, int rk ) const
    {
        int tr_rank;
        m_impl->translateRanks( *other_com.m_impl, 1, &rk, &tr_rank );
        return tr_rank;
    }
    // .............................................................................
    std::vector<int> communicator::translateRanks( const communicator &other_com,
                                                   const std::vector<int> &ranksToTranslate )
    {
        std::vector<int> tr_ranks( ranksToTranslate.size() );
        m_impl->translateRanks( *other_com.m_impl, ranksToTranslate.size(), ranksToTranslate.data(),
                                tr_ranks.data() );
        return tr_ranks;
    }
    // =============================================================================
    void communicator::barrier() const { m_impl->barrier(); }
    // ========================================================================
    status communicator::probe( int source, int tag ) const { return m_impl->probe( source, tag ); }
    // ------------------------------------------------------------------------------------
    bool communicator::iprobe( status &status, int source, int tag ) const { return m_impl->iprobe( source, tag, status ); }
}  // namespace xcore

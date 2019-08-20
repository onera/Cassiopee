#include <cassert>

#include "xmpi/communicator.hpp"

#if defined( _MPI )
#include "xmpi/communicator_mpi_impl.hpp"
#else
#include "xmpi/communicator_stub_impl.hpp"
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
# if defined(_MPI)
    communicator::communicator( const MPI_Comm &excom )
        : m_impl( new communicator::Implementation( excom ) )
    {
        rank = m_impl->getRank();
        size = m_impl->getSize();
    }
# endif
    // .............................................................................
    communicator::~communicator() { delete m_impl; }
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
    status communicator::probe( int source, int tag ) { return m_impl->probe( source, tag ); }
    // ------------------------------------------------------------------------------------
    bool communicator::iprobe( status &status, int source, int tag ) { return m_impl->iprobe( source, tag, status ); }
}  // namespace xcore

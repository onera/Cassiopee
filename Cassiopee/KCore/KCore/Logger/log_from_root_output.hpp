# ifndef _KCORE_LOGGER_LOG_FROM_ROOT_OUTPUT_HPP_
#   define _KCORE_LOGGER_LOG_FROM_ROOT_OUTPUT_HPP_
#   include <cassert>
#   include "Memory/unique_ptr.hpp"
#   include "Logger/logger.hpp"

namespace K_LOGGER
{
  /*!
   * The output which display messages only for master root ( 0 by default )
   */
  template<typename L>
  class log_from_root_output : public logger::listener
  {
  public:
    /**
     * @brief      Logs output from exclusively root ( rank == 0 ).
     *
     * @param[in]  flags   The flags for the message to listen
     * @param[in]  lstner  The sequential listener to use for root. Must be created on fly.
     * @param[in]  rank    The rank of the process ( basically, if rank = 0, the process listen for messages)
     */
    log_from_root_output( int flags, const L* lstner, int rank ) :
      logger::listener(rank == 0 ? flags : K_LOGGER::logger::listener::listen_for_nothing),
      m_listener(lstner)
    { assert(lstner != NULL); }    
  private:
    virtual std::ostream& report() {
      return m_listener->report();
    }
    K_MEMORY::unique_ptr<L> m_listener;
  };
}

#endif
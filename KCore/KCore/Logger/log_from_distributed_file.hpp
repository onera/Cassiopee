#ifndef _KCORE_LOGGER_LOG_FROM_DISTRIBUTED_FILE_HPP_
# define _KCORE_LOGGER_LOG_FROM_DISTRIBUTED_FILE_HPP_
# include <fstream>
# include <cassert>
# include "Logger/logger.hpp"
# include "Logger/log_to_file.hpp"
# include "Memory/shared_ptr.hpp"

namespace K_LOGGER
{
  class log_from_distributed_file : public logger::listener
  {
  public:
    log_from_distributed_file( int flags, 
                               const std::string& basename, int rank );
    ~log_from_distributed_file();
    virtual std::ostream& report() {
      assert(m_pt_log.get() != NULL);
      return m_pt_log->report();
    }
  private:
    K_MEMORY::shared_ptr<log_to_file> m_pt_log;
  };
}


#endif

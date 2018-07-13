# include "Logger/log_from_distributed_file.hpp"
# include <sstream>
# include <iomanip>

namespace K_LOGGER
{
  log_from_distributed_file::log_from_distributed_file(int flags,
                         const std::string& basename, int rank ) try:
     logger::listener(flags),
     m_pt_log()
  {
        std::stringstream file_name;
        file_name << basename << std::setfill('0') << std::setw(5)
                    << rank << ".txt";
        K_MEMORY::shared_ptr<log_to_file> pt_log(new log_to_file(flags, file_name.str()));
      m_pt_log = pt_log;
  }
  catch(std::ios_base::failure) {
    std::cerr << "File creation failed. This listener will be unavailable." << std::endl;
  }
  // ...............................................................................................
  log_from_distributed_file::~log_from_distributed_file()
  {}
};

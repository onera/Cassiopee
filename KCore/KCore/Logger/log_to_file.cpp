#include "Logger/log_to_file.hpp"

namespace K_LOGGER {
    log_to_file::log_to_file( int flags, std::string const &filename, std::ios_base::openmode mode ) try : 
    logger::listener( flags ), m_fileName( filename ) {
            m_file.open( m_fileName.c_str(), mode );
        }
    catch ( std::ios_base::failure ) {
        std::cerr << "File creation failed. This listener will be unavailable."
                  << std::endl;
    }
    // ...............................................................................................
    log_to_file::~log_to_file( ) { m_file.close( ); }
    // ...............................................................................................
    std::ostream &log_to_file::report( ) { return m_file; }
}

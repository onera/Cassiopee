#include "Logger/log_to_std_error.hpp"

namespace K_LOGGER {
    log_to_std_error::log_to_std_error( int flags ) : logger::listener( flags ) {}
    // ...............................................................................................
    std::ostream &log_to_std_error::report( ) { return std::cerr; }
}

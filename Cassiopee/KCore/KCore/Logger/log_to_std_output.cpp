#include "Logger/log_to_std_output.hpp"

namespace K_LOGGER {
    log_to_std_output::log_to_std_output( int flags ) : logger::listener( flags ) {}
    // ...............................................................................................
    std::ostream &log_to_std_output::report( ) { return std::cout; }
}

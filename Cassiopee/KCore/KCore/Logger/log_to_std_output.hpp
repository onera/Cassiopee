#ifndef _KCORE_LOGGER_LOG_TO_STD_OUTPUT_HPP
#define _KCORE_LOGGER_LOG_TO_STD_OUTPUT_HPP
#include "Logger/logger.hpp"

namespace K_LOGGER {
    class log_to_std_output : public logger::listener {
    public:
        log_to_std_output( int flags );
        virtual std::ostream &report( );
    };
}

#endif

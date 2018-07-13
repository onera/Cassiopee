#ifndef _KCORE_LOGGER_LOG_TO_STD_ERROR_HPP
#define _KCORE_LOGGER_LOG_TO_STD_ERROR_HPP
#include "Logger/logger.hpp"

namespace K_LOGGER {
    class log_to_std_error : public logger::listener {
    public:
        log_to_std_error( int flags );

    private:
        virtual std::ostream &report( );
    };
}

#endif

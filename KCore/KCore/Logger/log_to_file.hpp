#ifndef _KCORE_LOGGER_LOG_TO_FILE_HPP
#define _KCORE_LOGGER_LOG_TO_FILE_HPP
#include "Logger/logger.hpp"
#include <fstream>

namespace K_LOGGER {
    /**
     * @brief      Class for print log in a file.
     */
    class log_to_file : public logger::listener {
    public:
        /**
         * @brief      Logs in a file.
         *
         * @param[in]  flags     The flags to know which kind of message this
         *                       listener must manage
         * @param      filename  The filename where write the messages
         */
        log_to_file( int flags, std::string const &filename, std::ios_base::openmode mode = std::ios_base::out );
        /**
         * @brief      Destroy the listener
         */
        ~log_to_file( );

        /**
         * @brief      Return the output stream managed by the listener
         *
         * @return     The output stream
         */
        virtual std::ostream &report( );

    private:
        std::string   m_fileName;
        std::ofstream m_file;
    };
}
#endif

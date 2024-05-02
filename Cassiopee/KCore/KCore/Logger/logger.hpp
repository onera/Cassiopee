#ifndef _KCORE_LOGGER_HPP_
#define _KCORE_LOGGER_HPP_
#include <iostream>
#include <list>
#include <sstream>
#include "Memory/unique_ptr.hpp"
#include "Memory/shared_ptr.hpp"

namespace K_LOGGER {
    /*!     \class  Logger
     *      \brief  This class manages a  logger which output some categorized
     * messages
     *              following listeners.
     *      \author Xavier JUVIGNY
     *
     *      A logger is a output system for different levels of messages (
     * Information, Warning, Debug,
     *      Trace, and so on. )
     *      By default, the logger has none listeners.
     *
     *      Exemple :
     *      --------
     *      K_LOGGER::logger log;
     *      // Add a listener for the console output for information only
     *      log.subscribe(new K_LOGGER::log_to_std_output(listener::listen_for_information));
     *      // Add a listener for the error output and assertion :
     *      log.subscribe(new K_LOGGER::log_file("ErrorLog", listener::listen_for_error |
     *                                                       listener::listen_for_assertion );
     *
     *      log_Information << "Bla bla blah" << std::endl;
     *      log_Error       << "Aïe Aïe Aïe"  << std::endl;
     */
    class logger {
    public:
        /*!    \class  Listener
         *     \brief  Base class for all listeners for the Logger ( as Console,
         *             File, Error output )
         *     \author Xavier JUVIGNY
         *
         *     A listener is a realisation of a stream output which writes only some type of messages.
         *     on a specific output. This class below is an abstract class. The concrete classes are
         *     defined in other files.
         *
         */
        class listener {
        public:
            /*! \enum Listener levels
             *
             *  Describe the level of the listener and the active channels for this listener.
             */
            enum {
                listen_for_nothing     = 0,
                listen_for_assertion   = 1,
                listen_for_error       = 2,
                listen_for_warning     = 4,
                listen_for_information = 8,
                listen_for_trace       = 16,

                listen_for_all = 0xFFFF
            };

            /*! \brief Construtor  : By default, build a listener who listen
             * none channel.
             */
            listener( int flags = listen_for_nothing );
            /*! \brief Destructor : Close the stream associated with the
             * listener.
             */
            virtual ~listener( );

            /*! \brief Return the ostream associated with the listener.
             */
            virtual std::ostream &report( ) = 0;
            /* ! \brief Test if the listener listen the channel mode given as
             * parameter.
             *
             *   \param mode  The channel mode ( Information, Debug, and so on.
             * )
             */
            bool toReport( int mo ) const { return m_flags & mo; }

        private:
            int m_flags;
        };

        class iterator;
        class const_iterator {
        public:
            const_iterator( const logger &log, bool end = false );
            const_iterator( const iterator &it );
            const_iterator( const const_iterator &it );
            ~const_iterator( );

            const_iterator &operator=( const const_iterator &it );
            bool operator!=( const const_iterator &it ) const;

            const listener *operator*( ) const;
            const_iterator &operator++( );

        private:
            struct Implementation;
            K_MEMORY::unique_ptr<Implementation> m_pt_impl;
        };

        class iterator {
        public:
            iterator( logger &log, bool end = false );
            iterator( const iterator &it );
            ~iterator( );

            iterator &operator=( const iterator &it );
            bool operator!=( const iterator &it ) const;

            listener *operator*( ) const;
            iterator &operator++( );
            friend class const_iterator;

        private:
            struct Implementation;
            K_MEMORY::unique_ptr<Implementation> m_pt_impl;
        };

        /*! \brief Constructor. Build a logger without listener.
             */
        logger( );

        /*! brief Destructor. Flush all ostream before deleting the logger
         */
        ~logger( );// { flush( ); }


        /*! \brief Add a new listerner for the log
         */
        bool subscribe( listener *listener );
        /*! \brief Remove a listener from the log.
         */
        bool unsubscribe( listener *listener );

        /*! \brief Change the channel mode of the messages.
         */
        logger &operator[]( int mo );

        logger &set_mode( int mo );
        int get_mode( ) const;

        /*! \brief Flush all ostream contained in the listeners.
         */
        logger &flush( ) {
            for ( logger::iterator it_listener = begin(); it_listener != end(); ++it_listener )
                (*it_listener)->report( ).flush( );
            return *this;
        }
        // ...............................................................................................
        iterator       begin( ) { return iterator( *this, false ); }
        const_iterator begin( ) const { return const_iterator( *this, false ); }
        iterator       end( ) { return iterator( *this, true ); }
        const_iterator end( ) const { return const_iterator( *this, true ); }
        // ..........................................................................
        /*! \brief Flux operator to print an object representation in each
         * listener
         *         who have the current channel mode available.
         */
        template <typename K>
        inline logger &operator<<( const K &obj ) {
            for ( logger::iterator it_listener = begin(); it_listener != end(); ++it_listener )
                if ( (*it_listener)->toReport( get_mode( ) ) ) { (*it_listener)->report( ) << obj; }
            return *this;
        }
        /*! \brief Channel modes
         */
        enum {
            nothing     = listener::listen_for_nothing,
            assertion   = listener::listen_for_assertion,
            error       = listener::listen_for_error,
            warning     = listener::listen_for_warning,
            information = listener::listen_for_information,
            trace       = listener::listen_for_trace,
            all         = listener::listen_for_all
        };

        /*!
         * Stream class to change the mode of output for the logger
         *
         * Example : logger << logger::mode(logger::information) << ...
         */
        struct mode {
        public:
            mode( int n ) : m_param( n ) {}

            int m_param;
        };

        /// \privatesection
        typedef logger &( *logger_manip )( logger & );
        logger &operator<<( logger_manip manip ) { return manip( *this ); }

        logger &operator<<( const mode &m ) {
            set_mode( m.m_param );
            return *this;
        }

        friend class iterator;
        friend class const_iterator;

        // Ansi codes
        static const char* Normal;
        static const char* Bright;
        static const char* Underline;
        static const char* Inverse;
        static const char* PrimaryFont;
        static const char* SecondFont;
        static const char* ThirdFont;
        static const char* FourthFont;
        static const char* FifthFont;
        static const char* SixthFont;
        static const char* SeventhFont;
        static const char* HeighthFont;
        static const char* NinthFont;
        static const char* TenthFont;
        static const char* NormalIntensity;
        static const char* NoUnderline;
        static const char* Black;
        static const char* Red;
        static const char* Green;
        static const char* Yellow;
        static const char* Blue;
        static const char* Magenta;
        static const char* Cyan;
        static const char* White;
        static const char* DefaultColor;
        static const char* BBlack;
        static const char* BRed;
        static const char* BGreen;
        static const char* BYellow;
        static const char* BBlue;
        static const char* BMagenta;
        static const char* BCyan;
        static const char* BWhite;
        static const char* DefaultBackground;
        static const char* Framed;
        static const char* Encircled;
        static const char* Overlined;
        static const char* NoFramed;

    private:
        logger( const logger &log );
        logger &operator=( const logger & );

        struct Implementation;
        K_MEMORY::shared_ptr<Implementation> m_pt_impl;
    };

#define LogAssert( cond )                                                                                          \
    K_LOGGER::logger::mode( ( ( cond ) ? K_LOGGER::logger::assertion : K_LOGGER::logger::nothing ) )                           \
        << "[\033[33m] [Assertion] " << std::string( __FILE__ ) << " in " << std::string( __FUNCTION__ ) << " at " \
        << __LINE__ << " ]\033[0m] : "

#define LogWarning                                                                                                     \
    K_LOGGER::logger::mode( K_LOGGER::logger::warning ) << "[ \033[31m[Warning] " << std::string( __FILE__ ) << " in "         \
                                                << std::string( __FUNCTION__ ) << " at " << __LINE__ \
                                                << " ]\033[0m] :"

#define LogError                                                                                                     \
    K_LOGGER::logger::mode( K_LOGGER::logger::error ) << "[\033[41;33m [Error] " << std::string( __FILE__ ) << " in "        \
                                              << std::string( __FUNCTION__ ) << " at " << __LINE__ \
                                              << " ]\033[0m] : "

#define LogInformation K_LOGGER::logger::mode( K_LOGGER::logger::information ) << "[\033[32;1mInformation\033[0m] "

#define LogTrace                                                                                                     \
    K_LOGGER::logger::mode( K_LOGGER::logger::trace ) << "[\033[32m [Trace] " << std::string( __FILE__ ) << " in "           \
                                              << std::string( __FUNCTION__ ) << " at " << __LINE__ \
                                              << " ]\033[0m : "
}

namespace std {
    inline K_LOGGER::logger &flush( K_LOGGER::logger &out ) {
        out.flush( );
        return out;
    }
    inline K_LOGGER::logger &endl( K_LOGGER::logger &out ) {
        out << "\n";
        out.flush( );
        return out;
    }
}

#endif

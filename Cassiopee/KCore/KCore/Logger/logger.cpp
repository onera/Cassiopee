#include "Logger/logger.hpp"
#include <algorithm>
#include <cassert>
#include <iostream>

namespace K_LOGGER {
    const char* logger::Normal            = "\033[0m";
    const char* logger::Bright            = "\033[1m";
    const char* logger::Underline         = "\033[4m";
    const char* logger::Inverse           = "\033[7m";
    const char* logger::PrimaryFont       = "\033[10m";
    const char* logger::SecondFont        = "\033[11m";
    const char* logger::ThirdFont         = "\033[12m";
    const char* logger::FourthFont        = "\033[13m";
    const char* logger::FifthFont         = "\033[14m";
    const char* logger::SixthFont         = "\033[15m";
    const char* logger::SeventhFont       = "\033[16m";
    const char* logger::HeighthFont       = "\033[17m";
    const char* logger::NinthFont         = "\033[18m";
    const char* logger::TenthFont         = "\033[18m";
    const char* logger::NormalIntensity   = "\033[22m";
    const char* logger::NoUnderline       = "\033[24m";
    const char* logger::Black             = "\033[30m";
    const char* logger::Red               = "\033[31m";
    const char* logger::Green             = "\033[32m";
    const char* logger::Yellow            = "\033[33m";
    const char* logger::Blue              = "\033[34m";
    const char* logger::Magenta           = "\033[35m";
    const char* logger::Cyan              = "\033[36m";
    const char* logger::White             = "\033[37m";
    const char* logger::DefaultColor      = "\033[49m";
    const char* logger::BBlack            = "\033[40m";
    const char* logger::BRed              = "\033[41m";
    const char* logger::BGreen            = "\033[42m";
    const char* logger::BYellow           = "\033[43m";
    const char* logger::BBlue             = "\033[44m";
    const char* logger::BMagenta          = "\033[45m";
    const char* logger::BCyan             = "\033[46m";
    const char* logger::BWhite            = "\033[47m";
    const char* logger::DefaultBackground = "\033[49m";
    const char* logger::Framed            = "\033[51m";
    const char* logger::Encircled         = "\033[52m";
    const char* logger::Overlined         = "\033[53m";
    const char* logger::NoFramed          = "\033[54m";
    // Global variable, invisible from linker
    /*namespace {
        std::mutex mutex_create, mutex_register, mutex_remove;
    }*/

    // logger Implementation :
    struct logger::Implementation {
        typedef std::list< listener* > Container;
        Implementation( ) : m_current_mode( logger::information ), m_listeners( ) {}
        int       m_current_mode;
        Container m_listeners;

        static K_MEMORY::shared_ptr< logger::Implementation > m_pt_shared_impl;
    };

    K_MEMORY::shared_ptr< logger::Implementation > logger::Implementation::m_pt_shared_impl((logger::Implementation*)NULL);
    // ...............................................................................................
    // logger iterator implementation
    struct logger::iterator::Implementation {
        typedef logger::Implementation::Container::iterator iterator_impl;
        iterator_impl                                       m_iterator;
    };
    // ...............................................................................................
    // logger const_iterator implementation
    struct logger::const_iterator::Implementation {
        typedef logger::Implementation::Container::const_iterator const_iterator_impl;
        const_iterator_impl                                       m_const_iterator;
    };
    // ===============================================================================================
    // logger iterator definition :
    logger::iterator::iterator( logger& log, bool end ) : m_pt_impl( new logger::iterator::Implementation ) {
        if ( end ) {
            m_pt_impl->m_iterator = log.m_pt_impl->m_listeners.end( );
        } else {
            m_pt_impl->m_iterator = log.m_pt_impl->m_listeners.begin( );
        }
    }
    // ...............................................................................................
    logger::iterator::iterator( const logger::iterator& it ) : m_pt_impl( new logger::iterator::Implementation ) {
        m_pt_impl->m_iterator = it.m_pt_impl->m_iterator;
    }
    // ...............................................................................................
    logger::iterator::~iterator( ) {}
    // ...............................................................................................
    logger::iterator& logger::iterator::operator=( const logger::iterator& it ) {
        if ( this != &it ) { m_pt_impl->m_iterator = it.m_pt_impl->m_iterator; }
        return *this;
    }
    // ...............................................................................................
    bool logger::iterator::operator!=( const logger::iterator& it ) const {
        return m_pt_impl->m_iterator != it.m_pt_impl->m_iterator;
    }
    // ...............................................................................................
    logger::listener* logger::iterator::operator*( ) const {
        return *m_pt_impl->m_iterator;
    }
    // ...............................................................................................
    logger::iterator& logger::iterator::operator++( ) {
        ( m_pt_impl->m_iterator )++;
        return *this;
    }
    // -----------------------------------------------------------------------------------------------
    logger::const_iterator::const_iterator( const logger& log, bool end )
        : m_pt_impl( new logger::const_iterator::Implementation ) {
        if ( end ) {
            m_pt_impl->m_const_iterator = log.m_pt_impl->m_listeners.end( );
        } else {
            m_pt_impl->m_const_iterator = log.m_pt_impl->m_listeners.begin( );
        }
    }
    // ...............................................................................................
    logger::const_iterator::const_iterator( const logger::iterator& it )
        : m_pt_impl( new logger::const_iterator::Implementation ) {
        m_pt_impl->m_const_iterator = it.m_pt_impl->m_iterator;
    }
    // ...............................................................................................
    logger::const_iterator::const_iterator( const logger::const_iterator& it )
        : m_pt_impl( new logger::const_iterator::Implementation ) {
        m_pt_impl->m_const_iterator = it.m_pt_impl->m_const_iterator;
    }
    // ...............................................................................................
    logger::const_iterator::~const_iterator( ) {}
    // ...............................................................................................
    logger::const_iterator& logger::const_iterator::operator=( const logger::const_iterator& it ) {
        if ( this != &it ) { m_pt_impl->m_const_iterator = it.m_pt_impl->m_const_iterator; }
        return *this;
    }
    // ...............................................................................................
    bool logger::const_iterator::operator!=( const logger::const_iterator& it ) const {
        return m_pt_impl->m_const_iterator != it.m_pt_impl->m_const_iterator;
    }
    // ...............................................................................................
    const logger::listener* logger::const_iterator::operator*( ) const {
        return *m_pt_impl->m_const_iterator;
    }
    // ...............................................................................................
    logger::const_iterator& logger::const_iterator::operator++( ) {
        ( m_pt_impl->m_const_iterator )++;
        return *this;
    }
    // -----------------------------------------------------------------------------------------------
    logger::logger( ) : m_pt_impl() {
        //std::lock_guard< std::mutex > lock( mutex_create );
        if ( logger::Implementation::m_pt_shared_impl == nullptr )
            logger::Implementation::m_pt_shared_impl = std::make_shared<logger::Implementation>();
        m_pt_impl                                    = logger::Implementation::m_pt_shared_impl;
    }
    // ...............................................................................................
    logger::~logger( ) { flush( ); }
    // ...............................................................................................
    bool logger::subscribe( logger::listener* listener ) {
        //std::lock_guard< std::mutex > lock( mutex_register );
        if ( listener == NULL ) return false;
        logger::Implementation::Container::iterator itL = std::find( m_pt_impl->m_listeners.begin( ), m_pt_impl->m_listeners.end( ), listener );
        if ( itL != m_pt_impl->m_listeners.end( ) ) return false;
        m_pt_impl->m_listeners.push_back( listener );
        return true;
    }
    // ...............................................................................................
    bool logger::unsubscribe( logger::listener* listener ) {
        //std::lock_guard< std::mutex > lock( mutex_remove );
        if ( listener == NULL ) return false;
        logger::Implementation::Container::iterator itL = std::find( m_pt_impl->m_listeners.begin( ), m_pt_impl->m_listeners.end( ), listener );
        if ( itL == m_pt_impl->m_listeners.end( ) ) return false;
        m_pt_impl->m_listeners.remove( listener );
        return true;
    }
    // ...............................................................................................
    logger& logger::operator[]( int mode ) {
        m_pt_impl->m_current_mode = mode;
        return *this;
    }
    // ...............................................................................................
    logger& logger::set_mode( int mode ) {
        m_pt_impl->m_current_mode = mode;
        return *this;
    }
    // ...............................................................................................
    int logger::get_mode( ) const {
        return m_pt_impl->m_current_mode;
    }
    // -----------------------------------------------------------------------------------------------
    logger::listener::listener( int flags ) : m_flags( flags ) {}
    // ...............................................................................................
    logger::listener::~listener( ) {}
}

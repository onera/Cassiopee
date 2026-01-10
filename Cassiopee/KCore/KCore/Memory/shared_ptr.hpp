/*
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#define CPPVERSION 2011

#ifndef _KCORE_SHARED_PTR_HPP_
#define _KCORE_SHARED_PTR_HPP_
#if  CPPVERSION >= 2011
# include <memory>
namespace K_MEMORY {
    template<typename K>
    using shared_ptr = std::shared_ptr<K>;
}
#else
namespace K_MEMORY
{
    template<typename K>
    class shared_ptr
    {
    public:
        typedef K element_type;

        shared_ptr() : m_pt_ref_counter(NULL), m_ptr(NULL)
        {}

        template<typename U>
        shared_ptr( U* pt ) : m_pt_ref_counter(new long(1)), m_ptr(static_cast<K*>(pt))
        {}

        template<typename U>
        shared_ptr( const shared_ptr<U>& pt ) : m_pt_ref_counter(pt.m_pt_ref_counter),
                                                m_ptr(static_cast<K*>(pt.m_ptr))
        {
            if ( m_pt_ref_counter )
                *m_pt_ref_counter += 1;
        }

        ~shared_ptr()
        {
            deallocate();
        }

        shared_ptr& operator = ( const shared_ptr& ptr )
        {
            if ( this != &ptr )
            {
                deallocate();
                this->m_pt_ref_counter = ptr.m_pt_ref_counter;
                this->m_ptr = ptr.m_ptr;
                if ( this->m_pt_ref_counter )
                    (*this->m_pt_ref_counter) += 1;
            }
            return *this;
        }

        template<typename U>
        shared_ptr& operator = ( const shared_ptr<U>& ptr )
        {
            deallocate();
            this->m_pt_ref_counter = ptr.m_pt_ref_counter;
            this->m_ptr = static_cast<K*>(ptr.m_ptr);
            if ( this->m_pt_ref_counter )
            (*this->m_pt_ref_counter) += 1;
            return *this;
        }

        void swap( shared_ptr& x )
        {
            long* i_tmp = this->m_pt_ref_counter;
            K*        x_tmp = this->m_ptr;
            this->m_pt_ref_counter = x.m_pt_ref_counter;
            this->m_ptr = x.m_ptr;
            x.m_pt_ref_counter = i_tmp;
            x.m_ptr = x_tmp;
        }

        void reset()
        {
            deallocate();
            this->m_pt_ref_counter = NULL;
            this->m_ptr            = NULL;
        }

        template<typename U>
        void reset( U* ptu )
        {
            deallocate();
            this->m_pt_ref_counter = new long(1);
            this->m_ptr = static_cast<K*>(ptu);
        }

        element_type* get() const
        {
            return this->m_ptr;
        }

        element_type& operator * () const
        {
            return *this->m_ptr;
        }

        element_type* operator ->() const
        {
            return this->m_ptr;
        }

        long use_count() const
        {
            if ( not this->m_pt_ref_counter ) return 0L;
            return *(this->m_pt_ref_counter);
        }

        bool unique() const
        {
            if ( (this->m_pt_ref_counter) && ( (*this->m_pt_ref_counter) == 1 ) )
                return true;
            return false;
        }

        operator bool() const
        {
            return this->m_ptr != NULL;
        }

        long* m_pt_ref_counter;
        K* m_ptr;
    private:
        void deallocate()
        {
            (*this->m_pt_ref_counter) -= 1;
            if ( 0 == *this->m_pt_ref_counter)
            {
                delete this->m_pt_ref_counter;
                delete this->m_ptr;
            }
        }
    };
}
#endif
#endif

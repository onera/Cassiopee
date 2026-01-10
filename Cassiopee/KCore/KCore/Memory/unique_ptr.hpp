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
#ifndef _KCORE_UNIQUE_PTR_HPP_
#define _KCORE_UNIQUE_PTR_HPP_
#if __cplusplus >= 201103L
# include <memory>
namespace K_MEMORY {
    template<typename K>
    using unique_ptr = std::unique_ptr<K>;
}
#else
#include <cassert>
#include <cstdlib>
#include <utility>
#include <algorithm>

namespace K_MEMORY {
    template<typename K>
    class unique_ptr {
    public:
        typedef K               value;
        typedef K&              reference;
        typedef const K&        const_reference;
        typedef K*              pointer;
        typedef const K*        const_pointer;

        unique_ptr() : m_pointer(NULL) {}
        unique_ptr( pointer p ) : m_pointer(p) {}
        unique_ptr( unique_ptr& p) : m_pointer(p.m_pointer) { p.m_pointer = NULL; }
        ~unique_ptr() { delete m_pointer; }

        unique_ptr& operator = ( unique_ptr& pt ) {
            if ( this != &pt ) {
                delete m_pointer;
                m_pointer = pt.m_pointer;
                pt.m_pointer = NULL;
            }
            return *this;
        }

        pointer release() { pointer pt = m_pointer; m_pointer = NULL; return pt; }
        void reset( pointer ptr ) { pointer old_pt = m_pointer; m_pointer = ptr; delete old_pt; }
        void swap( unique_ptr& ptr ) { std::swap(m_pointer, ptr.m_pointer); }

        pointer get() { return m_pointer; }
        const_pointer get() const { return m_pointer; }
        operator bool() const { return m_pointer != NULL; }
        reference operator * () { return *m_pointer; }
        const_reference operator * () const { return *m_pointer; }
        pointer operator ->() { return m_pointer; }
        const_pointer operator ->() const { return m_pointer; }
    private:
        pointer m_pointer;
    };
    template<typename K1, typename K2> bool
    operator == ( const unique_ptr<K1>& x, const unique_ptr<K2>& y ) { return x.get() == y.get(); }
    template<typename K1, typename K2> bool
    operator != ( const unique_ptr<K1>& x, const unique_ptr<K2>& y ) { return x.get() != y.get(); }
    template<typename K1, typename K2> bool
    operator < ( const unique_ptr<K1>& x, const unique_ptr<K2>& y ) { return x.get() < y.get(); }
    template<typename K1, typename K2> bool
    operator <= ( const unique_ptr<K1>& x, const unique_ptr<K2>& y ) { return x.get() <= y.get(); }
    template<typename K1, typename K2> bool
    operator > ( const unique_ptr<K1>& x, const unique_ptr<K2>& y ) { return x.get() > y.get(); }
    template<typename K1, typename K2> bool
    operator >= ( const unique_ptr<K1>& x, const unique_ptr<K2>& y ) { return x.get() >= y.get(); }

}
#endif
#endif

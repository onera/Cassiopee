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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_MEM_CHUNK_HXX
#define NUGA_MEM_CHUNK_HXX

#include "Nuga/include/allocator.hxx"

namespace NUGA
{

  ///
  template <typename T, typename allocator_t = NUGA::allocator<false>>
  struct mem_chunk
  {
    T*          data;
    E_Int       alloc_sz;

    mem_chunk(E_Int size) :data(allocator_t::template allocate<T>(size)), alloc_sz(size) {}
  };

}

#endif

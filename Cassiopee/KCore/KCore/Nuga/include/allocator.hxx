/*    
    Copyright 2013-2025 Onera.

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

#ifndef NUGA_ALLOCATOR_HXX
#define NUGA_ALLOCATOR_HXX

#include "Nuga/include/defs.h"
#include <iostream>

namespace NUGA
{

  template <int CALLOC>
  struct allocator
  {
    template <typename T> static T* allocate(E_Int size);
    template <typename T> static void deallocate(T*& ptr);
  };

  // CPP allocator
  template <>
  template <typename T>
  T* allocator<0>::allocate(E_Int size)
  {
    T* p{ nullptr };

    try
    {
      p = new T[size];// : (T*)malloc(size * sizeof(T));
      if (p == nullptr)
        throw "Memory error!!";
      else
        return p;
    }
    catch (T* s)
    {
      std::cout << "Memory problem in DynArray : " << s << std::endl;
    }

    return nullptr;
  }

  // C allocator
  template <>
  template <typename T>
  T* allocator<1>::allocate(E_Int size)
  {
    T* p{ nullptr };

    try
    {
      p = (T*)malloc(size * sizeof(T));
      if (p == nullptr)
        throw "Memory error!!";
      else
        return p;
    }
    catch (T* s)
    {
      std::cout << "Memory problem in DynArray : " << s << std::endl;
    }

    return nullptr;
  }

  // CPP deallocation
  template <>
  template <typename T>
  void allocator<0>::deallocate(T*& ptr)
  {
    if (ptr == nullptr) return;

    delete[] ptr;
    ptr = nullptr;
  }

  // C deallocation
  template <>
  template <typename T>
  void allocator<1>::deallocate(T*& ptr)
  {
    if (ptr == nullptr) return;

    free(ptr);
    ptr = nullptr;
  }

} //NUGA

#endif

/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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

/*



NUGA



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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
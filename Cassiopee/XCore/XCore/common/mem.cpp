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
#include "mem.h"

void *xmalloc(E_Int nbytes, const char *file, E_Int line)
{
    if (nbytes < 0) {
        fprintf(stderr, "trying to allocate a negative amount (%.2f) of bytes in file %s:%d\n",
            nbytes/1000000.f, file, line);
        abort();
    }
    void *ptr = malloc(nbytes);
    if (ptr == NULL && nbytes > 0) {
        fprintf(stderr, "\n    Failed to allocate %.2f MB in file %s, line " SF_D_ "\n\n",
            nbytes/1000000., file, line);
        abort();
    }
    return ptr;
}

void *xcalloc(E_Int count, E_Int size, const char *file, E_Int line)
{
    if (count < 0) {
        fprintf(stderr, "xcalloc: count (%d) is negative in file %s:%d\n", count, file, line);
        abort();
    }

    void *ptr = calloc(count, size);
    if (ptr == NULL && size > 0) {
        fprintf(stderr, "\n    Failed to allocate %.2f MB in file %s, line " SF_D_ "\n\n",
            count*size/1000000., file, line);
        abort();
    }
    return ptr;
}

void xfree(void *ptr, const char *file, E_Int line)
{
  if (ptr) {
    free(ptr);
    ptr = NULL;
  }
}

void *xresize(void *ptr, E_Int nbytes, const char *file, E_Int line)
{

    if (nbytes < 0) {
        fprintf(stderr, "trying to resize to a negative amount (%.2f) of bytes in file %s:%d\n",
            nbytes/1000000.f, file, line);
        abort();
    }

    ptr = realloc(ptr, nbytes);
    if (ptr == NULL && nbytes > 0) {
        fprintf(stderr, "\n    Failed to allocate %.2f MB in file %s, line " SF_D_ "\n\n",
            nbytes/1000000., file, line);
        abort();
    }
    return ptr;
}

Mem_arena::Mem_arena()
{
    arena = NULL;
    size = 0;
    capacity = 0;
    nresize = 0;
    nalloc = 0;
}

Mem_arena::Mem_arena(size_t nbytes)
{
    arena = XMALLOC(nbytes);
    size = 0;
    capacity = nbytes;
    nresize = 0;
    nalloc = 0;
}

void Mem_arena::reserve(size_t nbytes)
{
    assert(arena == NULL);
    arena = XMALLOC(nbytes);
    capacity = nbytes;
}

void *Mem_arena::alloc(size_t nbytes, const char *file, int line)
{
    void *ptr = NULL;
    
    nalloc++;
    
    if (size + nbytes > capacity) {
        /*
        nresize++;
        size_t new_size = size + nbytes;
        size_t new_capacity = (size_t)((E_Float)new_size * 1.5);
        arena = XRESIZE(arena, new_capacity);
        assert(size + nbytes < new_capacity);
        
        ptr = (void *)((char *)arena + size);
        
        size = new_size;
        capacity = new_capacity;
        */
        fprintf(stderr, "Mem_arena: out of memory at %s:%d\n", file, line);
        print_stats();
        exit(1);
    } else {
        ptr = (void *)((char *)arena + size);
        size += nbytes;
    }
    
    return ptr;
}

void Mem_arena::drop()
{
    XFREE(arena);
    size = capacity = nalloc = nresize = 0;
}

void Mem_arena::print_stats()
{
    printf("\n");
    printf("    /**** Mem_arena stats ****/\n");
    printf("      size:     %zu B\n", size);
    printf("      capacity: %zu B\n", capacity);
    printf("      usage:    %.2f %%\n", size * 100.0f / capacity);
    printf("      number of allocations:   %zu\n", nalloc);
    printf("      number of reallocations: %zu\n", nresize);
    printf("    /*************************/\n");
    printf("\n");
}
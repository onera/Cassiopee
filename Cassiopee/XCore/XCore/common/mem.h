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
#ifndef MEM_H
#define MEM_H

#include "xcore.h"
#include <stddef.h>

#define XMALLOC(nbytes) \
  xmalloc((nbytes), __FILE__, __LINE__)

#define XCALLOC(count, size) \
  xcalloc((count), (size), __FILE__, __LINE__)

#define XFREE(ptr) \
  (xfree((ptr), __FILE__, __LINE__), (ptr) = NULL)

#define XRESIZE(ptr, nbytes) \
  xresize((ptr), (nbytes), __FILE__, __LINE__)

void *xmalloc(E_Int, const char *, E_Int);
void *xcalloc(E_Int, E_Int, const char *, E_Int);
void xfree(void *, const char *, E_Int);
void *xresize(void *, E_Int, const char *, E_Int);

struct Mem_arena {
    void *arena;

    size_t size;
    size_t capacity;

    size_t nresize;
    size_t nalloc;

    Mem_arena();

    Mem_arena(size_t nbytes);

    void reserve(size_t nbytes);

    void *alloc(size_t nbytes, const char *file, int line);

    void drop();

    void print_stats();
};

#define ARENA_ALLOC(nbytes) arena.alloc((nbytes), __FILE__, __LINE__)

#endif

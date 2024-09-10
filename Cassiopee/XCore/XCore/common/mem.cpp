/*    
    Copyright 2013-2024 Onera.

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

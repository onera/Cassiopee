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

void *xmalloc(size_t , const char *, E_Int);
void *xcalloc(size_t, size_t, const char *, E_Int);
void xfree(void *, const char *, E_Int);
void *xresize(void *, size_t, const char *, E_Int);

#endif
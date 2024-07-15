#include "mem.h"

void *xmalloc(size_t nbytes, const char *file, E_Int line)
{
  assert(nbytes >= 0);
  void *ptr = malloc(nbytes);
  if (ptr == NULL && nbytes > 0) {
    fprintf(stderr, "\n    Failed to allocate %.2f MB in file %s, line " SF_D_ "\n\n",
      nbytes/1000000., file, line);
    assert(0);
    //abort();
  }
  return ptr;
}

void *xcalloc(size_t count, size_t size, const char *file, E_Int line)
{
  assert(count >= 0);
  assert(size >= 0);
  void *ptr = calloc(count, size);
  if (ptr == NULL && size > 0) {
    fprintf(stderr, "\n    Failed to allocate %.2f MB in file %s, line " SF_D_ "\n\n",
      count*size/1000000., file, line);
    assert(0);
    //abort();
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

void *xresize(void *ptr, size_t nbytes, const char *file, E_Int line)
{
  //assert(ptr);
  assert(nbytes >= 0);
  ptr = realloc(ptr, nbytes);
  if (ptr == NULL && nbytes > 0) {
    fprintf(stderr, "\n    Failed to allocate %.2f MB in file %s, line " SF_D_ "\n\n",
      nbytes/1000000., file, line);
    assert(0);
    //abort();
  }
  return ptr;
}

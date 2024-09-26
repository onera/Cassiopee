#include "Array.h"
#include "common/mem.h"

void ArrayI_free(ArrayI *arr)
{
    arr->count = 0;
    XFREE(arr->ptr);
}

void ArrayI_alloc(ArrayI *arr, E_Int nelem)
{
    arr->count = nelem;
    arr->ptr = (E_Int *)XMALLOC(nelem * sizeof(E_Int));
}
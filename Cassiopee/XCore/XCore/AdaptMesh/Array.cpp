#include "Array.h"
#include "common/mem.h"

void ArrayI_free(ArrayI *arr)
{
    XFREE(arr->ptr);
}
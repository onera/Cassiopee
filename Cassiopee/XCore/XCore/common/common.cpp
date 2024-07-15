#include "common.h"
#include <cstdio>

void merr(const char *fmt, ...)
{
	fprintf(stderr, "\n\t");
	va_list args;
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n\n");
}

void parray(Int *arr, Int n)
{
    for (Int i = 0; i < n; i++) printf("%d ", arr[i]);
    printf("\n");
}

/*
Int Get_pos(Int e, Int *pn, Int size)
{
    for (E_Int i = 0; i < size; i++) {
        if (pn[i] == e) return i;
    }
    
    return -1;
}

void Right_shift(E_Int *pn, E_Int pos, E_Int size)
{
    E_Int tmp[10];
    assert(size <= 10);
    for (E_Int i = 0; i < size; i++) tmp[i] = pn[i];
    for (E_Int i = 0; i < size; i++) pn[i] = tmp[(i+pos)%size];
}
*/
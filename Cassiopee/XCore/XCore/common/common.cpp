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
    for (Int i = 0; i < n; i++) printf(SF_D_ " ", arr[i]);
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
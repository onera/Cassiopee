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
#ifndef _COMMON_XCORE_H
#define _COMMON_XCORE_H

#include "float.h"
#include "xcore.h"

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

#define RAISE(error) \
    do { \
        char msg[1024] = {}; \
        strcat(msg, __func__); \
        strcat(msg, ": "); \
        strcat(msg, error); \
        PyErr_SetString(PyExc_ValueError, msg); \
    } while (0);

#define IntArray(n) (E_Int *)XCALLOC(n, sizeof(E_Int))
#define FloatArray(n) (E_Float *)XCALLOC(n, sizeof(E_Float))

#define E_FLOAT_MAX FLT_MAX
#define E_FLOAT_MIN -FLT_MAX

void merr(const char *fmt, ...);

void parray(E_Int *arr, E_Int n);

inline
E_Int Get_pos(E_Int e, E_Int *pn, E_Int size)
{
    for (E_Int i = 0; i < size; i++) {
        if (pn[i] == e) return i;
    }
    
    return -1;
}

inline
void Right_shift(E_Int *pn, E_Int pos, E_Int size)
{
    E_Int tmp[24];
    assert(size <= 24);
    for (E_Int i = 0; i < size; i++) tmp[i] = pn[i];
    for (E_Int i = 0; i < size; i++) pn[i] = tmp[(i+pos)%size];
}

#endif

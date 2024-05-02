#ifndef _COMMON_XCORE_H
#define _COMMON_XCORE_H

#include "float.h"
#include "xcore.h"

#define RAISE(error) PyErr_SetString(PyExc_ValueError, (error))

#define E_FLOAT_MAX FLT_MAX
#define E_FLOAT_MIN -FLT_MAX

E_Int Get_pos(E_Int e, E_Int *pn, E_Int size);
void Right_shift(E_Int *pn, E_Int pos, E_Int size);


#endif
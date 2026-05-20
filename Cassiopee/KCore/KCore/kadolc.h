// adolc definitions when using E_ADOUBLE
#ifdef E_ADOLC
#include <adolc/adolc.h>
#endif

#ifdef E_ADOUBLE
#define __AD(name) d__#name
#define E_Float adouble
#else
#define __AD(name) name
#endif

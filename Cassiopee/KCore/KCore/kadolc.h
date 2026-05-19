// adolc definitions when using E_ADOUBLE
#ifdef E_ADOLC
#include <adolc/adolc.h>
#endif

#ifdef E_ADOUBLE
#define AD__ d__
#define E_Float adouble
#else
#define AD__
#endif

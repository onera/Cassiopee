// To avoid pb when python is compiled with gcc and we compile with icc
#ifdef __INTEL_COMPILER
#define __attribute__(x) 
#endif
#include <math.h> // because of _hypot in python.h
#include "Python.h"

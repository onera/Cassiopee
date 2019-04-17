// To avoid pb when python is compiled with gcc and we compile with icc
#ifdef __INTEL_COMPILER
#define __attribute__(x) 
#endif
#include <math.h> // because of _hypot in python.h
#include "Python.h"
#if PY_VERSION_HEX >= 0x03000000
#define PyString_FromString PyBytes_FromString
#define PyString_AsString   PyBytes_AsString
#define PyString_Check PyBytes_Check
#define PyInt_Check PyLong_Check
#define PyInt_AsLong PyLong_AsLong
#define PyInt_FromLong PyLong_FromLong
#endif

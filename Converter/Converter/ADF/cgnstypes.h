#ifndef CGNSTYPES_H
#define CGNSTYPES_H

#define CG_BUILD_LEGACY 0
#define CG_BUILD_64BIT  0
#define CG_BUILD_SCOPE  0

#define CG_MAX_INT32 0x7FFFFFFF
#ifdef _WIN32
# define CG_LONG_T __int64
#else
// Modified Cassiopee (E_LONG)
#include "kcore.h"
# define CG_LONG_T E_LONG
#endif

#if CG_BUILD_LEGACY
# define CG_SIZEOF_SIZE    32
# define CG_SIZE_DATATYPE "I4"
# define cgerr_t  int
# define cgint_t  int
# define cgsize_t int
# define cgid_t   double
#else
# if CG_BUILD_64BIT
#  define CG_SIZEOF_SIZE    64
#  define CG_SIZE_DATATYPE "I8"
   typedef CG_LONG_T cgsize_t;
# else
#  define CG_SIZEOF_SIZE    32
#  define CG_SIZE_DATATYPE "I4"
   typedef int cgsize_t;
# endif
  typedef int cgerr_t;
  typedef int cgint_t;
  typedef double cgid_t;
#endif

typedef CG_LONG_T cglong_t;
typedef unsigned CG_LONG_T cgulong_t;

#endif

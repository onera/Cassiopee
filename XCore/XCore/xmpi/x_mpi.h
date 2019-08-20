#ifndef _XCORE_XMPI_X_MPI_H
#define _XCORE_XMPI_X_MPI_H

#  if defined( _MPI )
// Pour etre compatible avec Microsoft MPI sous MSys 2
#    if defined(_WIN64)
#    include <cstdint>
typedef __int64 int64_t;
#    endif
//
#  include <mpi.h>
#  endif

#endif
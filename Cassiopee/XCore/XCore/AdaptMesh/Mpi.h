#pragma once

#ifdef _MPI

#include <mpi.h>

#else

#define MPI_Comm_rank(a, b)
#define MPI_Comm_size(a, b)

#define MPI_Request char

#define MPI_Isend(a, b, c, d, e, f, g)
#define MPI_Irecv(a, b, c, d, e, f, g)
//#define MPI_Waitall(a) ((void)0)

#define MPI_Barrier(a)
#define MPI_Allreduce(a, b, c, d, e, f)
#define MPI_Waitall(a, b, c) ((void)0)

#define MPI_Scan(a, b, c, d, e, f)

#endif

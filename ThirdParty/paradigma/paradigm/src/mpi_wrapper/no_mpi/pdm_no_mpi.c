/*============================================================================
 * No MPI
 *============================================================================*/
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-value"
#endif
/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/


/*============================================================================
 * Definition des variables globales
 *============================================================================*/

/*============================================================================
 * Defintion des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Init(int *argc, char ***argv)
{
  argc;
  argv;

  return 0;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Finalize (void)
{
  return 0;
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm)
{
  pdm_mpi_comm;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_2_mpi_comm : Unavailable function with pdm_no_mpi library\n" );
  abort();
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_free_mpi_comm(void *pt_mpi_comm)
{
  pt_mpi_comm;

  return NULL;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_mpi_2_pdm_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm(void *pt_mpi_comm)
{
  pt_mpi_comm;

  return 0;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_File_open (wrapping de la fonction MPI_File_open)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_open(PDM_MPI_Comm comm, char *filename, int amode, PDM_MPI_File *fh)
{
  comm;
  filename;
  amode;
  fh;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_open : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_close (wrapping de la fonction MPI_File_close)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_close(PDM_MPI_File *fh)
{
  fh;

  return 0;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_seek (wrapping de la fonction MPI_File_seek)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_seek(PDM_MPI_File fh, PDM_MPI_Offset offset, int whence)
{
  fh;
  offset;
  whence;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_seek : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_size (wrapping de la fonction MPI_File_get_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_size(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  fh;
  offset;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_get_size : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_position (wrapping de la fonction MPI_File_get_position)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_position(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  fh;
  offset;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_get_position : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_set_view (wrapping de la fonction MPI_File_set_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_set_view(PDM_MPI_File fh, PDM_MPI_Offset disp, PDM_MPI_Datatype etype,
	              PDM_MPI_Datatype filetype, char *datarep)
{
  fh;
  disp;
  etype;
  filetype;
  datarep;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_set_view : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_view (wrapping de la fonction MPI_File_get_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_view(PDM_MPI_File fh, PDM_MPI_Offset *disp,
                      PDM_MPI_Datatype *etype, PDM_MPI_Datatype *filetype, char *datarep)
{
  fh;
  disp;
  etype;
  filetype;
  datarep;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_get_view : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at (wrapping de la fonction MPI_File_read_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                     int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  fh;
  offset;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_read_at : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at_all (wrapping de la fonction MPI_File_read_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  fh;
  offset;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_read_at_all : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at (wrapping de la fonction MPI_File_write_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                      int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  fh;
  offset;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_write_at : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at_all (wrapping de la fonction MPI_File_write_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  fh;
  offset;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_write_at_all : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read (wrapping de la fonction MPI_File_read)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read(PDM_MPI_File fh, void *buf, int count,
                  PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  fh;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_write_at : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_all (wrapping de la fonction MPI_File_read_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_all(PDM_MPI_File fh, void *buf, int count,
                      PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  fh;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_Read_all : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write (wrapping de la fonction MPI_File_write)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write(PDM_MPI_File fh, void *buf, int count,
                   PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  fh;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_write : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_all (wrapping de la fonction MPI_File_write_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_all(PDM_MPI_File fh, void *buf, int count,
                       PDM_MPI_Datatype datatype, int *n_octet_lus)

{
  fh;
  buf;
  count;
  datatype;
  n_octet_lus;

  PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_File_write_at_all : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcount;
  recvtype;
  root;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Gather : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Igather (wrapping de la fonction MPI_Igather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Igather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcount;
  recvtype;
  root;
  comm;
  request;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Igather : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gatherv (wrapping de la fonction MPI_Gatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int *recvcounts, int *displs,
                PDM_MPI_Datatype recvtype, int root, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcounts;
  displs;
  recvtype;
  root;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Gatherv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Recv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
             int tag, PDM_MPI_Comm comm)
{
  buf;
  count;
  datatype;
  source;
  tag;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Recv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Irecv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
              int tag, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  buf;
  count;
  datatype;
  source;
  tag;
  comm;
  request;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Irecv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Send (wrapping de la fonction MPI_Send)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Send(void *buf, int count, PDM_MPI_Datatype datatype, int dest,
             int tag, PDM_MPI_Comm comm)
{
  buf;
  count;
  datatype;
  dest;
  tag;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Send : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Issend (wrapping de la fonction MPI_Issend)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Issend(const void *buf, int count, PDM_MPI_Datatype datatype, int dest, int tag,
               PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  buf;
  count;
  datatype;
  dest;
  tag;
  comm;
  request;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Issend : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Wait (wrapping de la fonction MPI_Wait)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Wait(PDM_MPI_Request *request)

{
  request;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Wait : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_hindexed (wrapping de la fonction MPI_Type_hindexed)
 *
 *----------------------------------------------------------------------------*/


int PDM_MPI_Type_create_hindexed (int count,
                              const int array_of_blocklengths[],
                              const PDM_MPI_Aint array_of_displacements[],
                              PDM_MPI_Datatype oldtype,
                              PDM_MPI_Datatype *newtype)
{
  count;
  array_of_blocklengths;
  array_of_displacements;
  oldtype;
  newtype;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Type_create_hindexed : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_commit (wrapping de la fonction MPI_Type_commit)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_commit(PDM_MPI_Datatype *datatype)
{
  datatype;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Type_commit : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_free (wrapping de la fonction MPI_Type_free)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_free(PDM_MPI_Datatype *datatype)
{
  datatype;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Type_free : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_f2c (wrapping de la fonction MPI_comm_f2c)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_Comm_f2c(PDM_MPI_Fint comm)
{
  comm;

  return 0;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_c2f (wrapping de la fonction MPI_comm_c2f)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Fint PDM_MPI_Comm_c2f(PDM_MPI_Comm comm)
{
  comm;

  return 0;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Scatter (wrapping de la fonction MPI_Scatter)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scatter(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
                int root, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcount;
  recvtype;
  root;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Scatter : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Barrier (wrapping de la fonction MPI_Barrier)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Barrier(PDM_MPI_Comm comm)
{
  comm;

  return 0;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Bcast (wrapping de la fonction MPI_Bcast)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Bcast(void *buffer, int count, PDM_MPI_Datatype datatype,
              int root, PDM_MPI_Comm comm)
{
  buffer;
  count;
  datatype;
  root;
  comm;

  return 0;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                  void *recvbuf, int recvcount,
                  PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcount;
  recvtype;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Allgather : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgatherv (wrapping de la fonction MPI_Allgatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts,
                   int *displs, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcounts;
  displs;
  recvtype;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Allgatherv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce (wrapping de la fonction MPI_Reduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		   PDM_MPI_Datatype datatype, PDM_MPI_Op op,
		   int root, PDM_MPI_Comm comm)
{
  sendbuf;
  recvbuf;
  count;
  datatype;
  op;
  root;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Reduce : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allreduce (wrapping de la fonction MPI_Allreduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm)
{
  sendbuf;
  recvbuf;
  count;
  datatype;
  op;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Allreduce : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Scan (wrapping de la fonction MPI_Scan)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm)
{
  sendbuf;
  recvbuf;
  count;
  datatype;
  op;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Scan : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}


int PDM_MPI_Iscan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm,
             PDM_MPI_Request *request)
{
  sendbuf;
  recvbuf;
  count;
  datatype;
  op;
  comm;
  request;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Iscan : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoall (wrapping de la fonction MPI_Alltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                 void *recvbuf, int recvcount,
                 PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcount;
  recvtype;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Alltoall : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Ialltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                 void *recvbuf, int recvcount,
                 PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm,
                 PDM_MPI_Request *request)
{
  sendbuf;
  sendcount;
  sendtype;
  recvbuf;
  recvcount;
  recvtype;
  comm;
  request;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Ialltoall : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoallv (wrapping de la fonction MPI_Alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                  PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                  int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcounts;
  sdispls;
  rdispls;
  sendtype;
  recvbuf;
  recvcounts;
  recvtype;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Alltoallv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

int PDM_MPI_Alltoallv_l(void *sendbuf, int *sendcounts, size_t *sdispls,
                      PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                      size_t *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  sendbuf;
  sendcounts;
  sdispls;
  rdispls;
  sendtype;
  recvbuf;
  recvcounts;
  recvtype;
  comm;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Alltoallv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoallv (wrapping de la fonction MPI_Ialltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                  PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                  int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm,
                  PDM_MPI_Request *request)
{
  sendbuf;
  sendcounts;
  sdispls;
  rdispls;
  sendtype;
  recvbuf;
  recvcounts;
  recvtype;
  comm;
  request;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Ialltoallv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Error_string (wrapping de la fonction MPI_Error_string)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Error_string(int errorcode, char *string, int *resultlen)
{
  errorcode;
  string;
  resultlen;

  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_Allgatherv : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_rank(PDM_MPI_Comm comm, int *rank)
{
  comm;
  *rank=0;
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_size(PDM_MPI_Comm comm, int *size)
{
  comm;
  *size=1;
  return 1;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_get_max_error_string
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_get_max_error_string(void)
{
  PDM_error(__FILE__, __LINE__, 0, "PDM_MPI_get_max_error_string : Unavailable function with pdm_no_mpi library\n" );
  abort();
  return 1;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_free
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_free(PDM_MPI_Comm *comm)
{
  comm;

  return 0;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_split(PDM_MPI_Comm comm, int color, int key, PDM_MPI_Comm *newcomm)
{
  comm;
  color;
  key;

  *newcomm = PDM_MPI_COMM_NULL;

  return 0;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

#if defined(__clang__)
#pragma clang diagnostic pop
#endif

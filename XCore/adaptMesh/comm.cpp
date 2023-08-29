#include "proto.h"

void comm_waitall(mesh *M)
{
  MPI_Waitall(M->nreq, M->req, MPI_STATUSES_IGNORE);
  M->nreq = 0;
}

void comm_interface_data_d(mesh *M, E_Float *data, E_Int stride, E_Float **recv_data)
{
  E_Int *owner = M->owner;
  E_Int nppatches = M->nppatches;
  E_Float *pt; 
  for (E_Int i = 0; i < nppatches; i++) {
    E_Int nfaces = M->ppatches[i].nfaces;
    E_Int *faces = M->ppatches[i].faces;
    E_Int dest = M->ppatches[i].nei_proc;

    M->ppatches[i].send_buf_d = (E_Float *)realloc(
      M->ppatches[i].send_buf_d, stride*nfaces * sizeof(E_Float));

    E_Int l = 0; 
    for (E_Int j = 0; j < nfaces; j++) {
      assert(faces[j] < M->nfaces);
      E_Int own = owner[faces[j]];
      pt = &data[stride*own]; 
      for (E_Int k = 0; k < stride; k++)
        M->ppatches[i].send_buf_d[l++] = pt[k];
    }
    assert(l == stride*nfaces);

    MPI_Irecv(recv_data[i], stride*nfaces, MPI_DOUBLE, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);
    MPI_Isend(M->ppatches[i].send_buf_d, stride*nfaces, MPI_DOUBLE, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);

    assert(M->nreq < 2*M->npc);
  }

  comm_waitall(M);
}

void comm_interface_data_i(mesh *M, E_Int *data, E_Int stride, E_Int **recv_data)
{
  E_Int *owner = M->owner;
  E_Int nppatches = M->nppatches;
  E_Int *pt; 
  for (E_Int i = 0; i < nppatches; i++) {
    E_Int nfaces = M->ppatches[i].nfaces;
    E_Int *faces = M->ppatches[i].faces;
    E_Int dest = M->ppatches[i].nei_proc;

    M->ppatches[i].send_buf_i = (E_Int *)realloc(
      M->ppatches[i].send_buf_i, stride*nfaces * sizeof(E_Int));

    E_Int l = 0; 
    for (E_Int j = 0; j < nfaces; j++) {
      assert(faces[j] < M->nfaces);
      E_Int own = owner[faces[j]];
      pt = &data[stride*own]; 
      for (E_Int k = 0; k < stride; k++)
        M->ppatches[i].send_buf_i[l++] = pt[k];
    }
    assert(l == stride*nfaces);

    MPI_Irecv(recv_data[i], stride*nfaces, MPI_INT, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);
    MPI_Isend(M->ppatches[i].send_buf_i, stride*nfaces, MPI_INT, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);

    assert(M->nreq < 2*M->npc);
  }

  comm_waitall(M);
}

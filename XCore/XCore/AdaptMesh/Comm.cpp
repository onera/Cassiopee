#include "Proto.h"
#include "../common/mem.h"

static
void comm_waitall(AMesh *M)
{
  MPI_Waitall(M->nreq, M->req, MPI_STATUSES_IGNORE);
  M->nreq = 0;
}

void exchange_proc_data_d(AMesh *M, E_Float *data, E_Float ***rbuf)
{
  // Allocate
  if (*rbuf == NULL) {
    *rbuf = (E_Float **)XCALLOC(M->npatches, sizeof(E_Float *));
    
    for (E_Int i = 0; i < M->npatches; i++) {
      E_Int nf = M->patches[i].nfaces;
      
      (*rbuf)[i] = (E_Float *)XCALLOC(nf, sizeof(E_Float));
    }
  }

  for (E_Int i = 0; i < M->npatches; i++) {
    E_Int nf = M->patches[i].nfaces;
    E_Int *faces = M->patches[i].faces;
 
    E_Float *ptr = M->patches[i].sbuf_d;
    assert(ptr != NULL);

    for (E_Int j = 0; j < nf; j++) {
      E_Int cell = M->owner[faces[j]];
      *ptr++ = data[cell];
    }
    
    MPI_Irecv((*rbuf)[i], nf, MPI_DOUBLE, M->patches[i].nei_proc, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);
    MPI_Isend(M->patches[i].sbuf_d, nf, MPI_DOUBLE, M->patches[i].nei_proc, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);

    assert(M->nreq < 2*M->npc);
  }

  comm_waitall(M);
}
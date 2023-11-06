/*    
    Copyright 2013-2023 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "proto.h"

void comm_waitall(mesh *M)
{
  MPI_Waitall(M->nreq, M->req, MPI_STATUSES_IGNORE);
  M->nreq = 0;
}

void comm_interface_data_d(mesh *M, E_Float *data, E_Int stride,
  E_Float **recv_data)
{
  E_Int *owner = M->owner;
  E_Int nppatches = M->nppatches;
  E_Float *pt; 
  for (E_Int i = 0; i < nppatches; i++) {
    E_Int nfaces = M->ppatches[i].nfaces;
    E_Int *faces = M->ppatches[i].faces;
    E_Int dest = M->ppatches[i].nei_proc;

    M->ppatches[i].send_buf_d = (E_Float *)XRESIZE(
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

void comm_interface_data_i(mesh *M, E_Int *data, E_Int stride,
  E_Int **recv_data)
{
  E_Int *owner = M->owner;
  E_Int nppatches = M->nppatches;
  E_Int *pt; 
  for (E_Int i = 0; i < nppatches; i++) {
    E_Int nfaces = M->ppatches[i].nfaces;
    E_Int *faces = M->ppatches[i].faces;
    E_Int dest = M->ppatches[i].nei_proc;

    M->ppatches[i].send_buf_i = (E_Int *)XRESIZE(
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

    MPI_Irecv(recv_data[i], (int)stride*nfaces, XMPI_INT, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);
    MPI_Isend(M->ppatches[i].send_buf_i, (int)stride*nfaces, XMPI_INT, dest, 0,
      MPI_COMM_WORLD, &M->req[M->nreq++]);

    assert(M->nreq < 2*M->npc);
  }

  comm_waitall(M);
}
/*    
    Copyright 2013-2024 Onera.

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
#include "Mesh.h"
#include "common/mem.h"

PyObject *K_XCORE::AdaptMesh_AssignRefData(PyObject *self, PyObject *args)
{
    PyObject *MESH, *CREF;

    if (!PYPARSETUPLE_(args, OO_, &MESH, &CREF)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    Int *ptr = NULL;
    Int ret, nfld, size;
    ret = K_NUMPY::getFromNumpyArray(CREF, ptr, size, nfld, true);
    if (ret != 1 || size != M->nc || nfld != 1) {
        RAISE("Bad cref input.");
        return NULL;
    }

    M->cref = (Int *)XRESIZE(M->cref, M->nc * sizeof(Int));
    for (Int i = 0; i < M->nc; i++) M->cref[i] = ptr[i];

    Py_DECREF(CREF);

    // Allocate patch buffers

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];
        P->sbuf_i = (Int *)XRESIZE(P->sbuf_i, P->nf * sizeof(Int));
        P->rbuf_i = (Int *)XRESIZE(P->rbuf_i, P->nf * sizeof(Int));
        memset(P->sbuf_i, 0, P->nf * sizeof(Int));
    }

    // Smooth cell ref data

    if (M->pid == 0) puts("Smoothing cell refinement data...");

    Mesh_smooth_cref(M);

    if (M->pid == 0) puts("Assigning face refinement data...");

    // Allocate face ref data

    M->fref = (Int *)XRESIZE(M->fref, M->nf * sizeof(Int));
    memset(M->fref, 0, M->nf * sizeof(Int));

    // Do patch faces first

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];
        memset(P->sbuf_i, 0, P->nf * sizeof(Int));
    }

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int face = P->pf[j];
            Int own = M->owner[face];

            if (M->cref[own] > 0) {
                Int flvl = M->flevel[face];
                Int clvl = M->clevel[own];
                assert(flvl >= clvl);

                if (flvl == clvl) {
                    M->fref[face] = 1;
                    P->sbuf_i[j] = 1;
                }
            }
        }

        MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        MPI_Irecv(P->rbuf_i, P->nf, XMPI_INT, P->nei, P->nei,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
    }

    Mesh_comm_waitall(M);

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int rval = P->rbuf_i[j];
            if (rval > 0) {
                Int face = P->pf[j];
                //Int own = M->owner[face];
                //Int flvl = M->flevel[face];
                //Int clvl = M->clevel[own];
                //if (flvl == clvl) M->fref[face] = 1;
                M->fref[face] = 1;
            }
        }
    }

    // Do the ref_cells faces
    for (Int cid = 0; cid < M->nc; cid++) {
        if (M->cref[cid] == 0) continue;

        Int *cell = Mesh_get_cell(M, cid);
        Int *crange = Mesh_get_crange(M, cid);
        Int cstride = M->cstride[cid];
        Int clvl = M->clevel[cid];

        for (Int i = 0; i < cstride; i++) {
            Int *pf = cell + 4*i;

            for (Int j = 0; j < crange[i]; j++) {
                Int face = pf[j];
                if (M->fref[face] == 1) continue;

                Int flvl = M->flevel[face];

                assert(flvl >= clvl);

                if (flvl == clvl) M->fref[face] = 1;
            }
        }
    }

    return Py_None;
}

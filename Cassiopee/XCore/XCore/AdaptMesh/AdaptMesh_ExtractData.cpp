/*    
    Copyright 2013-2026 ONERA.

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

PyObject *K_XCORE::AdaptMesh_ExtractOwners(PyObject *self, PyObject *args)
{
    PyObject *MESH;

    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    if (M->pid == 0) puts("Extracting face owners...");

    npy_intp dims[2];

    dims[1] = 1;
    dims[0] = (npy_intp)M->nf;
    PyArrayObject *OWN = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);

    E_Int *po = (E_Int *)PyArray_DATA(OWN);

    for (E_Int i = 0; i < M->nf; i++) {
        po[i] = M->owner[i];
    }

    return (PyObject *)OWN;
}

PyObject *K_XCORE::AdaptMesh_ExtractNeighbours(PyObject *self, PyObject *args)
{
    PyObject *MESH;

    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    if (M->pid == 0) puts("Extracting face neighbours...");

    npy_intp dims[2];

    dims[1] = 1;
    dims[0] = (npy_intp)M->nf;
    PyArrayObject *NEI = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);

    E_Int *pn = (E_Int *)PyArray_DATA(NEI);

    for (E_Int i = 0; i < M->nf; i++) {
        pn[i] = M->neigh[i];
    }

    return (PyObject *)NEI;
}

PyObject *K_XCORE::AdaptMesh_ExtractCellLevels(PyObject *self, PyObject *args)
{
    PyObject *MESH;

    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    if (M->pid == 0) puts("Extracting cell levels...");

    npy_intp dims[2];

    dims[1] = 1;
    dims[0] = (npy_intp)M->nc;
    PyArrayObject *LVL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);

    E_Int *pl = (E_Int *)PyArray_DATA(LVL);

    for (E_Int i = 0; i < M->nc; i++) {
        pl[i] = M->clevel[i];
    }

    return (PyObject *)LVL;
}


PyObject *K_XCORE::AdaptMesh_ExtractCellRanges(PyObject *self, PyObject *args)
{
    PyObject *MESH;

    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    if (M->pid == 0) puts("Extracting cell ranges...");

    PyObject *STR = K_NUMPY::buildNumpyArray(M->nc, 6, 1, 1);

    E_Int *ptr = K_NUMPY::getNumpyPtrI(STR);

    for (E_Int cid = 0; cid < M->nc; cid++) {
        E_Int *crange = Mesh_get_crange(M, cid);

        for (E_Int j = 0; j < M->cstride[cid]; j++) {
            ptr[cid + j*M->nc] = crange[j];
        }
    }

    return (PyObject *)STR;
}

PyObject *K_XCORE::AdaptMesh_ExtractHaloCellLevels(PyObject *self, PyObject *args)
{
    PyObject *MESH;

    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    if (M->pid == 0) puts("Extracting halo cell levels...");

    // Allocate

    PyObject *out = PyList_New(0);

    npy_intp dims[2];

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        P->sbuf_i = (E_Int *)XRESIZE(P->sbuf_i, P->nf * sizeof(E_Int));

        dims[1] = 1;
        dims[0] = (npy_intp)P->nf;

        PyArrayObject *LVL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);

        PyList_Append(out, (PyObject *)LVL);
        Py_DECREF(LVL);
    }

    // Exchange

    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        PyObject *LVL = PyList_GetItem(out, i);

        E_Int *pl = (E_Int *)PyArray_DATA((PyArrayObject *)LVL);

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int fid = P->pf[j];
            E_Int own = M->owner[fid];
            P->sbuf_i[j] = M->clevel[own];
        }

        MPI_Isend(P->sbuf_i, P->nf, XMPI_INT, P->nei, M->pid,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
        MPI_Irecv(pl, P->nf, XMPI_INT, P->nei, P->nei,
            MPI_COMM_WORLD, &M->reqs[M->nrq++]);
    }

    Mesh_comm_waitall(M);

    return out;
}
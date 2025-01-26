#include "Mesh.h"

PyObject *K_XCORE::AdaptMesh_ExtractTaggedFaces(PyObject *self, PyObject *args)
{
    PyObject *AMESH;
    if (!PYPARSETUPLE_(args, O_, &AMESH)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(AMESH, "AdaptMesh")) {
        RAISE("Bad AdaptMesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(AMESH, "AdaptMesh");

    assert(M->ftag);

    E_Int count = 0;

    for (E_Int i = 0; i < M->nf; i++) {
        count += (M->ftag[i] == 1);
    }

    npy_intp dims[2];
    dims[1] = 1;
    dims[0] = (npy_intp)count;
    
    PyArrayObject *out = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *pf = (E_Int *)PyArray_DATA(out);
    E_Int *ptr = pf;

    for (E_Int fid = 0; fid < M->nf; fid++) {
        if (M->ftag[fid] == 1)
            *ptr++ = fid;
    }
    
    return (PyObject *)out;
}

PyObject *K_XCORE::AdaptMesh_TagFaces(PyObject *self, PyObject *args)
{
    PyObject *AMESH, *FACES;
    if (!PYPARSETUPLE_(args, OO_, &AMESH, &FACES)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(AMESH, "AdaptMesh")) {
        RAISE("Bad AdaptMesh hook.");
        return NULL;
    }

    E_Int *faces;
    E_Int NF;
    E_Int ret = K_NUMPY::getFromNumpyArray(FACES, faces, NF, true);
    if (ret != 1) {
        RAISE("Bad face list.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(AMESH, "AdaptMesh");

    assert(M->ftag);
    memset(M->ftag, 0, M->nf * sizeof(E_Int));
    
    for (E_Int i = 0; i < NF; i++) {
        M->ftag[faces[i]] = 1;
    }

    return Py_None;
}

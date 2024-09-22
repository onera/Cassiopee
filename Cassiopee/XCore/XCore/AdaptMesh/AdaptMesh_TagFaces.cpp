#include "Mesh.h"

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

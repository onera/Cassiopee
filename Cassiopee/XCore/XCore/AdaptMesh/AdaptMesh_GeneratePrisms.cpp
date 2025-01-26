#include "Mesh.h"

PyObject *K_XCORE::AdaptMesh_GeneratePrisms(PyObject *self, PyObject *args)
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
    E_Int nf;
    E_Int ret = K_NUMPY::getFromNumpyArray(FACES, faces, nf, true);
    if (ret != 1) {
        RAISE("Bad face list.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(AMESH, "AdaptMesh");

    Mesh_generate_prisms(M, faces, nf);

    return Py_None;
}

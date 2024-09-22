#include "Mesh.h"

PyObject *K_XCORE::AdaptMesh_TriangulateFaces(PyObject *self, PyObject *args)
{
    PyObject *AMESH, *FACES;

    if (!PYPARSETUPLE_(args, OO_, &AMESH, &FACES)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(AMESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
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

    Mesh_triangulate_faces(M, faces, NF);

    return Py_None;
}

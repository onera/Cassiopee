#include "Mesh.h"

PyObject *K_XCORE::AdaptMesh_ExtractMesh(PyObject *self, PyObject *args)
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

    if (M->pid == 0) puts("Extracting mesh...");

    PyObject *karray = Mesh_export_karray(M);

    return karray;
}

#include "Mesh.h"

PyObject *K_XCORE::AdaptMesh_LoadBalance(PyObject *self, PyObject *args)
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

    // Make global cell-cell connectivity
    Mesh_make_cell_cells(M);

    // Load balance
    Int ret = Mesh_load_balance(M);

    if (ret != 0) {
        RAISE("Failed to load-balance the mesh.");
        return NULL;
    }

    if (M->pid == 0) puts("OK LoadBalance");

    return Py_None;
}

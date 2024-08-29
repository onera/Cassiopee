#include "mesh.h"

PyObject *K_XCORE::IntersectMesh_ExtractMesh(PyObject *self, PyObject *args)
{
    PyObject *MESH;
  
    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "IntersectMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    IMesh *M = (IMesh *)PyCapsule_GetPointer(MESH, "IntersectMesh");

    puts("Extracting mesh...");

    PyObject *master_out = M->export_karray();

    return master_out;
}
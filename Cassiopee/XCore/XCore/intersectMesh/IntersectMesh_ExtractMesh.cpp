#include "mesh.h"

PyObject *K_XCORE::IntersectMesh_ExtractMesh(PyObject *self, PyObject *args)
{
    PyObject *MESH;
    E_Int remove_periodic;
  
    if (!PYPARSETUPLE_(args, O_ I_, &MESH, &remove_periodic)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "IntersectMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    IMesh *M = (IMesh *)PyCapsule_GetPointer(MESH, "IntersectMesh");

    puts("Extracting mesh...");

    PyObject *master_out = M->export_karray(remove_periodic);

    return master_out;
}
#include "mesh.h"

PyObject *K_XCORE::IntersectMesh_Exit(PyObject *self, PyObject *args)
{
    PyObject *MASTER;
  
    if (!PYPARSETUPLE_(args, O_, &MASTER)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MASTER, "IntersectMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    IMesh *M = (IMesh *)PyCapsule_GetPointer(MASTER, "IntersectMesh");

    delete M;

    return Py_None;
}
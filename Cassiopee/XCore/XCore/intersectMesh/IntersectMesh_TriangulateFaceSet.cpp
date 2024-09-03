#include "mesh.h"

PyObject *K_XCORE::IntersectMesh_TriangulateFaceSet(PyObject *self,
    PyObject *args)
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

    puts("Triangulating master projection faces...");

    M->triangulate_face_set();

    return Py_None;
}
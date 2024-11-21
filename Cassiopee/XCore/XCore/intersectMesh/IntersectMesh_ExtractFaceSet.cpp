#include "mesh.h"

PyObject *K_XCORE::IntersectMesh_ExtractFaceSet(PyObject *self, PyObject *args)
{
    PyObject *IMESH;
    if (!PYPARSETUPLE_(args, O_, &IMESH)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(IMESH, "IntersectMesh")) {
        RAISE("Bad IntersectMesh hook.");
        return NULL;
    }

    IMesh *M = (IMesh *)PyCapsule_GetPointer(IMESH, "IntersectMesh");
    
    // Allocate
    npy_intp dims[2];
    dims[1] = 1;
    dims[0] = (npy_intp)M->faces_to_tri.size();
    
    PyArrayObject *out = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *pf = (E_Int *)PyArray_DATA(out);
    E_Int *ptr = pf;

    for (const E_Int fid : M->faces_to_tri)
        *ptr++ = fid;

    return (PyObject *)out;
}

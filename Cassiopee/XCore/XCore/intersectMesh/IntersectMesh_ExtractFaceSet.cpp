/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
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

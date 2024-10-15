/*    
    Copyright 2013-2024 Onera.

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
#include "xcore.h"
#include "common/Karray.h"
#include "common/common.h"
#include "mesh.h"
#include "smesh.h"
#include "dcel.h"
#include "vertex.h"
#include "face.h"
#include "hedge.h"
#include "io.h"
#include "cycle.h"
#include "triangle.h"

PyObject *K_XCORE::prepareMeshesForIntersection(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVE, *TAG;
  
    if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVE, &TAG)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MASTER, "IntersectMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    IMesh &M = *(IMesh *)PyCapsule_GetPointer(MASTER, "IntersectMesh");

    M.make_skin();
    assert(M.patch.empty());
    for (E_Int fid : M.skin) {
        assert(fid < M.nf);
        M.patch.insert(fid);
    }

    Karray sarray;

    E_Int ret;

    ret = Karray_parse_ngon(SLAVE, sarray);

    if (ret != 0) {
        RAISE("Bad slave mesh.");
        return NULL;
    }

    puts("Preparing meshes for intersection...");

    // Check slave point tags
    E_Float *tag = NULL;
    E_Int tag_size = -1;
    ret = K_NUMPY::getFromNumpyArray(TAG, tag, tag_size, true);
    if (ret != 1) {
        Karray_free_ngon(sarray);
        RAISE("Bad slave points tag.");
        return NULL;
    }
 
    // Extract Mf and Sf, the planar surfaces to intersect
    // TODO(Imad): quasi-planar surfaces

    // Init master/slave meshes
    IMesh S(*sarray.cn, sarray.x, sarray.y, sarray.z, sarray.npts);
    S.make_skin();

    for (E_Int fid : S.skin) {
        const auto &pn = S.F[fid];
        size_t stride = pn.size();
        assert(stride == 3 || stride == 4);

        E_Int keep = 1;

        for (size_t j = 0; j < stride; j++) {
            E_Int point = pn[j];
            if (tag[point] == 0) {
                keep = 0;
                break;
            }
        }

        if (keep) S.patch.insert(fid);
    }

    puts("Adapting intersection zones...");

    ret = meshes_mutual_refinement(M, S);
    if (ret != 0) {
        Karray_free_ngon(sarray);
        RAISE("Bad mesh mutual refinement.");
        return NULL;
    }

    puts("Extracting conformized meshes...");

    S = S.extract_conformized();
    M = M.extract_conformized();
    
    // Export
    puts("Exporting to CGNS format...");

    PyObject *Sout = S.export_karray();

    PyObject *Out = PyList_New(0);

    PyList_Append(Out, Sout);

    Py_DECREF(Sout);

    Karray_free_ngon(sarray);

    // Extract master and slave patches
    puts("Saving slave intersection patch...");
    npy_intp dims[2];
    dims[1] = 1;

    dims[0] = (npy_intp)S.patch.size();
    PyArrayObject *SP = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *sptr = (E_Int *)PyArray_DATA(SP);
    E_Int *ptr = sptr;
    for (E_Int face : S.patch) *ptr++ = face + 1;

    PyList_Append(Out, (PyObject *)SP);

    Py_DECREF(SP);
    Py_DECREF(TAG);

    return Out;
}

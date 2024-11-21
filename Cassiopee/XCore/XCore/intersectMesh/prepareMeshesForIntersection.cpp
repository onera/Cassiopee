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
#include "io.h"
#include "triangle.h"

PyObject *K_XCORE::prepareMeshesForIntersection(PyObject *self, PyObject *args)
{
    /*
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
    IMesh S(*sarray.cn, sarray.X(), sarray.Y(), sarray.Z(), sarray.npts);
    S.make_skin();
    S.patch.clear();

    assert(tag_size == S.np);

    for (E_Int fid : S.skin) {
        const auto &pn = S.F[fid];
        size_t stride = pn.size();

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

    printf("S.patch size before triangulation: %lu\n", S.patch.size());

    S.faces_to_tri.clear();
    for (auto fid : S.patch) S.faces_to_tri.insert(fid);
    S.triangulate_face_set(false);

    printf("S.patch size after triangulation: %lu\n", S.patch.size());

    //M.project_patch(S);

    puts("Adapting intersection zones...");

    printf("M.patch size before refinement: %lu\n", M.patch.size());

    ret = meshes_mutual_refinement(M, S);
    if (ret != 0) {
        Karray_free_ngon(sarray);
        RAISE("Bad mesh mutual refinement.");
        return NULL;
    }

    printf("M.patch size after refinement: %lu\n", M.patch.size());
    printf("S.patch size after refinement: %lu\n", S.patch.size());

    puts("Extracting conformized meshes...");

    S = S.extract_conformized();
    M = M.extract_conformized();

    // Export
    puts("Exporting to CGNS format...");

    PyObject *Mout = M.export_karray();
    PyObject *Sout = S.export_karray();

    Karray_free_ngon(sarray);

    // Extract master and slave patches
    puts("Saving intersection patches...");
    npy_intp dims[2];
    dims[1] = 1;

    dims[0] = (npy_intp)M.patch.size();
    PyArrayObject *MP = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *mptr = (E_Int *)PyArray_DATA(MP);
    E_Int *ptr = mptr;
    for (E_Int fid : M.patch) *ptr++ = fid + 1;

    dims[0] = (npy_intp)S.patch.size();
    PyArrayObject *SP = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *sptr = (E_Int *)PyArray_DATA(SP);
    ptr = sptr;
    for (E_Int fid : S.patch) *ptr++ = fid + 1;

    PyObject *Out = PyList_New(0);    

    PyList_Append(Out, Mout);
    PyList_Append(Out, (PyObject *)MP);
    PyList_Append(Out, Sout);
    PyList_Append(Out, (PyObject *)SP);

    Py_DECREF(Mout);
    Py_DECREF(Sout);
    Py_DECREF(MP);
    Py_DECREF(SP);

    Py_DECREF(TAG);

    return Out;
    */
    return Py_None;
}

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
#include "karray.h"
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
    PyObject *MASTER, *SLAVE, *PATCH, *TAG;
  
    if (!PYPARSETUPLE_(args, OOOO_, &MASTER, &SLAVE, &PATCH, &TAG)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    Karray sarray;

    E_Int ret;

    ret = Karray_parse_ngon(MASTER, marray);

    if (ret != 0) return NULL;

    ret = Karray_parse_ngon(SLAVE, sarray);

    if (ret != 0) {
        Karray_free_ngon(marray);
        return NULL;
    }

    // Init and orient master/slave meshes
    IMesh M(*marray.cn, marray.X, marray.Y, marray.Z, marray.npts);
    IMesh S(*sarray.cn, sarray.X, sarray.Y, sarray.Z, sarray.npts);

    // Check intersection patch
    E_Int *mpatch = NULL;
    E_Int mpatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(PATCH, mpatch, mpatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad master patch.");
        return NULL;
    }

    printf("Master patch: " SF_D_ " faces\n", mpatch_size);

    for (E_Int i = 0; i < mpatch_size; i++) M.patch.insert(mpatch[i]-1);

    // Check slave point tags
    E_Float *tag = NULL;
    E_Int tag_size = -1;
    ret = K_NUMPY::getFromNumpyArray(TAG, tag, tag_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad slave points tag.");
        return NULL;
    }
 
    // Extract Mf and Sf, the planar surfaces to intersect
    // TODO(Imad): quasi-planar surfaces
    for (E_Int i = 0; i < S.nf; i++) {
        const auto &pn = S.F[i];
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

        if (keep) S.patch.insert(i);
    }

    ret = meshes_mutual_refinement(M, S);
    if (ret != 0) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        return NULL;
    }

    M = M.extract_conformized();
    S = S.extract_conformized();
    
    // Export
    PyObject *Mout = M.export_karray();
    PyObject *Sout = S.export_karray();

    PyObject *Out = PyList_New(0);
    PyList_Append(Out, Mout);
    PyList_Append(Out, Sout);
    Py_DECREF(Mout);
    Py_DECREF(Sout);
    Karray_free_ngon(marray);
    Karray_free_ngon(sarray);


    // Extract master and slave patches
    npy_intp dims[2];
    dims[1] = 1;
    
    dims[0] = (npy_intp)M.patch.size();
    PyArrayObject *MP = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *mptr = (E_Int *)PyArray_DATA(MP);
    E_Int *ptr = mptr;
    for (E_Int face : M.patch) *ptr++ = face;

    dims[0] = (npy_intp)S.patch.size();
    PyArrayObject *SP = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *sptr = (E_Int *)PyArray_DATA(SP);
    ptr = sptr;
    for (E_Int face : S.patch) *ptr++ = face;

    PyList_Append(Out, (PyObject *)MP);
    PyList_Append(Out, (PyObject *)SP);
    Py_DECREF(MP);
    Py_DECREF(SP);


    return Out;
}

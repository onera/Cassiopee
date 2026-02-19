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
#include "icapsule.h"

PyObject *K_XCORE::icapsule_adapt2(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;

    if (!PYPARSETUPLE_(args, O_, &ICAPSULE)) {
        RAISE("Bad input");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");

    auto &M = icap->M;
    auto &Ss = icap->Ss;

    M.Mf = Smesh(M, M.skin, false);
    auto &Mf = M.Mf;

    for (size_t i = 0; i < Ss.size(); i++) {
        printf("Adapting for slave %zu\n", i);

        Mf.make_fcenters();
        Mf.make_fnormals();
        Mf.make_pnormals();
        Mf.make_point_faces();
        Mf.make_BVH();

        auto &S = Ss[i];

        Smesh Sf = Smesh::Smesh_from_point_tags(S, S.ptag.data(), true);

        Sf.make_fcenters();
        Sf.make_fnormals();
        Sf.make_pnormals();
        Sf.compute_min_distance_between_points();
        printf("Min dist: %f\n", Sf.min_pdist);

        //Sf.write_ngon("Sf_before_inter.im");

        // Locate Sf points on Mf faces
        auto plocs = Mf.locate2(Sf);
        std::vector<E_Int> spids(Sf.np);
        for (int i = 0; i < Sf.np; i++) spids[i] = i;
        Sf.replace_by_projections(spids, plocs);

        // Extract the initial Mf faces that cover Sf
        auto bfaces = Mf.extract_covering_faces(Sf, plocs);

        // Refinement loop
        ICapsule::refine(Mf, bfaces, Sf);

        // Reconstruct S
        Sf.reconstruct(S);

        // Reconstruct S skin
        S.make_skin();

        // Tag Sf faces
        Sf.tag_faces(S);
    }

    Mf.reconstruct(M);

    Py_INCREF(Py_None);
    return Py_None;

    /*
    PyObject *out = PyList_New(0);
    PyList_Append(out, M.export_karray());
    
    PyObject *slist = PyList_New(0);
    for (const auto &S : Ss) {
        PyList_Append(slist, S.export_karray());
    }
    PyList_Append(out, slist);
    Py_DECREF(slist);

    PyObject *tlist = PyList_New(0);
    for (const auto &S : Ss) {
        npy_intp dims[2];
        dims[0] = (npy_intp)S.ftag.size();
        dims[1] = 1;
        PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
        E_Int *ptr = (E_Int *)PyArray_DATA(arr);
        for (size_t i = 0; i < S.ftag.size(); i++) ptr[i] = S.ftag[i]+1;
        PyList_Append(tlist, (PyObject *)arr);
        Py_DECREF(arr);
    }
    PyList_Append(out, tlist);
    Py_DECREF(tlist);

    return out;
    */
}
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
#include <queue>
#include "icapsule.h"
#include "common/Karray.h"
#include "point.h"
#include "io.h"
#include "primitives.h"
#include "dcel.h"

PyObject *K_XCORE::icapsule_init2(PyObject *self, PyObject *args)
{
    ICapsule *ic = new ICapsule;
    PyObject *hook = PyCapsule_New((void *)ic, "ICAPSULE", NULL);
    return hook;
}

PyObject *K_XCORE::icapsule_set_master(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE, *MASTER, *CTAG;
    if (!PYPARSETUPLE_(args, OOO_, &ICAPSULE, &MASTER, &CTAG)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    Karray array;
    if (Karray_parse_ngon(MASTER, array) != 0) {
        RAISE("Failed to parse master NGon.");
        return NULL;
    }

    E_Float *ctag = NULL;
    E_Int size;
    E_Int ret = K_NUMPY::getFromNumpyArray(CTAG, ctag, size);
    if (ret != 1 || size != array.ncells()) 
    {
        RAISE("Bad cell tag array.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");
    icap->M = IMesh(array);

    auto &M = icap->M;
    M.set_tolerances(icap->NEAR_VERTEX_TOL, icap->NEAR_EDGE_TOL);
    M.make_skin();
    M.orient_skin(OUT);
    M.ctag.resize(M.nc);
    memcpy(M.ctag.data(), ctag, M.nc * sizeof(E_Float));

    Py_INCREF(Py_None);
    return Py_None;
}

PyObject *K_XCORE::icapsule_set_slaves(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE, *SLAVES, *PTAGS, *CTAGS;
    if (!PYPARSETUPLE_(args, OOOO_, &ICAPSULE, &SLAVES, &PTAGS, &CTAGS)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");

    if (!PyList_Check(SLAVES)) {
        RAISE("Slaves should be a list of NGons.");
        return NULL;
    }

    if (!PyList_Check(PTAGS)) {
        RAISE("Point tags should be a list of arrays.");
        return NULL;
    }

    if (!PyList_Check(CTAGS)) {
        RAISE("Cell tags should be a list of arrays.");
        return NULL;
    }

    E_Int nslaves = PyList_Size(SLAVES);

    // Meshes

    std::vector<Karray> sarrays(nslaves); 
    E_Int ok = 0;

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *SLAVE = PyList_GetItem(SLAVES, i);
        E_Int ret = Karray_parse_ngon(SLAVE, sarrays[i]);
        Py_DECREF(SLAVE);
        if (ret != 0) {
            RAISE("Slaves should be NGons.");
            break;
        }
        ok++;
    }

    if (ok != nslaves) {
        return NULL;
    }

    // Point tags

    std::vector<E_Float *> ptags(nslaves, NULL);
    
    if (PyList_Size(PTAGS) != nslaves) {
        RAISE("Ptags should be the same size as slaves.");
        return NULL;
    }

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *PTAG = PyList_GetItem(PTAGS, i);
        E_Int size = -1;
        E_Int ret = K_NUMPY::getFromNumpyArray(PTAG, ptags[i], size);
        //Py_DECREF(PTAG);
        if (ret != 1 || size != sarrays[i].npoints()) {
            RAISE("Ptags[i] should have size sarrays[i].npoints.");
            return NULL;
        }
    }

    // Cell tags

    std::vector<E_Float *> ctags(nslaves, NULL);
    
    if (PyList_Size(CTAGS) != nslaves) {
        RAISE("Cell tags list should be the same size as slaves.");
        return NULL;
    }

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *CTAG = PyList_GetItem(CTAGS, i);
        E_Int size = -1;
        E_Int ret = K_NUMPY::getFromNumpyArray(CTAG, ctags[i], size);
        //Py_DECREF(PTAG);
        if (ret != 1 || size != sarrays[i].ncells()) {
            RAISE("Ctags[i] should have size sarrays[i].npoints.");
            return NULL;
        }
    }

    icap->Ss.reserve(sarrays.size());
    for (size_t i = 0; i < sarrays.size(); i++) {
        icap->Ss.push_back(IMesh(sarrays[i]));
        auto &S = icap->Ss[i];
        S.set_tolerances(icap->NEAR_VERTEX_TOL, icap->NEAR_EDGE_TOL);
        S.make_skin();
        S.orient_skin(IN);
        S.ptag.resize(S.np);
        memcpy(S.ptag.data(), ptags[i], S.np*sizeof(E_Float));
        // Triangulate faces with tagged points
        std::vector<E_Int> fids;
        for (E_Int fid = 0; fid < S.nf; fid++) {
            const auto &pn = S.F[fid];
            bool triangulate = true;
            for (auto p : pn) {
                if (S.ptag[p] != 1.0) {
                    triangulate = false;
                    break;
                }
            }
            if (triangulate) fids.push_back(fid);
        }
        S.triangulate(fids);
        // Cell tags
        S.ctag.resize(S.nc);
        memcpy(S.ctag.data(), ctags[i], S.nc*sizeof(E_Float));
    }

    Py_INCREF(Py_None);
    return Py_None;
}

ICapsule::ICapsule(const Karray &marray, const std::vector<Karray> &sarrays,
    const std::vector<E_Float *> &ptags)
{
    E_Float NEAR_VERTEX_TOL = 1e-3;
    E_Float NEAR_EDGE_TOL = 1e-3;

    M = IMesh(marray);
    M.set_tolerances(NEAR_VERTEX_TOL, NEAR_EDGE_TOL);
    M.make_skin();
    M.orient_skin(OUT);

    Ss.reserve(sarrays.size());
    for (size_t i = 0; i < sarrays.size(); i++) {
        Ss.push_back(IMesh(sarrays[i]));
        Ss[i].set_tolerances(NEAR_VERTEX_TOL, NEAR_EDGE_TOL);
        Ss[i].make_skin();
        Ss[i].orient_skin(IN);
        Ss[i].ptag.resize(Ss[i].np);
        memcpy(Ss[i].ptag.data(), ptags[i], Ss[i].np*sizeof(E_Float));
        // Triangulate faces with tagged points
        std::vector<E_Int> fids;
        for (E_Int fid = 0; fid < Ss[i].nf; fid++) {
            const auto &pn = Ss[i].F[fid];
            bool triangulate = true;
            for (auto p : pn) {
                if (Ss[i].ptag[p] != 1.0) {
                    triangulate = false;
                    break;
                }
            }
            if (triangulate) fids.push_back(fid);
        }
        Ss[i].triangulate(fids);
    }
}

PyObject *K_XCORE::icapsule_init(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVES, *PTAGS;
    if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVES, &PTAGS)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    if (Karray_parse_ngon(MASTER, marray) != 0) {
        RAISE("Master should be an NGon.");
        return NULL;
    }

    E_Int nslaves = PyList_Size(SLAVES);
    std::vector<Karray> sarrays(nslaves); 
    E_Int ok = 0;

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *SLAVE = PyList_GetItem(SLAVES, i);
        E_Int ret = Karray_parse_ngon(SLAVE, sarrays[i]);
        Py_DECREF(SLAVE);
        if (ret != 0) {
            RAISE("Slaves should be NGons.");
            break;
        }
        ok++;
    }

    if (ok != nslaves) {
        Karray_free_ngon(marray);
        for (E_Int i = 0; i < ok; i++)
            Karray_free_ngon(sarrays[i]);
        return NULL;
    }

    std::vector<E_Float *> ptags(nslaves, NULL);
    
    if (PyList_Size(PTAGS) != nslaves) {
        RAISE("Ptags should be the same size as slaves.");
        Karray_free_ngon(marray);
        for (E_Int i = 0; i < nslaves; i++)
            Karray_free_ngon(sarrays[i]);
        return NULL;
    }

    for (E_Int i = 0; i < nslaves; i++) {
        PyObject *PTAG = PyList_GetItem(PTAGS, i);
        E_Int size = -1;
        E_Int ret = K_NUMPY::getFromNumpyArray(PTAG, ptags[i], size);
        //Py_DECREF(PTAG);
        if (ret != 1 || size != sarrays[i].npoints()) {
            RAISE("Ptag[i] should have size sarrays[i].npoints.");
            Karray_free_ngon(marray);
            for (E_Int j = 0; j < nslaves; j++)
                Karray_free_ngon(sarrays[j]);
            return NULL;
        }
    }

    ICapsule *icap = new ICapsule(marray, sarrays, ptags);

    PyObject *hook = PyCapsule_New((void *)icap, "ICAPSULE", NULL);

    return hook;
}

PyObject *K_XCORE::icapsule_adapt(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;
    if (!PYPARSETUPLE_(args, O_, &ICAPSULE)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE")) {
        RAISE("Bad capsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");

    auto &M = icap->M;
    auto &Ss = icap->Ss;

    Smesh Mf(M, M.skin, false);
    puts("Making surface mesh data");
    
    Mf.write_ngon("Mf_before_inter.im");

    //for (E_Int i = 25; i < 26; i++) {
    for (size_t i = 0; i < Ss.size(); i++) {
        printf("S%zu\n", i);

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
}

PyObject *K_XCORE::icapsule_intersect(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVES, *STAGS;
    if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVES, &STAGS)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    int ret = Karray_parse_ngon(MASTER, marray);
    if (ret != 0) {
        RAISE("Bad master mesh.");
        return NULL;
    }
    IMesh M(marray);

    int slave_count = PyList_Size(SLAVES);
    std::vector<IMesh> Ss;
    Ss.reserve(slave_count);

    assert(PyList_Size(SLAVES) == PyList_Size(STAGS));

    for (int i = 0; i < slave_count; i++) {
        Karray sarray;
        ret = Karray_parse_ngon(PyList_GetItem(SLAVES, i), sarray);
        assert(ret == 0);
        IMesh S(sarray);
        S.make_skin();
        S.orient_skin(IN);

        PyObject *STAG = PyList_GetItem(STAGS, i);
        E_Int *tag = NULL;
        E_Int tag_size = -1;
        ret = K_NUMPY::getFromNumpyArray(STAG, tag, tag_size);
        assert(ret == 1);
        S.ftag.reserve(tag_size);
        for (int j = 0; j < tag_size; j++) {
            E_Int fid = tag[j]-1;
            assert(fid < S.nf);
            S.ftag.push_back(fid);
        }

        Ss.push_back(S);
    }

    M.make_skin();
    M.orient_skin(OUT);
    Smesh Mf(M, M.skin, false);

    for (size_t i = 0; i < Ss.size(); i++) {
        
        printf("Intersecting slave %zu\n", i);

        //Mf.make_bbox();
        //Mf.hash_faces();
        Mf.make_fcenters();
        Mf.make_fnormals();
        Mf.make_pnormals();
        Mf.make_point_faces();
        Mf.make_BVH();

        //Mf.write_ngon("Mf_before_intersect.im");

        auto &S = Ss[i];

        Smesh Sf = Smesh::Smesh_from_tagged_faces(S, true);
        Sf.make_fcenters();
        Sf.make_fnormals();
        Sf.make_pnormals();

        //{
        //    char fname[32] = {0};
        //    sprintf(fname, "Sf_before_intersect_%d.im", i);
        //    Sf.write_ngon(fname);
        //}
        
        
        auto plocs = Mf.locate2(Sf);
        
        /*
        auto bfaces = Mf.extract_covering_faces(Sf, plocs);
        Smesh Bf(Mf, bfaces, false);
        Bf.write_ngon("Bf.im");
        Bf.make_fcenters();
        Bf.make_fnormals();
        Bf.make_pnormals();
        Bf.make_point_faces();
        Bf.make_BVH();
        plocs = Bf.locate2(Sf);
        */

        Dcel D = Dcel::intersect(Mf, Sf, plocs);
        //Dcel D = Dcel::intersect(Bf, Sf, plocs);

        //E_Int nf_before_intersect = Sf.nf;

        D.reconstruct(Mf, Dcel::RED);
        //D.reconstruct(Bf, Dcel::RED);

        
        //Mf.write_ngon("Mf_after_intersect.im");
        
        D.reconstruct(Sf, Dcel::BLACK);

        //{
        //    char fname[32] = {0};
        //    sprintf(fname, "Sf_after_intersect_%d.im", i);
        //    Sf.write_ngon(fname);
        //}

        //{
        //    char fname[32] = {0};
        //    sprintf(fname, "intersected_%d.im", i);
        //    D.write_inner_cycles(fname);
        //}

        Sf.reconstruct(S);

        // Tag Sf faces
        Sf.tag_faces(S);

        puts("Intersection OK\n");
        fflush(stdout);
    }
    //exit(0);
    
    Mf.reconstruct(M);

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
}

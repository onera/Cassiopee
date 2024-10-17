#include <queue>

#include "icapsule.h"
#include "common/Karray.h"
#include "point.h"
#include "io.h"
#include "primitives.h"
#include "dcel.h"

ICapsule::ICapsule(const Karray &marray, const std::vector<Karray> &sarrays,
    const std::vector<E_Float *> &ptags)
{
    M = IMesh(marray);
    M.set_tolerances(NEAR_VERTEX_TOL, NEAR_EDGE_TOL);
    M.make_skin();
    M.orient_skin(OUT);
    M.triangulate_skin();
    M.make_bbox();
    M.hash_skin();
    //M.make_skin_graph();
    Smesh Mf(M, M.skin, false);
    Mf.make_bbox();
    Mf.hash_faces();
    Mf.make_fcenters();
    Mf.make_fnormals();
    Mf.make_pnormals();

    Ss.reserve(sarrays.size());
    spatches.reserve(sarrays.size());

    std::vector<std::vector<PointLoc>> plocs_list;
    plocs_list.reserve(sarrays.size());
    std::vector<std::set<E_Int>> bfaces_list;
    bfaces_list.reserve(sarrays.size());

    for (size_t i = 0; i < sarrays.size(); i++) {

        printf("S%lu\n", i);

        // Create IMesh S
        IMesh S(sarrays[i]);
        S.set_tolerances(NEAR_VERTEX_TOL, NEAR_EDGE_TOL);
        S.make_skin();
        S.orient_skin(OUT);
        S.triangulate_skin();
        Ss.push_back(S);

        // Create SMesh Sf
        Smesh Sf = Smesh::Smesh_from_point_tags(S, ptags[i]);
        Sf.make_fcenters();
        Sf.make_fnormals();
        Sf.make_pnormals();
        Sf.compute_min_distance_between_points();

        // Locate Sf points on Mf faces
        auto plocs = Mf.locate(Sf);
        //Mf.correct_near_points_and_edges(Sf, plocs);

        // Correct AABB and Sf faces hash
        Sf.make_bbox();
        Sf.hash_faces();

        // Extract the initial Mf faces that bound Sf
        auto bfaces = Mf.extract_bounding_faces(Sf, plocs);
        //Mf.write_ngon("bounding_before", bfaces);

        //Sf.write_ngon("Sf");

        // Refinement loop
        refine(Mf, bfaces, Sf, plocs);

        bfaces = Mf.extract_bounding_faces(Sf, plocs);
        //Mf.write_ngon("bounding_after.im", bfaces);

        // Make a patch out of Mf bounding faces
        Smesh Bf(Mf, bfaces, false);
        Bf.make_fcenters();
        Bf.make_fnormals();
        Bf.make_pnormals();
        
        // Make two DCELs out of Bf and Sf
        Dcel Db(Bf, Dcel::RED);
        Db.write_inner_cycles("Db.im");

        Sf.write_ngon("refinement_Sf.im");

        Dcel Ds(Sf, Dcel::BLACK);
        Ds.write_inner_cycles("Ds_inner.im");
        Ds.write_outer_cycles("Ds_outer.im");
        Ds.write_hole_cycles("Ds_hole.im");

        // Add Sf
        spatches.push_back(Sf);

        // Add plocs
        plocs_list.push_back(plocs);

        puts("");
    }

    Mf.write_ngon("refined_Mf.im");


}


PyObject *K_XCORE::icapsule_extract_master(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;
    if (!PYPARSETUPLE_(args, O_, &ICAPSULE)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICapsule")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICapsule");

    auto Mout = icap->M.export_karray();

    return Mout;
}

PyObject *K_XCORE::icapsule_extract_slave(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;
    E_Int INDEX;
    if (!PYPARSETUPLE_(args, O_ I_, &ICAPSULE, &INDEX)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICapsule")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICapsule");

    if (INDEX >= (E_Int)icap->Ss.size()) {
        RAISE("Bad slave index.");
        return NULL;
    }

    auto Sout = icap->Ss[INDEX].export_karray();

    return Sout;
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
        E_Int ret = K_NUMPY::getFromNumpyArray(PTAG, ptags[i], size, true);
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

    PyObject *out = PyCapsule_New((void *)icap, "ICapsule", NULL);

    Karray_free_ngon(marray);
    for (E_Int i = 0; i < nslaves; i++)
        Karray_free_ngon(sarrays[i]);

    return out;
}

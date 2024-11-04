#include <queue>

#include "icapsule.h"
#include "common/Karray.h"
#include "point.h"
#include "io.h"
#include "primitives.h"
#include "dcel.h"

ICapsule ICapsule::do_it(const Karray &marray,
    const std::vector<Karray> &sarrays, const std::vector<E_Float *> &ptags)
{
    ICapsule icap;
    auto &M = icap.M;
    auto &Ss = icap.Ss;

    E_Float NEAR_VERTEX_TOL = 1e-3;
    E_Float NEAR_EDGE_TOL = 1e-3;

    puts("Creating mesh from karray");
    M = IMesh(marray);
    M.set_tolerances(NEAR_VERTEX_TOL, NEAR_EDGE_TOL);
    puts("Making skin");
    M.make_skin();
    puts("Orienting skin");
    M.orient_skin(OUT);
    puts("Triangulating skin");
    M.triangulate_skin();
    puts("Making edge data");
    M.make_edges();
    puts("Making bbox");
    M.make_bbox();
    puts("Hashing skin");
    M.hash_skin();
    puts("Creating surface mesh");
    Smesh Mf(M, M.skin, false);
    puts("Making surface mesh data");
    Mf.make_bbox();
    Mf.hash_faces();
    Mf.make_fcenters();
    Mf.make_fnormals();
    Mf.make_pnormals();

    Mf.write_ngon("Mf_raw.im");

    Ss.reserve(sarrays.size());

    // Adapt

    for (size_t i = 0; i < sarrays.size(); i++) {

        printf("S%lu\n", i);

        // Create IMesh S
        IMesh S(sarrays[i]);
        S.set_tolerances(NEAR_VERTEX_TOL, NEAR_EDGE_TOL);
        S.make_skin();
        S.orient_skin(IN);
        S.triangulate_skin();
        S.make_edges();

        // Create SMesh Sf
        Smesh Sf = Smesh::Smesh_from_point_tags(S, ptags[i], true);

        for (E_Int fid = 0; fid < Sf.nf; fid++) {
            E_Int gfid = Sf.l2gf.at(fid);
            assert(S.neigh[gfid] == -1);
        }

        Sf.make_fcenters();
        Sf.make_fnormals();
        Sf.make_pnormals();
        Sf.make_bbox();
        Sf.hash_faces();
        Sf.compute_min_distance_between_points();
        printf("Min dist: %f\n", Sf.min_pdist);

        Sf.write_ngon("Sf_raw.im");

        // Locate Sf points on Mf faces
        auto plocs = Mf.locate(Sf);
        std::vector<E_Int> spids(Sf.np);
        for (int i = 0; i < Sf.np; i++) spids[i] = i;
        Sf.replace_by_projections(spids, plocs);

        // Extract the initial Mf faces that cover Sf
        auto bfaces = Mf.extract_covering_faces(Sf, plocs);

        // Refinement loop
        refine(Mf, bfaces, Sf);

        // Reconstruct S
        Sf.reconstruct(S);

        // TODO(Imad): Tag Sf faces
        Sf.tag_faces(S);

        // Append
        Ss.push_back(S);

        puts("");
    }

    Mf.reconstruct(M);

    return icap;

    /*
    // Intersect

    for (size_t i = 0; i < sarrays.size(); i++) {

        // Make Sf
        Smesh Sf = Smesh::Smesh_from_tagged_faces(Ss[i], true);

        // Make Bf
        M.make_skin();
        Smesh Mf(M, M.skin, false);
        auto plocs = Mf.locate(Sf);

        auto bfaces = Mf.extract_covering_faces(Sf, plocs);
        Smesh Bf = Smesh::make_sub(Mf, bfaces, false);
        Bf.make_fcenters();
        Bf.make_fnormals();
        Bf.make_pnormals();
        Bf.make_bbox();
        Bf.hash_faces();

        plocs = Bf.locate(Sf);

        // Intersect
        Dcel D(Bf, Sf, plocs);

        // Reconstruct Bf and Sf
        D.reconstruct_smesh(Bf, Dcel::RED, false);
        D.reconstruct_smesh(Sf, Dcel::BLACK, true);

        {
            Smesh new_Bf = D.reconstruct(Bf, Dcel::RED, false);
            new_Bf.write_ngon("new_Bf.im");
            puts("reconstruted new Bf");
            Smesh new_Sf = D.reconstruct(Sf, Dcel::BLACK, true);
            new_Sf.write_ngon("new_Sf.im");
            puts("reconstruted new Sf");
        }

        // Reconstruct M and S
        Bf.reconstruct(M);
        Sf.reconstruct(Ss[i]);

        puts("");
    }

    Mf.write_ngon("refined_Mf.im");
    */
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

    ICapsule *icap = new ICapsule();
    *icap = ICapsule::do_it(marray, sarrays, ptags);

    Karray_free_ngon(marray);
    for (E_Int i = 0; i < nslaves; i++)
        Karray_free_ngon(sarrays[i]);
    
    PyObject *hook = PyCapsule_New((void *)icap, "ICapsule", NULL);

    return hook;
}

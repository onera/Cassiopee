#include "icapsule.h"
#include "common/Karray.h"
#include "point.h"
#include "io.h"
#include "primitives.h"

Smesh IMesh::make_smesh(const E_Float *ptag)
{
    patch.clear();

    for (E_Int fid : skin) {
        const auto &pn = F[fid];
        bool is_patch = true;
        for (E_Int pid : pn) {
            if (ptag[pid] != 1) {
                is_patch = false;
                break;
            }
        }
        if (is_patch)
            patch.insert(fid);
    }
    
    return Smesh(*this);
}

Smesh IMesh::make_smesh_from_skin(bool is_planar)
{
    patch.clear();

    for (E_Int fid : skin) {
        patch.insert(fid);
    }
    
    return Smesh(*this, is_planar);
}

std::vector<Point> epoints;
std::vector<Point> cpoints;

void Smesh::correct_near_points_and_edges(Smesh &Sf,
    std::vector<PointLoc> &plocs)
{
    E_Int on_vertex = 0, on_edge = 0;
    for (size_t i = 0; i < plocs.size(); i++) {
        auto &ploc = plocs[i];

        E_Int fid = ploc.fid;
        assert(fid < nf);
        const auto &pn = F[fid];

        if (ploc.v_idx != -1) {
            on_vertex++;
            E_Int p = pn[ploc.v_idx];
            E_Float dx = X[p]-Sf.X[i];
            E_Float dy = Y[p]-Sf.Y[i];
            E_Float dz = Z[p]-Sf.Z[i];
            E_Float dist = dx*dx + dy*dy + dz*dz;
            if (dist >= Sf.min_pdist_squared) {
                fprintf(stderr, "Tight near-vertex situation!\n");
                point_write("mpoint", X[p], Y[p], Z[p]);
                point_write("spoint", Sf.X[i], Sf.Y[i], Sf.Z[i]);
                assert(0);
            } else {
                Sf.X[i] = X[p];
                Sf.Y[i] = Y[p];
                Sf.Z[i] = Z[p];
            }
        } else if (ploc.e_idx != -1) {
            on_edge++;
            E_Int zero_crd = (ploc.e_idx+2)%3;
            E_Float U = ploc.bcrd[zero_crd];
            assert(Sign(U, NEAR_EDGE_TOL) == 0);
            E_Int i1 = (zero_crd+1)%3;
            E_Int i2 = (zero_crd+2)%3;
            E_Float V = ploc.bcrd[i1];
            E_Float W = ploc.bcrd[i2];
            E_Int a = pn[i1];
            E_Int b = pn[i2];
            V += U;
            assert(Sign(V+W-1) == 0);
            epoints.push_back(Point(Sf.X[i], Sf.Y[i], Sf.Z[i]));
            Sf.X[i] = V*X[a] + W*X[b];
            Sf.Y[i] = V*Y[a] + W*Y[b];
            Sf.Z[i] = V*Z[a] + W*Z[b];
            cpoints.push_back(Point(Sf.X[i], Sf.Y[i], Sf.Z[i]));
        }
    }
    printf("on vertex: %d - on edge: %d\n", on_vertex, on_edge);
}

std::vector<Point> chain_points;

Smesh Smesh::extract_bounding_smesh(const Smesh &Sf,
    const std::vector<PointLoc> &plocs)
{
    // Get boundary edges from spatch
    std::set<E_Int> bedges;
    for (size_t i = 0; i < Sf.E2F.size(); i++) {
        const auto &pf = Sf.E2F[i];
        assert(pf[0] != -1);
        if (pf[1] == -1) bedges.insert(i);
    }
    size_t nbedges = bedges.size();

    // Make the boundary point chain
    std::vector<E_Int> pchain;

    const auto &E = Sf.E;

    E_Int first_edge = *bedges.begin();

    pchain.push_back(E[first_edge].p);
    pchain.push_back(E[first_edge].q);

    bedges.erase(first_edge);

    E_Int current_point = pchain[1];

    while (pchain.size() < nbedges) {
        E_Int to_delete = -1;
        for (auto e : bedges) {
            if (E[e].p == current_point) {
                pchain.push_back(E[e].q);
                current_point = pchain.back();
                to_delete = e;
                break;
            } else if (E[e].q == current_point) {
                pchain.push_back(E[e].p);
                current_point = pchain.back();
                to_delete = e;
                break;
            }
        }
        assert(to_delete != -1);
        bedges.erase(to_delete);
    }

    assert(pchain.size() == nbedges);

    for (auto p : pchain) {
        chain_points.push_back(Point(Sf.X[p], Sf.Y[p], Sf.Z[p]));
    }

    /*
    std::set<E_Int> bfids;
    
    for (size_t i = 0; i < pchain.size(); i++) {
        E_Int p = pchain[i];
        E_Int q = pchain[(i+1)%pchain.size()];

        E_Float px = Mf.X[p], py = Mf.Y[p], pz = Mf.Z[p];
        E_Float qx = Mf.X[q], qy = Mf.Y[q], qz = Mf.Z[q];

        E_Float D[3] = {qx-px, qy-py, qz-pz};
        E_Float NORM = K_MATH::norm(D, 3);
        D[0] /= NORM, D[1] /= NORM, D[2] /= NORM;

        std::vector<E_Int> orig_faces;
        std::vector<E_Int> tail_faces;

        E_Int last_vertex = -1, last_edge = -1, dummy;

        Mf.get_shared_faces(plocs[p], orig_faces, last_vertex, last_edge); 
        Mf.get_shared_faces(plocs[q], tail_faces, dummy, dummy); 

        E_Int starting_face = Mf.deduce_face(orig_faces, px, py, pz,
            D, last_vertex, last_edge);
        assert(starting_face != -1);

        bool found_tail = false;
        E_Int current_fid = starting_face;
        E_Float current_pos[3] = {px, py, pz};

        E_Int walk = 0;
        E_Int max_walks = 20;

        while (!found_tail && walk <= max_walks) {
            found_tail = true;
            walk++;
        }

        assert(found_tail);
        assert(walk <= max_walks);


    }
    */

    return Smesh();
}

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
    M.make_skin_graph();
    Smesh Mf = M.make_smesh_from_skin(false);
    Mf.make_bbox();
    Mf.hash_faces();
    Mf.make_fnormals();

    Ss.reserve(sarrays.size());
    spatches.reserve(sarrays.size());

    std::vector<std::vector<PointLoc>> plocs;
    plocs.reserve(sarrays.size());
    mpatches.reserve(sarrays.size());

    for (size_t i = 0; i < sarrays.size(); i++) {
        Ss.push_back(IMesh(sarrays[i]));
        Ss[i].set_tolerances(NEAR_VERTEX_TOL, NEAR_EDGE_TOL);
        Ss[i].make_skin();
        Ss[i].orient_skin(IN);
        Ss[i].triangulate_skin();
        spatches.push_back(Ss[i].make_smesh(ptags[i]));
        spatches[i].make_bbox();
        spatches[i].hash_faces();
        spatches[i].compute_min_distance_between_points();
        plocs.push_back(Mf.locate(spatches[i]));
        Mf.correct_near_points_and_edges(spatches[i], plocs[i]);
        mpatches.push_back(Mf.extract_bounding_smesh(spatches[i], plocs[i]));
    }

    point_write("epoints", epoints);
    point_write("cpoints", cpoints);

    point_write("chain", chain_points);
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
        Py_DECREF(PTAG);
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

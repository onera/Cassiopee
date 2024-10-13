#include "icapsule.h"
#include "common/Karray.h"
#include "point.h"
#include "io.h"
#include "primitives.h"
#include <queue>

E_Int ray_point_orient(const E_Float o[3], const E_Float d[3],
    const E_Float fN[3], E_Float px, E_Float py, E_Float pz)
{
    E_Float w[3] = {px-o[0], py-o[1], pz-o[2]};
    E_Float c[3];
    K_MATH::cross(d, w, c);
    E_Float dp = K_MATH::dot(c, fN, 3);
    // TODO(Imad): needs FEA + float128
    E_Int cmp = Sign(dp);
    if (cmp > 0) return 1;
    if (cmp < 0) return -1;
    return 0;
}

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
            Sf.X[i] = V*X[a] + W*X[b];
            Sf.Y[i] = V*Y[a] + W*Y[b];
            Sf.Z[i] = V*Z[a] + W*Z[b];
        }
    }
    //printf("on vertex: %d - on edge: %d\n", on_vertex, on_edge);
}

std::set<E_Int> ewalls;
std::set<E_Int> fwalls;
std::vector<Point> pchains;

Smesh Smesh::extract_bounding_smesh(const Smesh &Sf,
    const std::vector<PointLoc> &plocs) const
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

    E_Int first_edge = *bedges.begin();

    pchain.push_back(Sf.E[first_edge].p);
    pchain.push_back(Sf.E[first_edge].q);

    bedges.erase(first_edge);

    E_Int current_point = pchain[1];

    while (pchain.size() < nbedges) {
        E_Int to_delete = -1;
        for (auto e : bedges) {
            if (Sf.E[e].p == current_point) {
                pchain.push_back(Sf.E[e].q);
                current_point = pchain.back();
                to_delete = e;
                break;
            } else if (Sf.E[e].q == current_point) {
                pchain.push_back(Sf.E[e].p);
                current_point = pchain.back();
                to_delete = e;
                break;
            }
        }
        assert(to_delete != -1);
        bedges.erase(to_delete);
    }

    assert(pchain.size() == nbedges);

    // Sort the pchain counterclockwise
    E_Int a = pchain[0], b = pchain[1], c = pchain[2];
    E_Float ux = Sf.X[b] - Sf.X[a];
    E_Float uy = Sf.Y[b] - Sf.Y[a];
    E_Float uz = Sf.Z[b] - Sf.Z[a];
    E_Float vx = Sf.X[c] - Sf.X[b];
    E_Float vy = Sf.Y[c] - Sf.Y[b];
    E_Float vz = Sf.Z[c] - Sf.Z[b];
    E_Float cp[3] = {uy*vz - uz*vy, uz*vx - ux*vz, ux*vy - uy*vx};
    const E_Float *N_b = &fnormals[3*plocs[b].fid];
    E_Float dp = K_MATH::dot(cp, N_b, 3);
    E_Int cmp = Sign(dp);
    assert(cmp != 0);
    if (cmp < 0)
        std::reverse(pchain.begin(), pchain.end());

    Sf.write_points("pchain", pchain);
    for (E_Int p : pchain) pchains.push_back({Sf.X[p], Sf.Y[p], Sf.Z[p]});

    std::set<E_Int> wfids;
    std::set<E_Int> weids;
 
    for (size_t i = 0; i < pchain.size(); i++) {
        E_Int p = pchain[i];
        E_Int q = pchain[(i+1)%pchain.size()];

        E_Float px = Sf.X[p], py = Sf.Y[p], pz = Sf.Z[p];
        E_Float qx = Sf.X[q], qy = Sf.Y[q], qz = Sf.Z[q];

        E_Float D[3] = {qx-px, qy-py, qz-pz};
        E_Float NORM = K_MATH::norm(D, 3);
        D[0] /= NORM, D[1] /= NORM, D[2] /= NORM;

        std::vector<E_Int> orig_faces;
        std::vector<E_Int> tail_faces;

        E_Int last_vertex = -1, last_edge = -1, dummy;

        get_shared_faces(plocs[p], orig_faces, last_vertex, last_edge); 
        get_shared_faces(plocs[q], tail_faces, dummy, dummy); 

        E_Int starting_face = deduce_face(orig_faces, px, py, pz,
            D, last_vertex, last_edge);
        assert(starting_face != -1);

        bool found_tail = false;
        E_Int cur_fid = starting_face;
        E_Float cur_pos[3] = {px, py, pz};

        E_Int walk = 0;
        E_Int max_walks = 20;

        while (!found_tail && walk <= max_walks) {
            
            wfids.insert(cur_fid);

            E_Float proj[3];
            get_unit_projected_direction(cur_fid, D, proj);

            const auto &pn = F[cur_fid];
            const auto &pe = F2E[cur_fid];
            const E_Float *fN = &fnormals[3*cur_fid];

            // First pass: define the wall data
            for (size_t i = 0; i < pn.size(); i++) {
                E_Int p = pn[i];
                E_Int q = pn[(i+1)%pn.size()];
                E_Int e = pe[i];
                E_Float px = X[p], py = Y[p], pz = Z[p];
                E_Float qx = X[q], qy = Y[q], qz = Z[q];
                if (ray_point_orient(cur_pos, proj, fN, px, py, pz) <= 0 ||
                    ray_point_orient(cur_pos, proj, fN, qx, qy, qz) <= 0) {
                    weids.insert(e);
                }
            }

            for (auto fid : tail_faces) {
                if (fid == cur_fid) {
                    found_tail = true;
                    break;
                }
            }

            if (found_tail) break;

            E_Int next_fid = -1;
            E_Float next_pos[3] = {EFLOATMAX, EFLOATMAX, EFLOATMAX};

            bool hit = false;
            
            for (size_t i = 0; i < pn.size(); i++) {
                E_Int p = pn[i];
                E_Int q = pn[(i+1)%pn.size()];
                E_Int e = pe[i];


                if (p == last_vertex || q == last_vertex || e == last_edge)
                    continue;
                
                E_Float px = X[p], py = Y[p], pz = Z[p];
                E_Float qx = X[q], qy = Y[q], qz = Z[q];
            
                E_Float t, s;
                hit = ray_edge_intersect(
                    cur_pos[0], cur_pos[1], cur_pos[2],
                    proj[0], proj[1], proj[2],
                    px, py, pz, qx, qy, qz,
                    t, s
                );

                if (hit) {
                    if (s > TOL && s < 1 - TOL) {
                        const auto &pe = F2E[cur_fid];
                        E_Int eid = pe[i];
                        last_edge = eid;
                        last_vertex = -1;
                        if (E2F[eid][0] == cur_fid) next_fid = E2F[eid][1];
                        else next_fid = E2F[eid][0];

                        next_pos[0] = cur_pos[0] + t * proj[0];
                        next_pos[1] = cur_pos[1] + t * proj[1];
                        next_pos[2] = cur_pos[2] + t * proj[2];
                    } else {
                        bool hit_p = (s <= TOL);
                        bool hit_q = (s >= 1 - TOL);
                        assert(!(hit_p && hit_q));
                        last_edge = -1;
                        if (hit_p) last_vertex = p;
                        else last_vertex = q;
                        next_pos[0] = X[last_vertex];
                        next_pos[1] = Y[last_vertex];
                        next_pos[2] = Z[last_vertex];
                        const auto &pf = P2F[last_vertex];
                        next_fid = deduce_face(pf,
                            next_pos[0], next_pos[1], next_pos[2],
                            D, last_vertex, last_edge
                        );
                        assert(next_fid != -1);
                    }
                    break;
                }
            }

            assert(hit);
            assert(next_fid != cur_fid);
            cur_fid = next_fid;
            cur_pos[0] = next_pos[0];
            cur_pos[1] = next_pos[1];
            cur_pos[2] = next_pos[2];
            walk++;
        }

        assert(found_tail);
        assert(walk <= max_walks);
    }

    write_edges("weids", weids);
    

    // TODO(Imad): project wpids on best-fit plane and jarvis march

    // BFS to get the smesh mpids
    std::queue<E_Int> Q;
    for (E_Int fid : wfids) Q.push(fid);

    printf("wfids before: %lu", wfids.size());

    while (!Q.empty()) {
        E_Int fid = Q.front();
        Q.pop();

        const auto &neis = F2F[fid];
        const auto &pe = F2E[fid];

        for (size_t i = 0; i < pe.size(); i++) {
            E_Int eid = pe[i];
            if (weids.find(eid) != weids.end()) continue;
            E_Int nei = neis[i];
            if (wfids.find(nei) == wfids.end()) {
                wfids.insert(nei);
                Q.push(nei);
            }
        }
    }

    write_ngon("wfids", wfids);

    printf(" - after: %lu\n", wfids.size());
    
    for (E_Int eid : weids) ewalls.insert(eid);
    for (E_Int fid : wfids) fwalls.insert(fid);
    
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
    Mf.make_point_faces();
    Mf.make_pnormals();
    Mf.make_point_edges();

    for (E_Int pid = 0; pid < Mf.np; pid++) {
        const auto &pe = Mf.P2E[pid];
        for (E_Int eid : pe) {
            const auto &e = Mf.E[eid];
            assert(e.p == pid || e.q == pid);
        }
    }

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
        spatches[i].make_fnormals();
        spatches[i].make_point_faces();
        spatches[i].make_pnormals();
        spatches[i].compute_min_distance_between_points();
        plocs.push_back(Mf.locate(spatches[i]));
        Mf.correct_near_points_and_edges(spatches[i], plocs[i]);
        printf("slave %d\n", i);
        mpatches.push_back(Mf.extract_bounding_smesh(spatches[i], plocs[i]));
    }

    Mf.write_edges("ewalls", ewalls);
    Mf.write_ngon("fwalls", fwalls);
    point_write("pchains", pchains);
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

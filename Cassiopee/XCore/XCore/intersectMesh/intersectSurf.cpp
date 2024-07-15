#include "xcore.h"
#include "karray.h"
#include "../common/common.h"
#include "mesh.h"
#include "smesh.h"
#include "dcel.h"
#include "vertex.h"
#include "face.h"
#include "hedge.h"
#include "io.h"
#include "cycle.h"

IMesh reconstruct_mesh(IMesh &M, const Dcel &D, E_Int color)
{
    // Isolate patch faces
    std::set<E_Int> pfset(M.patch);

    // Isolte patch points
    std::set<E_Int> ppset;

    for (E_Int face : pfset) {
        const auto &pn = M.F[face];
        for (E_Int p : pn) ppset.insert(p);
    }

    // Faces and point not belonging to pfset are renumbered first
    E_Int np = 0, nf = 0;

    std::map<E_Int, E_Int> new_fids;
    for (E_Int i = 0; i < M.nf; i++) {
        if (pfset.find(i) != pfset.end()) continue;

        new_fids[i] = nf++;
    }

    std::map<E_Int, E_Int> new_pids;
    std::map<E_Int, E_Int> new_ppids;

    for (E_Int i = 0; i < M.np; i++) {
        if (ppset.find(i) == ppset.end())
            new_pids[i] = np++;
    }

    //E_Int nop = np;

    for (E_Int p : ppset) new_ppids[p] = np++;

    assert(np == M.np);

    // Renumber ppoints and pfaces

    // Maps ofaces to locally indexed spawned faces
    std::map<E_Int, std::vector<E_Int>> ofid_to_ofids;

    std::map<E_Int, E_Int> new_ifids;

    std::map<E_Int, E_Int> new_ipids;

    for (size_t i = 0; i < D.F.size(); i++) {
        Face *f = D.F[i];

        if (D.C[i]->inout != Cycle::OUTER) continue;

        // Keep the faces whose oid of color is not -1

        E_Int ofid = f->oid[color];

        if (ofid == -1) continue;

        // ofid must have been a patch face
        assert(pfset.find(ofid) != pfset.end());

        // Renumber current iface
        new_ifids[i] = nf++;

        // Add current iface to oface spawn
        ofid_to_ofids[ofid].push_back(i);

        // Vertices of current iface must be kept
        auto vertices = Dcel::get_face_vertices(f);

        for (Vertex *v : vertices) {
            
            // Vertex could be a ppoint
            E_Int oid = v->oid[color];

            if (oid != -1) {
                assert(new_ppids.find(oid) != new_ppids.end());
                
                assert(ppset.find(oid) != ppset.end());


            } else {

                // This is a completely new point to M
                E_Int vid = v->id;
                assert(vid != -1);

                if (new_ipids.find(vid) == new_ipids.end()) {
                    new_ipids[vid] = np++;
                }
            }
        }
    }

    // All ofaces must be present in the map
    assert(ofid_to_ofids.size() == pfset.size());

    // Write untouched points
    std::vector<E_Float> new_X(np, -1), new_Y(np, -1), new_Z(np, -1);
    for (const auto &pt : new_pids) {
        E_Int opid = pt.first;
        E_Int npid = pt.second;


        assert(npid < M.np - (E_Int)ppset.size());

        assert(new_X[npid] == -1);
        assert(new_Y[npid] == -1);
        assert(new_Z[npid] == -1);

        new_X[npid] = M.X[opid];
        new_Y[npid] = M.Y[opid];
        new_Z[npid] = M.Z[opid];
    }

    // Add the ppids and ipids coordinates

    for (const auto &ppids : new_ppids) {
        E_Int opid = ppids.first;
        E_Int npid = ppids.second;

        assert(npid >= M.np - (E_Int)ppset.size() && npid < M.np);

        assert(new_X[npid] == -1);
        assert(new_Y[npid] == -1);
        assert(new_Z[npid] == -1);

        assert(opid < M.np);
        assert(npid < M.np);
        
        new_X[npid] = M.X[opid];
        new_Y[npid] = M.Y[opid];
        new_Z[npid] = M.Z[opid];
    }

    for (size_t i = 0; i < D.V.size(); i++) {
        assert((size_t)D.V[i]->id == i);
    }

    for (const auto &ipids : new_ipids) {
        E_Int opid = ipids.first;
        E_Int npid = ipids.second;
        
        assert(new_X[npid] == -1);
        assert(new_Y[npid] == -1);
        assert(new_Z[npid] == -1);

        Vertex *v = D.V[opid];

        new_X[npid] = v->x;
        new_Y[npid] = v->y;
        new_Z[npid] = 0;
    }

    // Add the kept faces

    std::vector<std::vector<E_Int>> new_F(nf);

    for (const auto &fids : new_fids) {
        E_Int ofid = fids.first;
        E_Int nfid = fids.second;


        assert(nfid < nf);

        const auto &pn = M.F[ofid];

        auto &new_face = new_F[nfid];

        for (E_Int opid : pn) {
            assert(opid < M.np);
            auto it = ppset.find(opid);

            if (it != ppset.end()) {
                new_face.push_back(new_ppids.at(opid));
            } else {
                new_face.push_back(new_pids.at(opid));
            }
        }
    }

    // Add the intersection faces

    for (const auto &fdata : ofid_to_ofids) {
        const auto &ofids = fdata.second;

        for (E_Int ofid : ofids) {
            Face *f = D.F[ofid];
            assert(f->oid[color] == fdata.first);

            auto vertices = Dcel::get_face_vertices(f);
            E_Int nfid = new_ifids.at(ofid);

            auto &new_face = new_F[nfid];
            assert(new_face.size() == 0);

            for (Vertex *v : vertices) {
                E_Int oid = v->oid[color];
               
                if (oid != -1) {
                    assert(ppset.find(oid) != ppset.end());
                    new_face.push_back(new_ppids.at(oid));
                } else {
                    assert(new_ipids.find(v->id) != new_ipids.end());
                    new_face.push_back(new_ipids.at(v->id));
                }
            }

        }
    }

    // Add the cell connectivity

    std::vector<std::vector<E_Int>> new_C(M.nc);

    for (E_Int i = 0; i < M.nc; i++) {
        const auto &pf = M.C[i];

        auto &new_cell = new_C[i];

        for (size_t j = 0; j < pf.size(); j++) {
            E_Int face = pf[j];

            auto it = ofid_to_ofids.find(face);

            if (it == ofid_to_ofids.end()) {
                assert(pfset.find(face) == pfset.end());
                new_cell.push_back(new_fids.at(face));
            } else {
                assert(pfset.find(face) != pfset.end());
                const auto &ifids = it->second;
                for (E_Int ifid : ifids) new_cell.push_back(new_ifids.at(ifid));
            }
        }
    }


    IMesh new_M;
    new_M.np = np;
    new_M.nf = nf;
    new_M.nc = M.nc;
    new_M.X = new_X;
    new_M.Y = new_Y;
    new_M.Z = new_Z;
    new_M.F = new_F;
    new_M.C = new_C;

    //assert(M.np == (E_Int)ppset.size() + nop);
    
    /*
    std::vector<E_Int> IP;
    for (E_Int i = M.np; i < np; i++) IP.push_back(i);

    point_write("IPOINTS", new_X.data(), new_Y.data(), new_Z.data(), IP);

    std::vector<E_Int> OP;
    for (E_Int i = nop; i < M.np; i++) OP.push_back(i);

    point_write("OPOINTS", new_X.data(), new_Y.data(), new_Z.data(), OP);
    */

    return new_M;
}

PyObject *K_XCORE::intersectSurf(PyObject *self, PyObject *args)
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

    for (E_Int i = 0; i < mpatch_size; i++) M.patch.insert(mpatch[i]);

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
    
    IMesh new_M = M.extract_conformized();
    IMesh new_S = S.extract_conformized();

    new_M.orient_skin(OUT);
    new_S.orient_skin(IN);

    // Extract surface meshes
    Smesh Mf(new_M);
    Smesh Sf(new_S);

    Dcel D(Mf, Sf);
 
    D.find_intersections();

    IMesh M_inter = reconstruct_mesh(new_M, D, Dcel::RED);
    
    IMesh S_inter = reconstruct_mesh(new_S, D, Dcel::BLACK);

    // Export
    PyObject *Mout = M_inter.export_karray();
    PyObject *Sout = S_inter.export_karray();

    PyObject *Out = PyList_New(0);
    PyList_Append(Out, Mout);
    PyList_Append(Out, Sout);
    Py_DECREF(Mout);
    Py_DECREF(Sout);
    Karray_free_ngon(marray);
    Karray_free_ngon(sarray);

    return Out;
}
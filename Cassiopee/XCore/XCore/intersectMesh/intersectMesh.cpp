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
#include "primitives.h"

static
IMesh reconstruct_mesh(IMesh &M, Smesh &Mf, const Dcel &D, E_Int color)
{
    // Isolate patch faces
    std::set<E_Int> pfset(M.patch);

    // Isolte patch points
    std::set<E_Int> ppset;

    for (E_Int face : pfset) {
        const auto &pn = M.F[face];
        for (E_Int p : pn) ppset.insert(p);
    }

    // Faces and points not belonging to pfset are renumbered first
    E_Int np = 0, nf = 0;

    std::map<E_Int, E_Int> new_fids;
    for (E_Int i = 0; i < M.nf; i++) {
        if (pfset.find(i) == pfset.end())
            new_fids[i] = nf++;
    }

    std::map<E_Int, E_Int> new_pids;
    std::map<E_Int, E_Int> new_ppids;

    for (E_Int i = 0; i < M.np; i++) {
        if (ppset.find(i) == ppset.end())
            new_pids[i] = np++;
    }

    for (E_Int p : ppset) new_ppids[p] = np++;

    assert(np == M.np);

    // Renumber ppoints and pfaces

    // Maps ofaces to locally indexed spawned faces
    std::map<E_Int, std::vector<E_Int>> ofid_to_ofids;

    std::map<E_Int, E_Int> new_ifids;

    std::map<E_Int, E_Int> new_ipids;

    std::set<E_Int> missing;
    for (E_Int i = 0; i < Mf.nf; i++) missing.insert(i);

    for (size_t i = 0; i < D.F.size(); i++) {
        Face *f = D.F[i];

        if (D.C[i]->inout != Cycle::OUTER) continue;

        // Keep the faces whose oid of color is not -1

        E_Int ofid = f->oid[color];

        if (ofid == -1) continue;

        missing.erase(ofid);

        ofid = Mf.l2gf.at(ofid);

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
        new_Z[npid] = v->z;
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
            //assert(f->oid[color] == fdata.first);

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
    new_M.ctag = M.ctag;

    return new_M;
}

PyObject *K_XCORE::intersectMesh(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVE, *SPATCH;
  
    if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVE, &SPATCH)) {
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
        return NULL;
    }

    // Init and orient master/slave meshes
    IMesh S(*sarray.cn, sarray.X, sarray.Y, sarray.Z, sarray.npts);

    M.make_skin();
    S.make_skin();

    M.orient_skin(OUT);
    S.orient_skin(IN);

    M.patch.clear();
    for (E_Int fid : M.skin) M.patch.insert(fid);

    printf("Master patch: %zu faces\n", M.patch.size());

    // Check slave intersection patch (zero-based)
    E_Int *spatch = NULL;
    E_Int spatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(SPATCH, spatch, spatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(sarray);
        RAISE("Bad slave patch.");
        return NULL;
    }

    printf("Slave patch: " SF_D_ " faces\n", spatch_size);

    for (E_Int i = 0; i < spatch_size; i++) S.patch.insert(spatch[i]-1);

    // Extract surface meshes
    Smesh Mf(M);
    Smesh Sf(S);
    
    //Mf.write_ngon("Mf");
    //Sf.write_ngon("Sf");

    puts("Making point edges...");
    Mf.make_point_edges();
    Sf.make_point_edges();

    puts("Making point faces...");
    Mf.make_point_faces_all();
    Sf.make_point_faces_all();

    puts("Making point/face normals...");
    Mf.make_pnormals();
    Sf.make_pnormals();

    puts("Hashing master faces...");
    Mf.make_bbox();
    Mf.hash_faces();

    puts("Initaliazing...");

    Dcel D(Mf, Sf);

    puts("Locating points...");

    D.locate_spoints(Mf, Sf);

    puts("Computing intersections...");

    D.find_intersections_3D(Mf, Sf);

    puts("Resolving hedges...");

    D.resolve_hedges(Mf, Sf);

    puts("Reconstructing meshes...");

    D.reconstruct(Mf, Sf);

    for (Vertex *v : D.V) {
        E_Int oid = v->oid[0];
        if (oid != -1) v->oid[0] = Mf.l2gp[oid];

        oid = v->oid[1];
        if (oid != -1) v->oid[1] = Sf.l2gp[oid];
    }
    
    M = reconstruct_mesh(M, Mf, D, Dcel::RED);

    IMesh S_inter = reconstruct_mesh(S, Sf, D, Dcel::BLACK);

    // Export
    printf("Exporting... ");

    PyObject *Sout = S_inter.export_karray();

    printf("Done.\n");

    Karray_free_ngon(sarray);

    Py_DECREF(SPATCH);

    return Sout;
}

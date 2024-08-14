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
IMesh reconstruct_mesh(IMesh &M, Smesh &Mf, const Dcel &D, Int color)
{
    // Isolate patch faces
    std::set<Int> pfset(M.patch);

    // Isolte patch points
    std::set<Int> ppset;

    for (Int face : pfset) {
        const auto &pn = M.F[face];
        for (Int p : pn) ppset.insert(p);
    }

    // Faces and points not belonging to pfset are renumbered first
    Int np = 0, nf = 0;

    std::map<Int, Int> new_fids;
    for (Int i = 0; i < M.nf; i++) {
        if (pfset.find(i) == pfset.end())
            new_fids[i] = nf++;
    }

    std::map<Int, Int> new_pids;
    std::map<Int, Int> new_ppids;

    for (Int i = 0; i < M.np; i++) {
        if (ppset.find(i) == ppset.end())
            new_pids[i] = np++;
    }

    for (Int p : ppset) new_ppids[p] = np++;

    assert(np == M.np);

    // Renumber ppoints and pfaces

    // Maps ofaces to locally indexed spawned faces
    std::map<Int, std::vector<Int>> ofid_to_ofids;

    std::map<Int, Int> new_ifids;

    std::map<Int, Int> new_ipids;

    std::set<Int> missing;
    for (Int i = 0; i < Mf.nf; i++) missing.insert(i);

    for (size_t i = 0; i < D.F.size(); i++) {
        Face *f = D.F[i];

        if (D.C[i]->inout != Cycle::OUTER) continue;

        // Keep the faces whose oid of color is not -1

        Int ofid = f->oid[color];

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
            Int oid = v->oid[color];

            if (oid != -1) {
                assert(new_ppids.find(oid) != new_ppids.end());
                
                assert(ppset.find(oid) != ppset.end());


            } else {

                // This is a completely new point to M
                Int vid = v->id;
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
    std::vector<Float> new_X(np, -1), new_Y(np, -1), new_Z(np, -1);
    for (const auto &pt : new_pids) {
        Int opid = pt.first;
        Int npid = pt.second;


        assert(npid < M.np - (Int)ppset.size());

        assert(new_X[npid] == -1);
        assert(new_Y[npid] == -1);
        assert(new_Z[npid] == -1);

        new_X[npid] = M.X[opid];
        new_Y[npid] = M.Y[opid];
        new_Z[npid] = M.Z[opid];
    }

    // Add the ppids and ipids coordinates

    for (const auto &ppids : new_ppids) {
        Int opid = ppids.first;
        Int npid = ppids.second;

        assert(npid >= M.np - (Int)ppset.size() && npid < M.np);

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
        Int opid = ipids.first;
        Int npid = ipids.second;
        
        assert(new_X[npid] == -1);
        assert(new_Y[npid] == -1);
        assert(new_Z[npid] == -1);

        Vertex *v = D.V[opid];

        new_X[npid] = v->x;
        new_Y[npid] = v->y;
        new_Z[npid] = v->z;
    }

    // Add the kept faces

    std::vector<std::vector<Int>> new_F(nf);

    for (const auto &fids : new_fids) {
        Int ofid = fids.first;
        Int nfid = fids.second;


        assert(nfid < nf);

        const auto &pn = M.F[ofid];

        auto &new_face = new_F[nfid];

        for (Int opid : pn) {
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

        for (Int ofid : ofids) {
            Face *f = D.F[ofid];
            //assert(f->oid[color] == fdata.first);

            auto vertices = Dcel::get_face_vertices(f);
            Int nfid = new_ifids.at(ofid);

            auto &new_face = new_F[nfid];
            assert(new_face.size() == 0);

            for (Vertex *v : vertices) {
                Int oid = v->oid[color];
               
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

    std::vector<std::vector<Int>> new_C(M.nc);

    for (Int i = 0; i < M.nc; i++) {
        const auto &pf = M.C[i];

        auto &new_cell = new_C[i];

        for (size_t j = 0; j < pf.size(); j++) {
            Int face = pf[j];

            auto it = ofid_to_ofids.find(face);

            if (it == ofid_to_ofids.end()) {
                assert(pfset.find(face) == pfset.end());
                new_cell.push_back(new_fids.at(face));
            } else {
                assert(pfset.find(face) != pfset.end());
                const auto &ifids = it->second;
                for (Int ifid : ifids) new_cell.push_back(new_ifids.at(ifid));
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

    //assert(M.np == (Int)ppset.size() + nop);
    
    /*
    std::vector<Int> IP;
    for (Int i = M.np; i < np; i++) IP.push_back(i);

    point_write("IPOINTS", new_X.data(), new_Y.data(), new_Z.data(), IP);

    std::vector<Int> OP;
    for (Int i = nop; i < M.np; i++) OP.push_back(i);

    point_write("OPOINTS", new_X.data(), new_Y.data(), new_Z.data(), OP);
    */

    return new_M;
}

PyObject *K_XCORE::intersectMesh(PyObject *self, PyObject *args)
{
    PyObject *MASTER, *SLAVE, *MPATCH, *SPATCH;
  
    if (!PYPARSETUPLE_(args, OOOO_, &MASTER, &SLAVE, &MPATCH, &SPATCH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray marray;
    Karray sarray;

    Int ret;

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

    // Check master intersection patch (zero-based)
    Int *mpatch = NULL;
    Int mpatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(MPATCH, mpatch, mpatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad master patch.");
        return NULL;
    }

    printf("Master patch: " SF_D_ " faces\n", mpatch_size);

    // Check slave intersection patch (zero-based)
    Int *spatch = NULL;
    Int spatch_size = -1;
    ret = K_NUMPY::getFromNumpyArray(SPATCH, spatch, spatch_size, true);
    if (ret != 1) {
        Karray_free_ngon(marray);
        Karray_free_ngon(sarray);
        RAISE("Bad slave patch.");
        return NULL;
    }

    printf("Slave patch: " SF_D_ " faces\n", spatch_size);

    for (Int i = 0; i < mpatch_size; i++) M.patch.insert(mpatch[i]);
    for (Int i = 0; i < spatch_size; i++) S.patch.insert(spatch[i]);

    M.orient_skin(OUT);
    S.orient_skin(IN);

    // Extract surface meshes
    Smesh Mf(M);
    Smesh Sf(S);
    
    //for (Float &z : Mf.Z) z = 0;
    //for (Float &z : Sf.Z) z = 0;

    Mf.write_ngon("Mf");
    Sf.write_ngon("Sf");

    Mf.make_point_edges();
    Sf.make_point_edges();

    Mf.make_point_faces_all();
    Sf.make_point_faces_all();

    Mf.make_pnormals();
    Sf.make_pnormals();

    Dcel D(Mf, Sf);

    D.locate_spoints(Mf, Sf);

    D.init_Cp(Mf, Sf);

    D.find_intersections_3D(Mf, Sf);

    D.resolve_hedges(Mf, Sf);

    D.reconstruct();

    for (Vertex *v : D.V) {
        Int oid = v->oid[0];
        if (oid != -1) v->oid[0] = Mf.l2gp[oid];

        oid = v->oid[1];
        if (oid != -1) v->oid[1] = Sf.l2gp[oid];
    }

    IMesh S_inter = reconstruct_mesh(S, Sf, D, Dcel::BLACK);

    IMesh M_inter = reconstruct_mesh(M, Mf, D, Dcel::RED);
    

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

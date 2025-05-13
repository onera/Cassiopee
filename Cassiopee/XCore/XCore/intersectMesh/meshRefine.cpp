/*    
    Copyright 2013-2025 Onera.

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
#include <stack>

#include "mesh.h"
#include "common/common.h"

#include "triangle.h"
#include "io.h"
#include "primitives.h"

E_Int meshes_mutual_refinement(IMesh &M, IMesh &S)
{
    size_t refM = 0, refS = 0;
    E_Int iter = 0;

    S.init_adaptation_data();
    M.init_adaptation_data();

    E_Int S_np_before = S.np;

    do {
        iter++;

        // Refine M wrt S

        S.make_point_faces();
        M.make_bbox();
        M.hash_skin(); // TODO(Imad): hash_patch!
    
        refM = M.refine(S);
        printf("Refined mf: %zu\n", refM);

        // Refine S wrt M
        /*if (iter == 1 || (iter > 1 && refM > 0)) {
            M.make_point_faces();
            S.make_bbox();
            S.hash_skin(); // TODO(Imad): hash_patch!

            refS = S.refine_slave(M);
            printf("Refined sf: %zu\n", refS);
        }*/
        
    } while (refS > 0);

    S.make_point_faces();

    // Project all the new points from S onto M faces
    // TODO(Imad): shouldn't this step be done after conformizing?

    for (E_Int spid = S_np_before; spid < S.np; spid++) {
        const auto &pf = S.P2F[spid];
        std::vector<E_Int> stids;
        for (auto stid : pf) {
            if (S.F[stid].size() == 3) stids.push_back(stid);
        }

        // Compute the normal at spid := sum of the normals of stids
        E_Float N[3] = {0};
        for (auto stid : stids) {
            const auto &pn = S.F[stid];
            E_Int a = pn[0], b = pn[1], c = pn[2];
            E_Float v0[3] = {S.X[b]-S.X[a], S.Y[b]-S.Y[a], S.Z[b]-S.Z[a]};
            E_Float v1[3] = {S.X[c]-S.X[a], S.Y[c]-S.Y[a], S.Z[c]-S.Z[a]};
            E_Float fN[3];
            K_MATH::cross(v0, v1, fN);
            N[0] += fN[0];
            N[1] += fN[1];
            N[2] += fN[2];
        }

        E_Float NORM = K_MATH::norm(N, 3);
        assert(Sign(NORM) != 0);
        N[0] /= NORM;
        N[1] /= NORM;
        N[2] /= NORM;

        TriangleIntersection TI;
        E_Int hit = 0;

        if (M.is_point_inside(S.X[spid], S.Y[spid], S.Z[spid])) {
            hit = M.project_point(S.X[spid], S.Y[spid], S.Z[spid],
                -N[0], -N[1], -N[2], TI, spid - S_np_before);
        } else {
            hit = M.project_point(S.X[spid], S.Y[spid], S.Z[spid],
                N[0], N[1], N[2], TI, spid - S_np_before);
        }

        assert(hit);
        assert(M.patch.find(TI.face) != M.patch.end());

        printf("Spid: %f %f %f -> Proj: %f %f %f (t = %f)\n",
            S.X[spid], S.Y[spid], S.Z[spid], TI.x, TI.y, TI.z, TI.t);
        
        // Replace the point by its projection
        S.X[spid] = TI.x;
        S.Y[spid] = TI.y;
        S.Z[spid] = TI.z;
    }

    return 0;
}

struct Fidn {
    E_Int fid;
    E_Float N[3];
};

size_t IMesh::refine_slave(const IMesh &master)
{
    const auto &mpatch = master.patch;
    std::set<E_Int> mpids;
    for (const auto mfid : mpatch) {
        assert(master.face_is_active(mfid));
        const auto &pn = master.F[mfid];
        for (const auto p : pn) mpids.insert(p);
    }

    std::vector<Fidn> spatch;
    spatch.reserve(patch.size());
    for (const auto fid : patch) {
        const auto &pn = F[fid];
        assert(face_is_tri(fid));
        assert(face_is_active(fid));
        E_Int a = pn[0], b = pn[1], c = pn[2];
        E_Float v0[3] = {X[b]-X[a], Y[b]-Y[a], Z[b]-Z[a]};
        E_Float v1[3] = {X[c]-X[a], Y[c]-Y[a], Z[c]-Z[a]};
        Fidn fidn;
        fidn.fid = fid;
        K_MATH::cross(v0, v1, fidn.N);
        E_Float NORM = K_MATH::norm(fidn.N, 3);
        assert(Sign(NORM) != 0);
        fidn.N[0] /= NORM;
        fidn.N[1] /= NORM;
        fidn.N[2] /= NORM;
        spatch.push_back(fidn);
    }

    const auto &mX = master.X;
    const auto &mY = master.Y;
    const auto &mZ = master.Z;

    // Master points to Slave triangles
    std::map<E_Int, std::vector<E_Int>> mpids_to_stids;

    // TODO(Imad): how to accelerate this?
    for (const auto &fidn : spatch) {
        const E_Int fid = fidn.fid;
        const E_Float *N = fidn.N;
        const auto &pn = F[fid];
        E_Int a = pn[0];
        E_Int b = pn[1];
        E_Int c = pn[2];

        for (auto mpid : mpids) {
            // Does the projection of mpid along N live in the triangle fid?

            E_Float V[3] = {mX[mpid]-X[a], mY[mpid]-Y[a], mZ[mpid]-Z[a]};
            E_Float dp = K_MATH::dot(V, N, 3);

            E_Float Proj[3];
            Proj[0] = mX[mpid] - dp * N[0];
            Proj[1] = mY[mpid] - dp * N[1];
            Proj[2] = mZ[mpid] - dp * N[2];

            bool inside = Triangle::is_point_inside(Proj[0], Proj[1], Proj[2],
                X[a], Y[a], Z[a],
                X[b], Y[b], Z[b],
                X[c], Y[c], Z[c]);
            
            if (inside) {
                mpids_to_stids[mpid].push_back(fid);
            }
        }
    }

    // Invert the map: S triangles to all the M points (projection) within them
    std::map<E_Int, std::vector<E_Int>> stids_to_mpids;

    for (const auto &m2s : mpids_to_stids) {
        E_Int mpid = m2s.first;
        const auto &stids = m2s.second;
        for (const auto tid : stids) {
            stids_to_mpids[tid].push_back(mpid);
        }
    }

    // Keep the S triangles containing three of more M points
    std::map<E_Int, std::vector<E_Int>> filtered_stids_map;

    for (const auto &s2m : stids_to_mpids) {
        if (s2m.second.size() >= 3) {
            filtered_stids_map.insert({s2m.first, s2m.second});
        }
    }

    // Sensor
    std::map<E_Int, std::vector<E_Int>> ref_stids_to_mtids;

    for (const auto &s2m : filtered_stids_map) {
        auto stid = s2m.first;
        const auto &mpids = s2m.second;

        // For every sface, how many points are contained within mface?
        std::map<E_Int, size_t> mcontained;

        for (auto mpid : mpids) {
            // M tris sharing mpid
            const auto &mtris = master.P2F[mpid];
 
            for (auto mtri : mtris) {
                assert(master.face_is_active(mtri));

                // M tri should belong to patch
                if (master.patch.find(mtri) == master.patch.end()) continue;

                // S tri contains another point of M tri
                mcontained[mtri]++;
            }
        }

        for (const auto &mfdat : mcontained) {
            E_Int mtid = mfdat.first;

            // The number of points of mtid contained within stid...
            size_t npts = mfdat.second;

            // ... should be at most equal to the stride of mtid!
            assert(npts <= master.F[mtid].size());

            // Discard mtid if it is not completely enclosed by stid
            if (npts < master.F[mtid].size()) continue;

            // Discard sface if it is a duplicate of mface
            if (faces_are_dups(stid, mtid, master)) continue;

            // Add mtid to the set of faces enclosed by stid
            ref_stids_to_mtids[stid].push_back(mtid);
        }
    }

    //printf("Faces to refine: %zu\n", ref_mfaces_to_sfaces.size());

    E_Int iter = 0;
    size_t ret = 0;

    while (!ref_stids_to_mtids.empty()) {
        iter++;

        auto ref_data = smooth_ref_data(ref_stids_to_mtids);

        auto ref_faces = prepare_for_refinement(ref_data);

        ret += ref_faces.size();

        refine_faces(ref_faces);

        std::map<E_Int, std::vector<E_Int>> new_sensor;

        for (const auto &fdata : ref_stids_to_mtids) {
            E_Int parent = fdata.first;
            const auto &mtids = fdata.second;

            const auto &children = fchildren[parent];

            for (E_Int child : children) {
                for (E_Int mtid : mtids) {
                    if (face_contains_sface(child, mtid, master)) {
                        // TODO(Imad): I think we can safely discard this check
                        assert(!faces_are_dups(child, mtid, master));
                        new_sensor[child].push_back(mtid);
                    }
                }
            }

            // Update mpatch
            patch.erase(parent);
            for (E_Int child : children) patch.insert(child);
        }

        ref_stids_to_mtids = new_sensor;
    }

    return ret;
}

size_t IMesh::refine(const IMesh &S)
{
    // Hash mpatch
    hash_patch();

    // Isolate spatch points
    const auto &spatch = S.patch;
    std::set<E_Int> spoints;

    for (E_Int sface : spatch) {
        assert(S.face_is_active(sface));
        const auto &spn = S.F[sface];
        for (E_Int sp : spn)
            spoints.insert(sp);
    }

    // Locate spatch points within mpatch
    std::map<E_Int, std::vector<E_Int>> spoints_to_mfaces;


    for (E_Int spt : spoints) {
        E_Int voxel_x = floor((S.X[spt] - xmin) / HX);
        E_Int voxel_y = floor((S.Y[spt] - ymin) / HY);
        E_Int voxel_z = floor((S.Z[spt] - zmin) / HZ);
        E_Int spt_bin = voxel_x + NX * voxel_y + NXY * voxel_z;

        const auto &pf = bin_faces[spt_bin];

        assert(pf.size() > 0);

        for (E_Int fid : pf) {

            const auto &pn = F[fid];

            assert(pn.size() == 3 || pn.size() == 4);

            E_Int a = pn[0], b = pn[1], c = pn[2];

            if (Triangle::is_point_inside(S.X[spt], S.Y[spt], S.Z[spt],
                X[a], Y[a], Z[a], X[b], Y[b], Z[b], X[c], Y[c], Z[c])) {
                
                spoints_to_mfaces[spt].push_back(fid);
            } else if (pn.size() == 4) {

                E_Int d = pn[3];

                if (Triangle::is_point_inside(S.X[spt], S.Y[spt], S.Z[spt],
                    X[c], Y[c], Z[c], X[d], Y[d], Z[d], X[a], Y[a], Z[a])) {
                
                    spoints_to_mfaces[spt].push_back(fid);
                }
            }
        }
    }

    // Mface to spoints located within it
    std::map<E_Int, std::vector<E_Int>> mfpoints;

    for (const auto &data : spoints_to_mfaces) {
        E_Int spt = data.first;
        const auto &mfaces = data.second;

        for (const E_Int mf : mfaces) {
            mfpoints[mf].push_back(spt);
        }
    }

    // Keep the mfaces which contain 3 points or more
    std::map<E_Int, std::vector<E_Int>> filtered_mfaces_map;

    for (const auto &mfpts : mfpoints) {
        if (mfpts.second.size() >= 3) {
            filtered_mfaces_map.insert({mfpts.first, mfpts.second});
        }
    }

    // Sensor
    std::map<E_Int, std::vector<E_Int>> ref_mfaces_to_sfaces;

    for (const auto &mfdata : filtered_mfaces_map) {
        E_Int mface = mfdata.first;
        const auto &spoints = mfdata.second;

        // For every sface, how many points are contained within mface?
        std::map<E_Int, size_t> scontained;

        for (E_Int sp : spoints) {
            // Sfaces sharing sp
            const auto &sfaces = S.P2F[sp];
 
            for (E_Int sface : sfaces) {
                assert(S.face_is_active(sface));

                // Sface should belong to spatch
                if (spatch.find(sface) == spatch.end()) continue;

                // mface contains another point of sface
                scontained[sface]++;
            }
        }

        for (const auto &sfdat : scontained) {
            E_Int sface = sfdat.first;

            // The number of points of sface contained within mface...
            size_t npts = sfdat.second;

            // ... should be at most equal to the stride of the sface!
            assert(npts <= S.F[sface].size());

            // Discard sface if it is not completely enclosed by mface
            if (npts < S.F[sface].size()) continue;

            // Discard sface if it is a duplicate of mface
            if (faces_are_dups(mface, sface, S)) continue;

            // Add sface to the set of faces enclosed by mface
            ref_mfaces_to_sfaces[mface].push_back(sface);
        }
    }

    //printf("Faces to refine: %zu\n", ref_mfaces_to_sfaces.size());

    E_Int iter = 0;
    size_t ret = 0;

    while (!ref_mfaces_to_sfaces.empty()) {
        iter++;

        auto ref_data = smooth_ref_data(ref_mfaces_to_sfaces);

        auto ref_faces = prepare_for_refinement(ref_data);

        ret += ref_faces.size();

        refine_faces(ref_faces);

        std::map<E_Int, std::vector<E_Int>> new_sensor;

        for (const auto &fdata : ref_mfaces_to_sfaces) {
            E_Int parent = fdata.first;
            const auto &sfaces = fdata.second;

            const auto &children = fchildren[parent];

            for (E_Int child : children) {
                for (E_Int sface : sfaces) {
                    if (face_contains_sface(child, sface, S)) {
                        assert(!faces_are_dups(child, sface, S));
                        new_sensor[child].push_back(sface);
                    }
                }
            }

            // Update mpatch
            patch.erase(parent);
            for (E_Int child : children) patch.insert(child);
        }

        ref_mfaces_to_sfaces = new_sensor;
    }

    return ret;
}


std::vector<E_Int> IMesh::smooth_ref_data(
    const std::map<E_Int, std::vector<E_Int>> &sensor)
{
    // TODO(Imad): shouldn't the size be factive.size() instead of nf?
    std::vector<E_Int> ref_data(nf, 0);
    std::stack<E_Int> stk;

    for (const auto &fdata : sensor) {
        assert(face_is_active(fdata.first));
        ref_data[fdata.first] = 1; // Warning: size must be nf for this to be ok
        stk.push(fdata.first);
    }

    /*
    while (!stk.empty()) {
        E_Int face = stk.top();
        stk.pop();

        E_Int face_incr = ref_data[face] + flevel[face];

        auto neis = get_active_neighbours(face);

        for (auto nei : neis) {
            assert(nei) != -1;

            E_Int nei_incr = ref_data[nei] + flevel[nei];

            E_Int diff = abs(face_incr - nei_incr);

            if (diff <= 1) continue;

            E_Int fid = face_incr > nei_incr ? nei : face;

            ref_data[fid] += 1;

            stk.push(fid);
        }
    }
    */

    return ref_data;
}

std::vector<E_Int> IMesh::prepare_for_refinement(const std::vector<E_Int> &ref_data)
{
    std::vector<E_Int> ref_faces;
    
    for (size_t i = 0; i < ref_data.size(); i++) {
        if (ref_data[i] > 0) ref_faces.push_back(i);
    }

    // Refine the lower-level faces first
    std::sort(ref_faces.begin(), ref_faces.end(),
        [&] (E_Int i, E_Int j) { return flevel[i] < flevel[j]; });
    
    // Resize data structures
    resize_point_data(ref_faces.size());
    resize_face_data(ref_faces.size());

    return ref_faces;
}

/*
std::vector<E_Int> IMesh::get_active_neighbours(E_Int face)
{
    return std::vector<E_Int>();
}
*/

void IMesh::refine_faces(const std::vector<E_Int> &ref_faces)
{
    for (E_Int ref_face : ref_faces) {
        if (face_is_tri(ref_face)) refine_tri(ref_face);
        else refine_quad(ref_face);
    }
}

void IMesh::resize_point_data(size_t nref_faces)
{
    size_t nnew_points = np + nref_faces * 5;
    X.resize(nnew_points);
    Y.resize(nnew_points);
    Z.resize(nnew_points);
}

void IMesh::resize_face_data(size_t nref_faces)
{
    size_t nnew_faces = nf + nref_faces * 4;
    F.resize(nnew_faces);
    flevel.resize(nnew_faces, -1);
}

void IMesh::refine_tri(E_Int tri)
{
    // Refine the edges
    const auto &pn = F[tri];
    E_Int ec[3];

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int p = pn[i];
        E_Int q = pn[(i+1)%pn.size()];
        UEdge edge(p, q);

        auto it = ecenter.find(edge);

        if (it == ecenter.end()) {
            X[np] = 0.5 * (X[p] + X[q]);
            Y[np] = 0.5 * (Y[p] + Y[q]);
            Z[np] = 0.5 * (Z[p] + Z[q]);
            ecenter[edge] = np;
            ec[i] = np;
            np++;
        } else {
            ec[i] = it->second;
        }
    }

    // New face points
    E_Int nf0 = nf, nf1 = nf+1, nf2 = nf+2, nf3 = nf+3;

    F[nf0] = { pn[0], ec[0], ec[2] };
    F[nf1] = { ec[0], pn[1], ec[1] };
    F[nf2] = { ec[2], ec[1], pn[2] };
    F[nf3] = { ec[0], ec[1], ec[2] };

    // Disable quad and enable its children
    factive.erase(tri);
    factive.insert(nf0);
    factive.insert(nf1);
    factive.insert(nf2);
    factive.insert(nf3);

    // Set quad children pointers
    fchildren[tri] = { nf0, nf1, nf2, nf3 };
    flevel[nf0] = flevel[nf1] = flevel[nf2] = flevel[nf3] = flevel[tri] + 1;

    nf += 4;
}

void IMesh::refine_quad(E_Int quad)
{
    // Refine the edges
    const auto &pn = F[quad];
    E_Int ec[4];

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int p = pn[i];
        E_Int q = pn[(i+1)%pn.size()];
        UEdge edge(p, q);

        auto it = ecenter.find(edge);

        if (it == ecenter.end()) {
            X[np] = 0.5 * (X[p] + X[q]);
            Y[np] = 0.5 * (Y[p] + Y[q]);
            Z[np] = 0.5 * (Z[p] + Z[q]);
            ecenter[edge] = np;
            ec[i] = np;
            np++;
        } else {
            ec[i] = it->second;
        }
    }

    // Face centroid

    X[np] = Y[np] = Z[np] = 0.0;
    for (E_Int i = 0; i < 4; i++) {
        E_Int p = pn[i];
        X[np] += X[p];
        Y[np] += Y[p];
        Z[np] += Z[p];
    }
    X[np] *= 0.25;
    Y[np] *= 0.25;
    Z[np] *= 0.25;

    // New face points
    E_Int nf0 = nf, nf1 = nf+1, nf2 = nf+2, nf3 = nf+3;

    F[nf0] = { pn[0], ec[0], np,    ec[3] };
    F[nf1] = { ec[0], pn[1], ec[1], np    };
    F[nf2] = { np   , ec[1], pn[2], ec[2] };
    F[nf3] = { ec[3], np   , ec[2], pn[3] };

    // Disable quad and enable its children
    factive.erase(quad);
    factive.insert(nf0);
    factive.insert(nf1);
    factive.insert(nf2);
    factive.insert(nf3);

    // Set quad children pointers
    fchildren[quad] = { nf0, nf1, nf2, nf3 };
    flevel[nf0] = flevel[nf1] = flevel[nf2] = flevel[nf3] = flevel[quad] + 1;

    np += 1;
    nf += 4;
}

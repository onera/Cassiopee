#include "icapsule.h"
#include "ray.h"
#include "BVH.h"
#include "io.h"
#include "primitives.h"

// Returns a list of all the intersections, forwards and backwards
// Up to the caller to parse the data
void Smesh::ray_BVH_intersect(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz, BVH_node *node,
    std::vector<PointLoc> &plocs) const
{
    bool hit = ray_AABB_intersect(ox, oy, oz, dx, dy, dz, node->box);
    if (!hit) return;

    if (!node->left && !node->right) {
        for (E_Int i = node->start; i < node->end; i++) {
            E_Int fid = bvh_indices[i];
            const auto &pn = Fc[fid];
            const E_Float *fc = &fcenters[3*fid];

            for (size_t j = 0; j < pn.size(); j++) {
                E_Int p = pn[j];
                E_Int q = pn[(j+1)%pn.size()];

                E_Float u, v, w, t, x, y, z;
                
                bool hit = MollerTrumboreAnyDir(
                    ox, oy, oz, dx, dy, dz,
                    X[p], Y[p], Z[p],
                    X[q], Y[q], Z[q],
                    fc[0], fc[1], fc[2],
                    u, v, w, t, x, y, z
                );
                
                if (hit) {
                    PointLoc ploc;
                    ploc.fid = fid;
                    ploc.sub = j;
                    ploc.bcrd[0] = u;
                    ploc.bcrd[1] = v;
                    ploc.bcrd[2] = w;
                    ploc.t = t;
                    ploc.x = x;
                    ploc.y = y;
                    ploc.z = z;

                    // on p
                    if      (Sign(1-u, NEAR_VERTEX_TOL) == 0)
                        ploc.v_idx = j;
                    // on q
                    else if (Sign(1-v, NEAR_VERTEX_TOL) == 0)
                        ploc.v_idx = (j+1)%pn.size();
                    // on edge {p, q}
                    else if (Sign(w, NEAR_EDGE_TOL) == 0)
                        ploc.e_idx = j;

                    plocs.push_back(ploc);
                    
                    break;
                }
            }
        }
        return;
    }

    ray_BVH_intersect(ox, oy, oz, dx, dy, dz, node->left, plocs);
    ray_BVH_intersect(ox, oy, oz, dx, dy, dz, node->right, plocs);
}

std::vector<Point> projections;

void Smesh::project(const Smesh &Mf, const std::vector<E_Int> &mpids,
    std::vector<PointLoc> &plocs) const
{
    for (size_t i = 0; i < mpids.size(); i++) {
        E_Int mpid = mpids[i];
        const E_Float *N = &Mf.pnormals[3*mpid];
        E_Float mx = Mf.X[mpid];
        E_Float my = Mf.Y[mpid];
        E_Float mz = Mf.Z[mpid];
        std::vector<PointLoc> mlocs; // to parse
        ray_BVH_intersect(mx, my, mz, N[0], N[1], N[2], bvh_root, mlocs);
        auto &ploc = plocs[i];
        E_Float min_abs_t = EFLOATMAX;
        for (const auto &mloc : mlocs) {
            if (fabs(mloc.t) < min_abs_t) {
                min_abs_t = fabs(mloc.t);
                ploc = mloc;
            }
        }
        if (ploc.fid != -1) {
            projections.push_back({ploc.x, ploc.y, ploc.z});
        }
    }
}

std::vector<E_Int> Smesh::deduce_ref_faces(const std::vector<E_Int> &mpids,
    const std::vector<PointLoc> &plocs_m, const Smesh &Mf,
    std::vector<E_Int> &ref_faces)
{
    // Sfids to enclosed mpids
    std::map<E_Int, std::vector<E_Int>> smap;
    for (size_t i = 0; i < mpids.size(); i++) {
        const auto &ploc_m = plocs_m[i];
        if (ploc_m.fid == -1) continue;

        std::vector<E_Int> sfids;
        E_Int dummy;
        get_shared_faces(ploc_m, sfids, dummy, dummy);

        for (E_Int sfid : sfids)
            smap[sfid].push_back(mpids[i]);
    }

    std::map<E_Int, std::vector<E_Int>> sfid_to_mfids;

    for (const auto &spdata : smap) {
        E_Int sfid = spdata.first;

        // Partly enclosed mfids to point count
        std::map<E_Int, E_Int> mfid_pcount;

        const auto &mpids = spdata.second;
        for (E_Int mpid : mpids) {
            const auto &pf = Mf.P2F[mpid];
            for (E_Int mfid : pf) mfid_pcount[mfid] += 1;
        }

        // Deduce the mfids that are completely enclosed
        for (const auto &mpdata : mfid_pcount) {
            E_Int mfid = mpdata.first;
            size_t count = mpdata.second;
            assert(count <= Mf.Fc[mfid].size());
            if (count == Mf.Fc[mfid].size()) {
                sfid_to_mfids[sfid].push_back(mfid);
            }
        }
    }

    ref_faces.clear();
    ref_faces.reserve(sfid_to_mfids.size());
    for (auto it = sfid_to_mfids.begin(); it != sfid_to_mfids.end(); it++)
        ref_faces.push_back(it->first);

    return ref_faces;
}

void ICapsule::refine(Smesh &Mf, std::set<E_Int> &mfids, Smesh &Sf,
    std::vector<PointLoc> &plocs_s)
{
    size_t ref_M, ref_S;
    
    E_Int iter = 0;

    do {
        iter++;

        std::vector<E_Int> fat_sfids;
        std::vector<E_Int> fat_mfids;
        
        Sf.make_BVH();
        Mf.make_BVH();

        /*********************** Sf refinement ***********************/

        // Isolate the points to project
        std::set<E_Int> mpids_set;
        for (E_Int fid : mfids) {
            const auto &pn = Mf.F[fid];
            for (E_Int p : pn) mpids_set.insert(p);
        }
        std::vector<E_Int> mpids;
        for (E_Int p : mpids_set) mpids.push_back(p);

        // Project mpids on Sf
        std::vector<PointLoc> plocs_m(mpids.size());
        Sf.project(Mf, mpids, plocs_m);

        // Deduce sfids to refine
        std::vector<E_Int> sref_faces;
        Sf.deduce_ref_faces(mpids, plocs_m, Mf, sref_faces);
        printf("Fat sfids: %lu\n", sref_faces.size());

        // Refine
        ref_S = sref_faces.size();
        if (ref_S > 0) {
            Sf.refine(sref_faces);
            Sf.conformize();
            Sf.make_pnormals();
        }


        /*********************** Mf refinement ***********************/

        // Test all the Sf points
        std::vector<E_Int> spids;
        spids.reserve(Sf.np);
        for (E_Int i = 0; i < Sf.np; i++) spids.push_back(i);

        // Reproject spids on mfaces
        plocs_s.resize(Sf.np);
        Mf.project(Sf, spids, plocs_s);
        Sf.replace_by_projections(spids, plocs_s);
        
        // Deduce mfids to refine
        std::vector<E_Int> mref_faces;
        Mf.deduce_ref_faces(spids, plocs_s, Sf, mref_faces);
        printf("Fat mfids: %lu\n", mref_faces.size());

        ref_M = mref_faces.size();
        if (ref_M > 0) {
            Mf.refine(mref_faces);
            Mf.conformize();
            Mf.make_pnormals();
            // update mfids
            for (E_Int fparent : mref_faces) {
                const auto &children = Mf.fchildren[fparent].back();
                for (E_Int child : children) mfids.insert(child);
            }
        }

        Mf.destroy_BVH(Mf.bvh_root);
        Sf.destroy_BVH(Sf.bvh_root);
    } while (ref_M > 0 || ref_S > 0);

    Mf.hash_faces();
    Sf.hash_faces();

    //point_write("projections", projections);
}

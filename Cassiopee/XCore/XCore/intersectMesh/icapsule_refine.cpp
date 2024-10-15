#include "icapsule.h"
#include "ray.h"
#include "BVH.h"
#include "io.h"

// Returns a list of all the intersections, forwards and backwards
// Up to the caller to parse the data
void Smesh::ray_BVH_intersect(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz, BVH_node *node,
    std::vector<PointLoc> &plocs)
{
    bool hit = ray_AABB_intersect(ox, oy, oz, dx, dy, dz, node->box);
    if (!hit) return;

    if (!node->left && !node->right) {
        for (E_Int i = node->start; i < node->end; i++) {
            E_Int fid = bvh_indices[i];
            const auto &pn = F[fid];
            E_Int a = pn[0], b = pn[1], c = pn[2];
            E_Float u, v, w, t, x, y, z;
            bool hit = MollerTrumboreAnyDir(ox, oy, oz, dx, dy, dz,
                X[a], Y[a], Z[a], X[b], Y[b], Z[b], X[c], Y[c], Z[c],
                u, v, w, t, x, y, z);
            if (hit) {
                PointLoc ploc;
                ploc.fid = fid;
                ploc.bcrd[0] = u;
                ploc.bcrd[1] = v;
                ploc.bcrd[2] = w;
                ploc.t = t;
                ploc.x = x;
                ploc.y = y;
                ploc.z = z;
                plocs.push_back(ploc);
            }
        }
        return;
    }

    ray_BVH_intersect(ox, oy, oz, dx, dy, dz, node->left, plocs);
    ray_BVH_intersect(ox, oy, oz, dx, dy, dz, node->right, plocs);
}

std::vector<Point> projections;

void Smesh::project(const Smesh &Mf, const std::vector<E_Int> &mpids,
    std::vector<PointLoc> &plocs)
{
    for (size_t i = 0; i < mpids.size(); i++) {
        E_Int mpid = mpids[i];
        assert(mpid >= 0);
        assert(mpid < Mf.np);
        assert(Mf.pnormals.size() == (size_t)(3*Mf.np));
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
                ploc.fid = mloc.fid;
                ploc.t = mloc.t;
                ploc.x = mloc.x;
                ploc.y = mloc.y;
                ploc.z = mloc.z;
            }
        }
        if (ploc.fid != -1) {
            projections.push_back({ploc.x, ploc.y, ploc.z});
        }
    }
}

static
std::vector<E_Int> deduce_ref_faces(const std::vector<E_Int> &mpids,
    const std::vector<PointLoc> &plocs_m, const Smesh &Mf,
    std::vector<E_Int> &ref_faces)
{
    // Sfids to enclosed mpids
    std::map<E_Int, std::vector<E_Int>> smap;
    for (size_t i = 0; i < mpids.size(); i++) {
        const auto &ploc_m = plocs_m[i];
        if (ploc_m.fid == -1) continue;
        smap[ploc_m.fid].push_back(mpids[i]);
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

        /*********************** Sf refinement ***********************/

        // Construct the BVH of Sf
        Sf.make_BVH();

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
        Sf.destroy_BVH(Sf.bvh_root);

        // Deduce sfids to refine
        std::vector<E_Int> sref_faces;
        deduce_ref_faces(mpids, plocs_m, Mf, sref_faces);

        for (const E_Int sfid : sref_faces) {
            fat_sfids.push_back(sfid);
        }
        printf("Fat sfids: %lu\n", fat_sfids.size());

        ref_S = sref_faces.size();
        if (ref_S > 0) {
            Sf.refine(sref_faces);
            Sf.conformize();
            Sf.hash_faces();
            Sf.make_pnormals();
        }

        /*********************** Mf refinement ***********************/

        // Test all the Sf points
        std::vector<E_Int> spids;
        spids.reserve(Sf.np);
        for (E_Int i = 0; i < Sf.np; i++) spids.push_back(i);

        // Deduce mfids to refine
        std::vector<E_Int> mref_faces;
        deduce_ref_faces(spids, plocs_s, Sf, mref_faces);

        for (const E_Int mfid : mref_faces) {
            fat_mfids.push_back(mfid);
        }
        printf("Fat mfids: %lu\n", fat_mfids.size());

        ref_M = mref_faces.size();
        if (ref_M > 0) {
            Mf.refine(mref_faces);
            Mf.conformize();
            Mf.hash_faces();
            Mf.make_pnormals();
            // update mfids
            for (E_Int fparent : mref_faces) {
                const auto &children = Mf.fchildren[fparent].back();
                for (E_Int child : children) mfids.insert(child);
            }
            plocs_s = Mf.locate(Sf);
        }
    } while (ref_M > 0 || ref_S > 0);

    point_write("projections", projections);
}
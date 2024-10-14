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
        const E_Float *N = &Mf.pnormals[3*mpid];
        E_Float mx = Mf.X[mpid];
        E_Float my = Mf.Y[mpid];
        E_Float mz = Mf.Z[mpid];
        std::vector<PointLoc> mlocs; // to parse
        ray_BVH_intersect(mx, my, mz, N[0], N[1], N[2], bvh_root, mlocs);
        auto &ploc = plocs[i];
        E_Float min_abs_t = EFLOATMAX;
        for (const auto &mloc : mlocs) {
            if (fabs(mloc.t) < EFLOATMAX) {
                min_abs_t = mloc.t;
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

void ICapsule::refine(Smesh &Mf, const std::set<E_Int> &mfids, Smesh &Sf,
    std::vector<PointLoc> &plocs_s)
{
    E_Int ref_M, ref_S;
    do {
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
        std::vector<PointLoc> plocs_m(mpids.size());

        // Project mpids on Sf
        Sf.project(Mf, mpids, plocs_m);

        ref_M = 0;
        ref_S = 0;

    } while (ref_M > 0 || ref_S > 0);

    point_write("projections", projections);
}
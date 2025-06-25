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
#include "icapsule.h"
#include "ray.h"
#include "BVH.h"
#include "io.h"
#include "primitives.h"

std::vector<PointLoc> Smesh::project(const Smesh &Mf,
    const std::vector<E_Int> &mpids) const
{
    std::vector<PointLoc> plocs;
    plocs.reserve(mpids.size());

    for (size_t i = 0; i < mpids.size(); i++) {
        E_Int mpid = mpids[i];
        const E_Float *N = &Mf.pnormals[3*mpid];
        E_Float mx = Mf.X[mpid];
        E_Float my = Mf.Y[mpid];
        E_Float mz = Mf.Z[mpid];
        std::vector<PointLoc> mlocs; // to parse
        ray_intersect_BVH(mx, my, mz, N[0], N[1], N[2], root_node_idx, mlocs);
        PointLoc ploc;
        E_Float min_abs_t = EFLOATMAX;
        for (const auto &mloc : mlocs) {
            if (fabs(mloc.t) < min_abs_t) {
                min_abs_t = fabs(mloc.t);
                ploc = mloc;
            }
        }
        plocs.push_back(ploc);
    }

    return plocs;
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

std::vector<PointLoc> ICapsule::refine(Smesh &Mf, std::set<E_Int> &mfids,
    Smesh &Sf)
{
    size_t ref_M, ref_S;
    
    E_Int iter = 0;
    
    std::vector<PointLoc> plocs_s;

    do {
        iter++;

        std::vector<E_Int> fat_sfids;
        std::vector<E_Int> fat_mfids;

        Mf.make_fcenters();
        Mf.make_BVH();

        /*********************** Mf refinement ***********************/

        // Test all the Sf points
        std::vector<E_Int> spids;
        spids.reserve(Sf.np);
        for (E_Int i = 0; i < Sf.np; i++) spids.push_back(i);

        // Reproject spids on mfaces
        Mf.make_fcenters();
        plocs_s = Mf.locate2(Sf);
        Sf.replace_by_projections(spids, plocs_s);

        // Deduce mfids to refine
        std::vector<E_Int> mref_faces;
        Mf.deduce_ref_faces(spids, plocs_s, Sf, mref_faces);
        printf("Fat mfids: %zu\n", mref_faces.size());

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

        /*********************** Sf refinement ***********************/

        // Isolate the points to project
        std::set<E_Int> mpids_set;
        for (E_Int fid : mfids) {
            const auto &pn = Mf.Fc[fid];
            for (E_Int p : pn) mpids_set.insert(p);
        }
        std::vector<E_Int> mpids;
        for (E_Int p : mpids_set) mpids.push_back(p);

        // Project mpids on Sf
        Sf.make_fcenters();
        Sf.make_BVH();
        auto plocs_m = Sf.project(Mf, mpids);

        // Deduce sfids to refine
        std::vector<E_Int> sref_faces;
        Sf.deduce_ref_faces(mpids, plocs_m, Mf, sref_faces);
        printf("Fat sfids: %zu\n", sref_faces.size());

        // Refine
        ref_S = sref_faces.size();
        if (ref_S > 0) {
            Sf.refine(sref_faces);
            Sf.conformize();
            Sf.make_pnormals();
        }

    } while (ref_M > 0 || ref_S > 0);

    return plocs_s;
}

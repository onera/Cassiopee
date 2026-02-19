/*    
    Copyright 2013-2026 ONERA.

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
#include <queue>
#include "smesh.h"
#include "primitives.h"
#include "ray.h"
#include "io.h"

void Smesh::get_shared_faces(const PointLoc &loc, std::vector<E_Int> &ret,
    E_Int &pid, E_Int &eid) const
{
    ret.clear();

    E_Int fid = loc.fid;
    assert(fid != -1);

    // TODO(Imad): e_idx is no good!!!!

    if (loc.e_idx != -1) {
        assert(loc.v_idx == -1);
        const auto &pe = F2E[fid];
        eid = pe[loc.e_idx];
        const auto &pf = E2F[eid];
        assert(pf[0] == fid || pf[1] == fid);
        ret.push_back(pf[0]);
        // O could be on a boundary edge
        if (pf[1] != -1) ret.push_back(pf[1]);
    }
    else if (loc.v_idx != -1) 
    {
        assert(loc.e_idx == -1);
        const auto &pn = Fc[fid];
        pid = pn[loc.v_idx];
        const auto &pf = P2F[pid];
        // For consistency
#ifndef NDEBUG
        bool found_fid = false;
#endif
        for (auto face : pf) 
        {
            ret.push_back(face);
#ifndef NDEBUG
            if (face == fid)
            {
                found_fid = true;
            }
#endif
        }
        assert(found_fid == true);
    }
    else 
    {
        ret.push_back(fid);
    }
}

std::set<E_Int> ewalls;
std::set<E_Int> fwalls;
std::vector<Point> pchains;

std::set<E_Int> Smesh::extract_covering_faces(const Smesh &Sf,
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

    /*
    // Sort the pchain counterclockwise
    E_Int a = pchain[0], b = pchain[1], c = pchain[2];
    E_Float ux = Sf.X[b] - Sf.X[a];
    E_Float uy = Sf.Y[b] - Sf.Y[a];
    E_Float uz = Sf.Z[b] - Sf.Z[a];
    E_Float vx = Sf.X[c] - Sf.X[b];
    E_Float vy = Sf.Y[c] - Sf.Y[b];
    E_Float vz = Sf.Z[c] - Sf.Z[b];
    E_Float cp[3] = {uy*vz - uz*vy, uz*vx - ux*vz, ux*vy - uy*vx};
    // TODO(Imad): inherit fnormals
    const E_Float *N_b = &fnormals[3*plocs[b].fid];
    E_Float dp = K_MATH::dot(cp, N_b, 3);
    E_Int cmp = Sign(dp);
    assert(cmp != 0);
    if (cmp < 0) std::reverse(pchain.begin(), pchain.end());
    */

    //Sf.write_points("pchain.im", pchain);

    for (E_Int p : pchain) pchains.push_back({Sf.X[p], Sf.Y[p], Sf.Z[p]});

    std::set<E_Int> wfids;
    std::set<E_Int> weids;
 
    for (size_t i = 0; i < pchain.size(); i++) {
        E_Int p = pchain[i];
        E_Int q = pchain[(i+1)%pchain.size()];

        E_Float px = Sf.X[p], py = Sf.Y[p], pz = Sf.Z[p];
        E_Float qx = Sf.X[q], qy = Sf.Y[q], qz = Sf.Z[q];

        //point_write("p", px, py, pz);
        //point_write("q", qx, qy, qz);

        E_Float D[3] = {qx-px, qy-py, qz-pz};
        E_Float NORM = K_MATH::norm(D, 3);
        D[0] /= NORM, D[1] /= NORM, D[2] /= NORM;

        std::vector<E_Int> orig_faces;
        std::vector<E_Int> tail_faces;

        E_Int last_vertex = -1, last_edge = -1, dummy;

        get_shared_faces(plocs[p], orig_faces, last_vertex, last_edge); 
        get_shared_faces(plocs[q], tail_faces, dummy, dummy);

        //write_ngon("shared", orig_faces);

        E_Int starting_face = deduce_face(orig_faces, px, py, pz,
            D, last_vertex, last_edge, dummy);
        assert(starting_face != -1);

        bool found_tail = false;
        E_Int cur_fid = starting_face;
        E_Float cur_pos[3] = {px, py, pz};

        E_Int walk = 0;
        E_Int max_walks = 100;

        while (!found_tail && walk <= max_walks) {
            
            wfids.insert(cur_fid);

            //write_face("cur_fid", cur_fid);

            E_Float proj[3];
            get_unit_projected_direction(cur_fid, D, proj);

            const auto &pn = Fc[cur_fid];
            const auto &pe = F2E[cur_fid];
            assert(pe.size() == pn.size());
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

            //write_edges("wall", weids);

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
                    if (s > TOL && s < 1 - TOL) 
                    {
                        const auto &pe = F2E[cur_fid];
                        E_Int eid = pe[i];
                        last_edge = eid;
                        last_vertex = -1;
                        assert(E2F[eid][0] == cur_fid || E2F[eid][1] == cur_fid);
                        if (E2F[eid][0] == cur_fid) next_fid = E2F[eid][1];
                        else next_fid = E2F[eid][0];

                        next_pos[0] = cur_pos[0] + t * proj[0];
                        next_pos[1] = cur_pos[1] + t * proj[1];
                        next_pos[2] = cur_pos[2] + t * proj[2];
                    } 
                    else 
                    {
                        bool hit_p = (s <= TOL);
                        #ifndef NDEBUG
                        bool hit_q = (s >= 1 - TOL);
                        #endif
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
                            D, last_vertex, last_edge, dummy);
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

    //write_edges("wall", weids);

    // TODO(Imad): project wpids on best-fit plane and jarvis march

    // BFS to get the smesh mpids
    std::queue<E_Int> Q;
    for (E_Int fid : wfids) Q.push(fid);

    //write_ngon("fwall_before_BFS", wfids);

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

    //write_ngon("fwall_after_BFS", wfids);

    //write_ngon("fwall", wfids);
    //write_edges("ewall", weids);

    for (E_Int eid : weids) ewalls.insert(eid);
    for (E_Int fid : wfids) fwalls.insert(fid);
    
    return wfids;
}

Smesh Smesh::extract_smesh(const std::set<E_Int> &fids, bool check_Euler)
{
    assert(0);
    return Smesh();
}

Smesh Smesh::extract_conformized()
{
    assert(0);
    return Smesh();
}


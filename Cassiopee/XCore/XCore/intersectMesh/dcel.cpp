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
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "dcel.h"
#include "primitives.h"
#include "io.h"
#include "smesh.h"
#include "triangle.h"

void Dcel::init_vertices(const Smesh &Mf, const Smesh &Sf,
    const std::vector<PointLoc> &plocs)
{
    E_Int duplicate_vertices = 0;

    for (E_Int i = 0; i < Mf.np; i++) {
        Vertex *v = new Vertex(Mf.X[i], Mf.Y[i], Mf.Z[i]);
        v->oids[Dcel::RED] = i;
        vertex_set.insert(v);
    }

    for (E_Int i = 0; i < Sf.np; i++) {
        Vertex tmp(Sf.X[i], Sf.Y[i], Sf.Z[i]);
        auto it = vertex_set.find(&tmp);
        Vertex *v = NULL;
        if (it == vertex_set.end()) {
            v = new Vertex(Sf.X[i], Sf.Y[i], Sf.Z[i]);
            v->oids[Dcel::BLACK] = i;
            vertex_set.insert(v);
        } else {
            v = *it;
            v->oids[Dcel::BLACK] = i;
            duplicate_vertices++;
        }
        v->ploc = plocs[i];
    }

    V.reserve(vertex_set.size());
    for (Vertex *v : vertex_set) {
        v->id = V.size();
        V.push_back(v);
    }

    printf("Duplicate vertices: %d\n", duplicate_vertices);
}

std::vector<Dcel::Vertex *> Dcel::get_face_vertices(const Face *f) const
{
    assert(f);
    std::vector<Vertex *> ret;
    Hedge *h = f->rep;
    ret.push_back(h->orig);
    Hedge *w = h->next;
    while (w != h) {
        ret.push_back(w->orig);
        w = w->next;
    }
    return ret;
}

Smesh Dcel::reconstruct(const Smesh &Mf, int color, bool check_Euler) const
{
    Smesh ret;
    ret.check_Euler = check_Euler;
    auto &new_F = ret.F;
    std::map<Vertex *, E_Int> new_pids;
    ret.np = ret.nf = 0;

    std::vector<Face *> fids;

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];
        Hedge *h = f->rep;
        if (h->color != color) continue;
        Hedge *w = h->next;
        while (w != h) {
            if (w->color != h->color) break;
            w = w->next;
        }
        if (w == h) {
            fids.push_back(f);
        }
    }

    char fname[128] = {0};
    sprintf(fname, "single_color_%d.im", color);
    write_ngon(fname, fids);

    for (Face *f : F) {
        if (f->oids[color] == -1) continue;

        std::vector<E_Int> pn;

        std::vector<Vertex *> vertices = get_face_vertices(f);
        for (Vertex *v : vertices) {
            auto it = new_pids.find(v);
            if (it == new_pids.end()) {
                new_pids[v] = ret.np;
                pn.push_back(ret.np);
                ret.np++;
            } else {
                pn.push_back(it->second);
            }
        }

        new_F.push_back(pn);
        ret.nf++;
    }

    auto &new_X = ret.X;
    auto &new_Y = ret.Y;
    auto &new_Z = ret.Z;

    new_X.resize(ret.np), new_Y.resize(ret.np), new_Z.resize(ret.np);
    for (const auto &vdat : new_pids) {
        new_X[vdat.second] = vdat.first->x;
        new_Y[vdat.second] = vdat.first->y;
        new_Z[vdat.second] = vdat.first->z;
    }

    ret.Fc = ret.F;

    ret.make_edges();

    return ret;
}

Dcel::Dcel(const Smesh &Mf, const Smesh &Sf, const std::vector<PointLoc> &plocs)
{
    init_vertices(Mf, Sf, plocs);

    init_hedges_and_faces(Mf, RED);
    init_hedges_and_faces(Sf, BLACK);

    if (check_hedges(H) != 0) {
        fprintf(stderr, "Dcel: Inconsistent half-edge records!\n");
        abort();
    }
    
    make_cycles();
    set_cycles_inout();

    if (check_faces(H, F) != 0) {
        fprintf(stderr, "Dcel: Inconsistent face records!\n");
        abort();
    }


    std::vector<Point> xpoints;

    std::vector<Point> dpoints;

    // Register the intersections between Sf points and Mf edges
    E_Int v_on_e = 0;
    for (Vertex *v : V) {
        // Strictly the Sf points
        if (v->oids[0] == -1 && v->oids[1] != -1) {
            const auto &ploc = v->ploc;
            if (ploc.e_idx == -1) continue;

            // Set the reference edge for normal computation
            E_Int fid = ploc.fid;
            const auto &pe = Mf.F2E[fid];
            E_Int e = pe[ploc.e_idx];

            // Add the vertex to hedge intersections
            Hedge *h = H[2*e];
            assert(h->color == RED);
            Vertex *O = h->orig;
            Vertex *T = h->twin->orig;
            if (cmp_vtx(O, T) > 0) h = h->twin;
            hedge_intersections[h].push_back(v);

            dpoints.push_back({v->x,v->y,v->z});

            v_on_e++;
        }
    }
    printf("Sf points on Mf edges: %d\n", v_on_e);
    //point_write("dpoints.im", dpoints);

    // Trace
    for (E_Int eid_s = 0; eid_s < Sf.ne; eid_s++) {
        Hedge *hs = H[2*(Mf.ne + eid_s)];
        assert(hs->color == BLACK);
        Vertex *O = hs->orig;
        Vertex *T = hs->twin->orig;
        if (cmp_vtx(O, T) > 0) hs = hs->twin;

        E_Int p = hs->orig->oids[1];
        E_Int q = hs->twin->orig->oids[1];

        E_Float spx = Sf.X[p], spy = Sf.Y[p], spz = Sf.Z[p];
        E_Float sqx = Sf.X[q], sqy = Sf.Y[q], sqz = Sf.Z[q];

        E_Float D[3] = {sqx-spx, sqy-spy, sqz-spz};
        E_Float NORM = K_MATH::norm(D, 3);
        D[0] /= NORM, D[1] /= NORM, D[2] /= NORM;

        std::vector<E_Int> orig_faces;
        std::vector<E_Int> tail_faces;

        E_Int last_vertex = -1, last_edge = -1, dummy;

        Mf.get_shared_faces(plocs[p], orig_faces, last_vertex, last_edge); 
        Mf.get_shared_faces(plocs[q], tail_faces, dummy, dummy);

        //if (eid_s == 133) {
        //    point_write("O.im", O->x, O->y, O->z);
        //    point_write("T.im", T->x, T->y, T->z);
        //    Mf.write_ngon("ofaces.im", orig_faces);
        //    Mf.write_ngon("tfaces.im", tail_faces);
        //}

        E_Int starting_face = Mf.deduce_face(orig_faces, spx, spy, spz,
            D, last_vertex, last_edge, eid_s);
        assert(starting_face != -1);

        bool found_tail = false;
        E_Int cur_fid = starting_face;
        E_Float cur_pos[3] = {spx, spy, spz};

        E_Int walk = 0;
        E_Int max_walks = 20;

        while (!found_tail && walk <= max_walks) {

            for (auto fid : tail_faces) {
                if (fid == cur_fid) {
                    found_tail = true;
                    break;
                }
            }

            if (found_tail) break;

            E_Float proj[3];
            Mf.get_unit_projected_direction(cur_fid, D, proj);

            const auto &pn = Mf.Fc[cur_fid];
            const auto &pe = Mf.F2E[cur_fid];
            assert(pe.size() == pn.size());
            //const E_Float *fN = &Mf.fnormals[3*cur_fid];

            E_Int next_fid = -1;
            E_Float next_pos[3] = {EFLOATMAX, EFLOATMAX, EFLOATMAX};

            bool hit = false;

            for (size_t i = 0; i < pn.size(); i++) {
                E_Int p = pn[i];
                E_Int q = pn[(i+1)%pn.size()];
                E_Int e = pe[i];

                if (p == last_vertex || q == last_vertex || e == last_edge)
                    continue;
                
                E_Float px = Mf.X[p], py = Mf.Y[p], pz = Mf.Z[p];
                E_Float qx = Mf.X[q], qy = Mf.Y[q], qz = Mf.Z[q];
            
                E_Float t, s;
                hit = ray_edge_intersect(
                    cur_pos[0], cur_pos[1], cur_pos[2],
                    proj[0], proj[1], proj[2],
                    px, py, pz, qx, qy, qz,
                    t, s
                );

                if (hit) {
                    if (s > TOL && s < 1 - TOL) {
                        // Hit edge middle
                        const auto &pe = Mf.F2E[cur_fid];
                        E_Int eid_m = pe[i];
                        last_edge = eid_m;
                        last_vertex = -1;
                        assert(Mf.E2F[eid_m][0] == cur_fid || Mf.E2F[eid_m][1] == cur_fid);
                        if (Mf.E2F[eid_m][0] == cur_fid) next_fid = Mf.E2F[eid_m][1];
                        else next_fid = Mf.E2F[eid_m][0];

                        next_pos[0] = cur_pos[0] + t * proj[0];
                        next_pos[1] = cur_pos[1] + t * proj[1];
                        next_pos[2] = cur_pos[2] + t * proj[2];

                        // Create a new intersection vertex
                        Vertex tmp(next_pos[0], next_pos[1], next_pos[2]);
                        assert(vertex_set.find(&tmp) == vertex_set.end());
                        Vertex *x = new Vertex(next_pos[0], next_pos[1],
                            next_pos[2]);
                        
                        x->ploc.fid = cur_fid;
                        x->ploc.e_idx = i;

                        x->id = V.size();
                        V.push_back(x);
                        vertex_set.insert(x);

                        // Register the intersection
                        Hedge *hm = H[2*eid_m];
                        assert(hm->color == RED);
                        Vertex *O = hm->orig;
                        Vertex *T = hm->twin->orig;
                        if (cmp_vtx(O, T) > 0) hm = hm->twin;

                        hedge_intersections[hm].push_back(x);
                        hedge_intersections[hs].push_back(x);
                    } else {
                        // Hit an edge endpoint
                        bool hit_p = (s <= TOL);
                        bool hit_q = (s >= 1 - TOL);
                        assert(!(hit_p && hit_q));
                        last_edge = -1;
                        if (hit_p) last_vertex = p;
                        else last_vertex = q;
                        next_pos[0] = Mf.X[last_vertex];
                        next_pos[1] = Mf.Y[last_vertex];
                        next_pos[2] = Mf.Z[last_vertex];
                        const auto &pf = Mf.P2F[last_vertex];
                        next_fid = Mf.deduce_face(pf,
                            next_pos[0], next_pos[1], next_pos[2],
                            D, last_vertex, last_edge, eid_s
                        );
                        assert(next_fid != -1);

                        // Find Vertex corresponding to hit vertex
                        Vertex tmp(Mf.X[last_vertex], Mf.Y[last_vertex], Mf.Z[last_vertex]);
                        auto it = vertex_set.find(&tmp);
                        assert(it != vertex_set.end());
                        Vertex *x = *it;

                        // Register the intersection
                        hedge_intersections[hs].push_back(x);

                        xpoints.push_back({x->x, x->y, x->z});
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

    point_write("xpoints.im", xpoints);

    // Cut
    for (auto &h2x : hedge_intersections) {
        Hedge *h = h2x.first;
        auto it = hedge_intersections.find(h->twin);
        assert(it == hedge_intersections.end());

        auto &xs = h2x.second;

        Vertex *o = h->orig;
        Vertex *tail = h->twin->orig;
        assert(cmp_vtx(o, tail) < 0);

        // TODO(Imad): check that the intersections are 'sufficiently' spaced out
        for (Vertex *x : xs) {
            E_Float D[3] = {x->x-o->x, x->y-o->y, x->z-o->z};
            x->d2 = K_MATH::dot(D, D, 3);
        }

        std::sort(xs.begin(), xs.end(), [&] (const Vertex *a, const Vertex *b)
        {
            assert(Sign(a->d2-b->d2) != 0);
            return a->d2 < b->d2;
        });

        Hedge *current_h = h;
        Hedge *t = h->twin;

        for (Vertex *x : xs) {
            Hedge *e1 = new Hedge(x, current_h->color);
            Hedge *e2 = new Hedge(x, t->color);

            H.push_back(e1);
            H.push_back(e2);

            e1->left = current_h->left;
            e2->left = t->left;

            current_h->twin = e2;
            e2->twin = current_h;
            t->twin = e1;
            e1->twin = t;

            current_h->next = e1;
            e1->prev = current_h;
            t->next = e2;
            e2->prev = t;

            current_h = e1;
        }
    }

    // Resolve
    std::vector<std::vector<Hedge *>> list(V.size());
    for (Hedge *h : H) {
        Vertex *o = h->orig;
        list[o->id].push_back(h);
    }

    for (size_t vid = 0; vid < V.size(); vid++) {
        Vertex *v = V[vid];

        E_Float N[3] = {0, 0, 0};
        
        if (v->oids[0] != -1) {
            const E_Float *pN = &Mf.pnormals[3*v->oids[0]];
            for (E_Int i = 0; i < 3; i++) N[i] = pN[i];
        } else if (v->oids[1] != -1) {
            const E_Float *fN = &Mf.fnormals[3*v->ploc.fid];
            for (E_Int i = 0; i < 3; i++) N[i] = fN[i];
        } else {
            E_Int fid_m = v->ploc.fid;
            E_Int eid_m = Mf.F2E[fid_m][v->ploc.e_idx];
            const auto &pf = Mf.E2F[eid_m];
            assert(pf[0] == fid_m || pf[1] == fid_m);
            const E_Float *fN1 = &Mf.fnormals[3*pf[0]];
            const E_Float *fN2 = &Mf.fnormals[3*pf[1]];
            for (int i = 0; i < 3; i++) N[i] += fN1[i] + fN2[i];
            E_Float NORM = K_MATH::norm(N, 3);
            for (int i = 0; i < 3; i++) N[i] /= NORM;
        }

        auto &leaving = list[vid];

        sort_leaving_hedges(leaving, N);

        for (size_t i = 0; i < leaving.size(); i++) {
            Hedge *h = leaving[i];
            Hedge *w = leaving[(i+1)%leaving.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }
    }

    if (check_hedges(H) != 0) {
        fprintf(stderr, "Dcel: Inconsistent half-edge records!\n");
        assert(0);
        abort();
    }

    make_cycles();
    set_cycles_inout();

    write_hole_cycles("hole.im");
    write_degen_cycles("degen.im");
    write_inner_cycles("inner.im");

    auto new_F = make_cycle_faces(C);
    set_face_labels(new_F);
    update_hedge_faces(new_F);
    check_faces(H, new_F);

    for (Face *f : F) delete f;
    F = new_F;

    puts("ok");
}

void Dcel::update_hedge_faces(std::vector<Face *> &new_F)
{
    for (Hedge *h : H) h->left = NULL;

    for (Face *f : new_F) {
        assert(f);

        Hedge *h = f->rep;
        assert(h);
        Cycle *c = h->cycle;
        assert(c->inout == Cycle::INNER);

        h->left = f;
        assert(f->rep->left == f);

        Hedge *w = h->next;
        while (w != h) {
            assert(w->left == NULL);
            w->left = f;
            w = w->next;
        }
    }

    for (Face *f : new_F) {
        assert(f->rep->left == f);
    }

    for (Hedge *h : H) {
        Face *f = h->left;
        Hedge *w = h->next;
        while (w != h) {
            assert(w->left == f);
            w = w->next;
        }
    }
}

E_Int Dcel::check_faces(const std::vector<Hedge *> &H,
    const std::vector<Face *> &F)
{
    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];
        if (h->prev->left != h->left) { assert(0); return 1; }
        if (h->next->left != h->left) { assert(0); return 1; }
    }

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];
        //Hedge *h = f->rep;
        //Cycle *c = h->cycle;
        //assert(c->inout == Cycle::DEGEN || c->inout == Cycle::INNER);
        if (f->rep->left != f) { assert(0); return 1; }
    }

    puts("CHECK: FACES OK.");

    return 0;
}



void Dcel::set_face_labels(std::vector<Face *> &new_F)
{
    // Label each face with the ids of the original faces containing it

    for (size_t i = 0; i < new_F.size(); i++) {
        Face *f = new_F[i];

        // Get the first RED and BLACK half-edges in the face cycle
        Hedge *h = f->rep;

        Cycle *c = h->cycle;
        assert(c->inout == Cycle::INNER);

        Hedge *R = NULL;
        Hedge *B = NULL;
        E_Int RB = 0;

        if (h->color == Dcel::RED) {
            R = h;
            B = get_hedge_of_color(f, Dcel::BLACK);
            if (B) RB = 1;
        } else if (h->color == Dcel::BLACK) {
            B = h;
            R = get_hedge_of_color(f, Dcel::RED);
            if (R) RB = 1;
        } else {
            assert(0);
        }

        if (RB) {
            // RED face should always exist
            assert(R->left);
            f->oids[Dcel::RED] = R->left->oids[Dcel::RED];
            assert(f->oids[Dcel::RED] != -1);

            // BLACK face might not exist
            f->oids[Dcel::BLACK] = (B->left) ? B->left->oids[Dcel::BLACK] : -1;
        } else {
            // Only single color possible is RED
            if (!R && B) {
                write_hedge("black.im", B);
                point_write("orig.im", B->orig->x, B->orig->y, B->orig->z);
            }
            assert(R && !B);
            //Hedge *REP = (R != NULL) ? R : B;
            //if (REP != R) {
            //    hedge_write("black", REP);
            //}
            f->oids[Dcel::RED] = R->left->oids[Dcel::RED];
        }
    }
}

Dcel::Hedge *Dcel::get_hedge_of_color(Face *f, int color)
{
    Hedge *h = f->rep;
    if (h->color == color) return h;
    Hedge *w = h->next;
    while (w != h) {
        if (w->color == color) return w;
        w = w->next;
    }
    return NULL;
}

std::vector<Dcel::Face *> Dcel::make_cycle_faces(const std::vector<Cycle *> &C)
{
    std::vector<Face *> new_F;

    for (Cycle *c : C) {
        //if (c->inout == Cycle::HOLE) continue;
        if (c->inout != Cycle::INNER) continue;
        
        // Create a face record
        Face *f = new Face;

        // Set its rep hedge to some edge of the cycle
        Hedge *h = c->rep;
        f->rep = h;

        new_F.push_back(f);
    }

    return new_F;
}

void Dcel::init_hedges_and_faces(const Smesh &Mf, int color)
{
    printf("Doing color %d\n", color);

    size_t current_nh = H.size();

    H.reserve(current_nh + Mf.ne*2);

    std::vector<std::vector<Hedge *>> list(Mf.np);

    // Create hedge records

    const auto &E = Mf.E;

    for (E_Int i = 0; i < Mf.ne; i++) {
        const auto &e = E[i];

        E_Int p = e.p;
        E_Int q = e.q;

        Hedge *h = NULL;
        Hedge *t = NULL;

        {
            Vertex tmp(Mf.X[p], Mf.Y[p], Mf.Z[p]);
            auto it = vertex_set.find(&tmp);
            assert(it != vertex_set.end());
            Vertex *P = *it;
            h = new Hedge(P, color);
            list[p].push_back(h);
        }

        {
            Vertex tmp(Mf.X[q], Mf.Y[q], Mf.Z[q]);
            auto it = vertex_set.find(&tmp);
            assert(it != vertex_set.end());
            Vertex *Q = *it;
            t = new Hedge(Q, color);
            list[q].push_back(t);
        }

        h->twin = t;
        t->twin = h;

        H.push_back(h);
        H.push_back(t);
    }
    
    // Pair-up hedges

    const auto &pnormals = Mf.pnormals;

    for (E_Int pid = 0; pid < Mf.np; pid++) {
        auto &hedges = list[pid];

        const E_Float *N = &pnormals[3*pid];
        sort_leaving_hedges(hedges, N);

        for (size_t i = 0; i < hedges.size(); i++) {
            Hedge *h = hedges[i];
            Hedge *w = hedges[(i+1)%hedges.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }
    }

    // Create face records

    const auto &F2E = Mf.F2E;
    const auto &E2F = Mf.E2F;

    size_t current_nf = F.size();

    F.reserve(current_nf + Mf.nf);

    for (E_Int fid = 0; fid < Mf.nf; fid++) {
        const auto &pe = F2E[fid];
        E_Int first_edge = pe[0];
        E_Int where = 2*first_edge;
        Hedge *h = H[current_nh + where];
        Hedge *t = H[current_nh + where + 1];
        assert(h->twin == t);
        assert(t->twin == h);

        Face *f = new Face;
        f->oids[color] = fid;

        // Face must lie to the left of hedge
        assert(E2F[first_edge][0] == fid || E2F[first_edge][1] == fid);
        Hedge *REP = (E2F[first_edge][0] == fid) ? h : t;

        assert(REP->left == NULL);

        f->rep = REP;
        REP->left = f;
        Hedge *w = REP->next;
        while (w != REP) { w->left = f; w = w->next; }
        
        F.push_back(f);
    }
}

void Dcel::sort_leaving_hedges(std::vector<Hedge *> &leaving,
    const E_Float N[3]) const
{
    // Choose a vector that is not parallel to N

    E_Float ref_vec[3] = {0, N[2], -N[1]};
    E_Float NORM = K_MATH::norm(ref_vec, 3);
    if (Sign(NORM) == 0) {
        ref_vec[0] = -N[2];
        ref_vec[1] = 0;
        ref_vec[2] =  N[0];
        NORM = K_MATH::norm(ref_vec, 3);
        assert(Sign(NORM) != 0);
    }

    E_Float dp = K_MATH::dot(ref_vec, N, 3);
    for (E_Int i = 0; i < 3; i++) ref_vec[i] = ref_vec[i] - dp * N[i];
    NORM = K_MATH::norm(ref_vec, 3);
    for (E_Int i = 0; i < 3; i++) ref_vec[i] /= NORM;

    std::vector<E_Float> angles(leaving.size());

    for (size_t i = 0; i < leaving.size(); i++) {
        Hedge *h = leaving[i];
        Hedge *t = h->twin;

        Vertex *P = h->orig;
        Vertex *Q = t->orig;
        assert(P != Q);

        // Project the hedge onto the plane (pid, N)
        E_Float PQ[3] = {Q->x-P->x, Q->y-P->y, Q->z-P->z};
        E_Float dp = K_MATH::dot(PQ, N, 3);
        E_Float PQ_proj[3];
        for (E_Int j = 0; j < 3; j++) PQ_proj[j] = PQ[j] - dp * N[j];

        E_Float costheta = K_MATH::dot(ref_vec, PQ_proj, 3) / K_MATH::norm(PQ_proj, 3);
        costheta = std::min(costheta, 1.0);
        costheta = std::max(costheta, -1.0);
        E_Float angle = acos(costheta);
        
        // Determine the direction of the angle
        E_Float C[3];
        K_MATH::cross(ref_vec, PQ_proj, C);
        
        if (K_MATH::dot(N, C, 3) > 0) angle = 2*K_MATH::PI - angle;
        
        angles[i] = angle;
    }

    std::vector<E_Int> indices(leaving.size());
    for (size_t i = 0; i < leaving.size(); i++) indices[i] = i;

    std::sort(indices.begin(), indices.end(), [&](E_Int i, E_Int j)
    {
        if (angles[i] < angles[j]) return true;
        
        if (angles[i] > angles[j]) return false;
        
        Hedge *h = leaving[i];
        Hedge *w = leaving[j];

        assert(h->color != w->color);

        Vertex *O_h = h->orig;
        Vertex *T_h = h->twin->orig;

        Vertex *O_w = w->orig;
        Vertex *T_w = w->twin->orig;

        assert(O_w == O_h);
        assert(T_w == T_h);

        // If this is the origin: red comes first
        if (cmp_vtx(O_w, T_w) < 0) {
            if (h->color == RED) return true;
            return false;
        } else {
            if (h->color == RED) return false;
            return true;
        }
    });
    
    std::vector<Hedge *> tmp(leaving);
    for (size_t i = 0; i < leaving.size(); i++) leaving[i] = tmp[indices[i]];
}

E_Int Dcel::check_hedges(const std::vector<Hedge *> &H)
{
    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];
        if (h->prev->next != h) { assert(0); return 1; }
        if (h->next->prev != h) { assert(0); return 1; }
        if (h->twin->twin != h) { assert(0); return 1; }
        if (h->twin->next->orig != h->orig) { assert(0); return 1; }
        if (h->prev->twin->orig != h->orig) { assert(0); return 1; }
    }

    puts("CHECK: EDGES OK.");

    return 0;
}

void Dcel::make_cycles()
{
    C.clear();
    for (Hedge *h : H) h->cycle = NULL;

    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];

        if (h->cycle) continue;

        Cycle *c = new Cycle(h);
        C.push_back(c);

        h->cycle = c;

        Hedge *w = h->next;
        while (w != h) {
            w->cycle = c;
            w = w->next;
        }
    }
}

void Dcel::set_cycles_inout()
{
    inner = 0;
    degen = 0;
    hole = 0;

    for (Cycle *cycle : C) {

        // Hole cycle is a cycle where all hedges have a null left face
        bool is_hole = true;

        Hedge *h = cycle->rep;
        if (!h->left) {
            Hedge *w = h->next;
            while (w != h) {
                if (w->left) {
                    is_hole = false;
                    break;
                }
                w = w->next;
            }
        } else {
            is_hole = false;
        }

        if (is_hole) {
            cycle->inout = Cycle::HOLE;
            hole++;
            continue;
        }

        // Get the leftmost vertex in the cycle
        Vertex *b = h->orig;

        Hedge *e2 = h; // Half-edge starting at leftmost vertex

        Hedge *w = h->next;
        while (w != h) {
            Vertex *p = w->orig;
            E_Int cmp = cmp_vtx(p, b);
            if (cmp < 0) {
                b = p;
                e2 = w;
            }

            w = w->next;
        }

        Hedge *e1 = e2->prev;

        assert(e2->orig == b);
        assert(e1->twin->orig == b);

        Vertex *a = e1->orig;
        Vertex *c = e2->twin->orig;

        // Vectors ab and bc
        E_Float ab[3] = {b->x-a->x, b->y-a->y, b->z-a->z};
        E_Float bc[3] = {c->x-b->x, c->y-b->y, c->z-b->z};

        E_Float N[3];
        K_MATH::cross(ab, bc, N);

        if (Sign(K_MATH::norm(N, 3) == 0)) {
            cycle->inout = Cycle::DEGEN;
            degen++;
        } else {
            cycle->inout = Cycle::INNER;
            inner++;
        }
    }

    printf("Inner cycles: " SF_D_ "\n", inner);
    printf("Degen cycles: " SF_D_ "\n", degen);
    printf("Hole cycles: " SF_D_ "\n", hole);
}

Dcel::~Dcel()
{
    for (Vertex *v : V) delete v;
    for (Hedge *h : H) delete h;
    for (Face *f : F) delete f;
    for (Cycle *c : C) delete c;
}

Dcel::Dcel(const Smesh &Mf, int color)
{
    const auto &X = Mf.X;
    const auto &Y = Mf.Y;
    const auto &Z = Mf.Z;

    V.reserve(Mf.np);
    for (E_Int pid = 0; pid < Mf.np; pid++) {
        Vertex *v = new Vertex(X[pid], Y[pid], Z[pid]);
        v->oids[color] = pid;
        vertex_set.insert(v);
        V.push_back(v);
    }

    init_hedges_and_faces(Mf, color);

    if (check_hedges(H) != 0) {
        fprintf(stderr, "Dcel: Inconsistent half-edge records!\n");
        abort();
    }

    if (check_faces(H, F) != 0) {
        fprintf(stderr, "Dcel: Inconsistent face records!\n");
        abort();
    }

    make_cycles();
    set_cycles_inout();
}



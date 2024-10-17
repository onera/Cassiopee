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

Dcel::Dcel(const Smesh &Mf, E_Int color)
{
    const auto &X = Mf.X;
    const auto &Y = Mf.Y;
    const auto &Z = Mf.Z;

    V.reserve(Mf.np);
    for (E_Int pid = 0; pid < Mf.np; pid++) {
        Vertex *v = new Vertex(X[pid], Y[pid], Z[pid]);
        v->oids[color] = pid;
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
    set_cycles_inout(Mf, color);
}

void Dcel::init_hedges_and_faces(const Smesh &Mf, E_Int color)
{
    printf("Doing color %d\n", color);

    H.reserve(Mf.ne*2);

    std::vector<std::vector<Dcel::Hedge *>> list(Mf.np);

    // Create hedge records

    const auto &E = Mf.E;

    for (E_Int i = 0; i < Mf.ne; i++) {
        const auto &e = E[i];

        E_Int p = e.p;
        E_Int q = e.q;

        Vertex *P = V[p];
        Hedge *h = new Hedge(P, color);
        list[p].push_back(h);
       
        Vertex *Q = V[q];
        Hedge *t = new Hedge(Q, color);
        list[q].push_back(t);

        h->twin = t;
        t->twin = h;

        H.push_back(h);
        H.push_back(t);
    }
    
    // Pair-up hedges

    const auto &pnormals = Mf.pnormals;

    for (E_Int pid = 0; pid < Mf.np; pid++) {
        auto &hedges = list[pid];
        assert(hedges.size() >= 2);

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

    F.reserve(Mf.nf);

    for (E_Int fid = 0; fid < Mf.nf; fid++) {
        const auto &pe = F2E[fid];
        E_Int first_edge = pe[0];
        E_Int where = 2*first_edge;
        Hedge *h = H[where];
        Hedge *t = H[where + 1];
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

        Vertex *P = h->orig;
        Vertex *Q = h->twin->orig;
        if (cmp_vtx(P, Q) < 0) return true;
            
        return false;
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
        if (f->rep->left != f) { assert(0); return 1; }
    }

    puts("CHECK: FACES OK.");

    return 0;
}

void Dcel::make_cycles()
{
    C.clear();

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

void Dcel::set_cycles_inout(const Smesh &Mf, E_Int color)
{
    inner = 0;
    outer = 0;
    degen = 0;
    hole = 0;

    for (Cycle *c : C) {
        Hedge *h = c->rep;
        if (h->left == NULL) {
            c->inout = Cycle::HOLE;
            hole++;
            continue;
        }

        // Get the leftmost vertex in the cycle
        Vertex *leftmost = h->orig;

        Hedge *e2 = h; // Half-edge starting at leftmost vertex
        Hedge *e1 = h->prev; // Half-edge ending at leftmost vertex

        Hedge *w = h->next;
        while (w != h) {
            Vertex *p = w->orig;
            E_Int cmp = cmp_vtx(p, leftmost);
            if (cmp < 0) {
                leftmost = p;
                e2 = w;
                e1 = w->prev;
            }

            w = w->next;
        }

        assert(e2->orig == leftmost);
        assert(e1->twin->orig == leftmost);

        Vertex *a = e1->orig;
        Vertex *b = e2->twin->orig;
        
        E_Float px = leftmost->x - a->x;
        E_Float py = leftmost->y - a->y;
        E_Float pz = leftmost->z - a->z;
        E_Float nx = b->x - leftmost->x;
        E_Float ny = b->y - leftmost->y;
        E_Float nz = b->z - leftmost->z;

        E_Float cp[3] = {py*nz - pz*ny, pz*nx - px*nz, px*ny - py*nx};

        E_Int pid = leftmost->oids[color];

        const E_Float *N = &Mf.pnormals[3*pid];

        E_Float cmp = K_MATH::dot(N, cp, 3);

        if (cmp < 0) {
            c->inout = Cycle::OUTER;
            outer++;
        } else if (cmp == 0) {
            c->inout = Cycle::DEGEN;
            degen++;
        } else {
            c->inout = Cycle::INNER;
            inner++;
        }
    }

    printf("Inner cycles: " SF_D_ "\n", inner);
    printf("Outer cycles: " SF_D_ "\n", outer);
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

/*

void Dcel::init_vertices(const Smesh &M0, const Smesh &M1)
{
    assert(Q.empty());

    for (E_Int i = 0; i < M0.np; i++) {
        //Q.insert(M0.X[i], M0.Y[i], M0.Z[i], M0.l2gp.at(i), Dcel::RED);
        Q.insert(M0.X[i], M0.Y[i], M0.Z[i], i, Dcel::RED);
    }

    for (E_Int i = 0; i < M1.np; i++) {
        //Q.insert(M1.X[i], M1.Y[i], M1.Z[i], M1.l2gp.at(i), Dcel::BLACK);
        Q.insert(M1.X[i], M1.Y[i], M1.Z[i], i, Dcel::BLACK);
    }
}

void Dcel::update_hedge_faces(const std::vector<Face *> &F)
{
    for (Face *f : F) {
        Hedge *h = f->rep;
        h->left = f;
        Hedge *w = h->next;
        while (w != h) {
            w->left = f;
            w = w->next;
        }
    }
}

std::vector<Face *> Dcel::make_cycle_faces(const std::vector<Cycle *> &C)
{
    std::vector<Face *> new_F;

    for (Cycle *c : C) {
        
        // Create a face record
        Face *f = new Face;

        // Set its rep hedge to some edge of the cycle
        Hedge *h = c->rep;
        f->rep = h;

        new_F.push_back(f);
    }

    return new_F;
}

void Dcel::set_face_labels(std::vector<Face *> &F)
{
    // Label each face with the ids of the original faces containing it

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];

        // Get the first RED and BLACK half-edges in the face cycle.
        Hedge *h = f->rep;

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
            // First case: R and B both exist
            assert(R->left);
            assert(B->left);
            assert(R->color == Dcel::RED);
            assert(B->color == Dcel::BLACK);
            f->oid[Dcel::RED] = R->left->oid[Dcel::RED];
            f->oid[Dcel::BLACK] = B->left->oid[Dcel::BLACK];
        } else {
            // Second case: the face is single color
            // Only single color possible is RED, otherwise intersection problem
            // is ill-posed
            Hedge *REP = (R != NULL) ? R : B;
            if (REP != R) {
                hedge_write("black", REP);
            }
            assert(REP == R);
            assert(REP->color == Dcel::RED);
            assert(REP->left);
            f->oid[REP->color] = REP->left->oid[REP->color];
        }
    }
}

Hedge *Dcel::get_hedge_of_color(Face *f, E_Int color)
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

Dcel::Dcel(Smesh &M0, Smesh &M1)
{
    init_vertices(M0, M1);
    Q.inorder(V);
    size_t count = 0;
    for (size_t i = 0; i < V.size(); i++) {
        V[i]->id = i;
        if (V[i]->oid[0] != -1 && V[i]->oid[1] != -1) count++;
    }
    printf("Duplicate points: %lu\n", count);

    init_hedges_and_faces(M0, RED);
    
    init_hedges_and_faces(M1, BLACK);

    assert(check_hedges(H));

    assert(check_faces(H, F));
}

std::vector<Vertex *> Dcel::get_face_vertices(Face *f)
{
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

void Dcel::cut_hedge_at_vertex(Hedge *e, Vertex *x)
{
    // Create two new half-edge records with x as their origin
    Hedge *e1 = new Hedge(x);
    Hedge *e2 = new Hedge(x);

    e1->color = e->color;
    e2->color = e->color;

    Hedge *t = e->twin;

    // Copy the face record
    e1->left = e->left;
    e2->left = t->left;

    // Pair-up the new half-edges
    e->twin = e2;
    e1->twin = t;
    t->twin = e1;
    e2->twin = e;

    // Set prev and next pointers at the endpoints
    e1->next = e->next;
    e2->next = t->next;
    e->next->prev = e1;
    t->next->prev = e2;

    H.push_back(e1);
    H.push_back(e2);

    e->next = e1;
    t->next = e2;

    e1->prev = e;
    e2->prev = t;

    Cp[x].push_back(e1);
    Cp[x].push_back(e2);
}

std::vector<Point> points;

void Dcel::find_intersections_3D(const Smesh &M, const Smesh &S)
{
    puts("    Isolating s_hedges...");

    std::vector<Hedge *> s_hedges;

    for (E_Int i = 2*M.ne; i < 2*(M.ne + S.ne); i += 2) {
        Hedge *h = H[i];
        assert(h->twin == H[i+1]);
        assert(h->color == Dcel::BLACK);
        Hedge *t = h->twin;
        assert(h != NULL);
        assert(t != NULL);
        Vertex *p = h->orig;
        Vertex *q = t->orig;
        assert(p != NULL);
        assert(q != NULL);
        E_Int cmp = cmp_vtx(p, q);
        assert(cmp != 0);
        if (cmp_vtx(p, q) < 0) {
            s_hedges.push_back(h);
        } else {
            s_hedges.push_back(t);
        }
    }

    puts("    Sorting s_hedges...");

    std::sort(s_hedges.begin(), s_hedges.end(), [&] (Hedge *h, Hedge *w)
    {
        assert(h->twin != w);
        assert(w->twin != h);

        E_Int cmp = cmp_vtx(h->orig, w->orig);
        if (cmp < 0) return true;
        if (cmp > 0) return false;

        cmp = cmp_vtx(h->twin->orig, w->twin->orig);
        assert(cmp != 0);

        if (cmp < 0) return true;
        return false;
    });

    puts("    Registering endpoint-edge intersections...");

    init_mh_sv_intersections(M);

    puts("    Tracing edges...");

    size_t before = V.size();
    
    for (size_t hid = 0; hid < s_hedges.size(); hid++) {
        Hedge *sh = s_hedges[hid];

        //printf("Tracing hedge %zu / %zu\n", hid+1, s_hedges.size());

        //trace_hedge(sh, M, S, hid);
        trace_hedge_2(sh, M, S, hid);
    }

    size_t after = V.size();

    printf("Duplicate intersections: %d\n", dup_x);
    printf("Vertices crossed: %lu\n", vertices_crossed.size());

    std::vector<Point> xpoints;
    for (size_t i = before; i < after; i++) {
        Vertex *v = V[i];
        xpoints.push_back(Point(v->x, v->y, v->z));
    }
    point_write("xpoints", xpoints);

    cut_hedges();
}

void Dcel::resolve_hedges(const Smesh &M, const Smesh &S)
{
    assert(Up.empty());
    assert(Lp.empty());

    Up.clear();
    Lp.clear();

    for (Hedge *h : H) {
        
        Hedge *t = h->twin;

        Vertex *p = h->orig;
        Vertex *q = t->orig;

        if (cmp_vtx(p, q) <= 0) {

            Up[p].push_back(h);
            Up[q].push_back(t);

        }
    }

    puts("Resolving vertices...");

    for (size_t i = 0; i < V.size(); i++) {

        //printf("Resolving vertex %d / %zu\n", i+1, V.size());

        Vertex *v = V[i];
        
        std::vector<Hedge *> leaving;

        for (Hedge *h : Up[v]) {
            assert(h->orig == v);
            leaving.push_back(h);
        }

        E_Float N[3] = {0};

        // M point
        if (v->oid[0] != -1) {

            E_Int mpid = v->oid[0];

            const E_Float *pN = &M.pnormals[3*mpid];
            for (E_Int i = 0; i < 3; i++) N[i] = pN[i];

        }
        
        // S point
        else if (v->oid[1] != -1) {

            const auto &loc = v->loc;
            
            E_Int mfid = loc.fid;

            if (loc.e_idx != -1) {

                const auto &pe = M.F2E[mfid];
                E_Int eid = pe[loc.e_idx];
                const auto &pf = M.E2F[eid];
                assert(mfid == pf[0] || mfid == pf[1]);

                E_Int mf1 = pf[0];
                E_Int mf2 = pf[1];

                const E_Float *fN1 = &M.fnormals[3*mf1];
                const E_Float *fN2 = &M.fnormals[3*mf2];

                for (E_Int i = 0; i < 3; i++) {
                    N[i] += fN1[i];
                    N[i] += fN2[i];
                }

                E_Float NORM = K_MATH::norm(N, 3);
                for (E_Int i = 0; i < 3; i++) N[i] /= NORM;

            } else if (loc.v_idx != -1) {
                
                const auto &pn = M.F[mfid];
                E_Int mpid = pn[loc.v_idx];
                const E_Float *pN = &M.pnormals[3*mpid];
                for (E_Int i = 0; i < 3; i++) N[i] = pN[i];

            } else {

                const E_Float *fN = &M.fnormals[3*mfid];

                for (E_Int i = 0; i < 3; i++) N[i] = fN[i];

            }

        }

        // Intersection
        else {

            //Hedge *h = v->xhedge;
            //assert(h);
            //Face *f1 = h->left;
            //Face *f2 = h->twin->left;
            //E_Int mf1 = f1->oid[0];
            //E_Int mf2 = f2->oid[0];

            E_Int eid = v->meid;
            const auto &pf = M.E2F[eid];

            E_Int mf1 = pf[0];
            E_Int mf2 = pf[1];
            assert(mf2 != -1);
    
            const E_Float *fN1 = &M.fnormals[3*mf1];
            const E_Float *fN2 = &M.fnormals[3*mf2];

            for (E_Int i = 0; i < 3; i++) {
                N[i] += fN1[i];
                N[i] += fN2[i];
            }

            E_Float NORM = K_MATH::norm(N, 3);
            for (E_Int i = 0; i < 3; i++) N[i] /= NORM;
        }

        sort_leaving_hedges(leaving, N, M);

        for (size_t i = 0; i < leaving.size(); i++) {
            Hedge *h = leaving[i];
            Hedge *w = leaving[(i+1)%leaving.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }
    }
}

void Dcel::reconstruct(const Smesh &M, const Smesh &S)
{
    check_hedges(H);

    make_cycles();
    assert(0);
    set_cycles_inout(M, 0);

    auto new_F = make_cycle_faces(C);
    
    set_face_labels(new_F);

    update_hedge_faces(new_F);

    for (Face *f : F) delete f;
    F = new_F;

    check_faces(H, F);

    write_degen_faces("degen");
    write_inner_faces("inner");
    write_outer_faces("outer");
}
*/

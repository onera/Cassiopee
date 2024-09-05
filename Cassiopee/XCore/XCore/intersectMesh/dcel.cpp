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
#include "status.h"
#include "segment.h"
#include "primitives.h"
#include "io.h"
#include "hedge.h"
#include "smesh.h"
#include "event.h"
#include "face.h"
#include "cycle.h"
#include "triangle.h"

E_Int Dcel::RED = 0;
E_Int Dcel::BLACK = 1;
E_Int Dcel::NO_IDEA = 2;

void Dcel::write_degen_faces(const char *fname)
{
    auto degen_indices = extract_indices_of_type(Cycle::DEGEN);
    auto degen_faces = extract_faces_of_indices(degen_indices);
    write_ngon(fname, degen_faces);
}

void Dcel::write_outer_faces(const char *fname)
{
    auto outer_indices = extract_indices_of_type(Cycle::OUTER);
    auto outer_faces = extract_faces_of_indices(outer_indices);
    write_ngon(fname, outer_faces);
}

void Dcel::write_inner_faces(const char *fname)
{
    auto inner_indices = extract_indices_of_type(Cycle::INNER);
    auto inner_faces = extract_faces_of_indices(inner_indices);
    write_ngon(fname, inner_faces);
}

void Dcel::write_ngon(const char *fname, const std::vector<Face *> &faces) const
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    E_Int np = 0;
    E_Int ne = 0;
    E_Int nf = (E_Int)faces.size();

    std::map<Vertex *, E_Int> vmap;
    std::vector<Vertex *> new_pids;

    for (Face *f : faces) {
        Hedge *h = f->rep;
        ne++;
        Vertex *p = h->orig;
        if (vmap.find(p) == vmap.end()) {
            vmap[p] = np++;
            new_pids.push_back(p);
        }
        Hedge *w = h->next;
        while (w != h) {
            p = w->orig;
            if (vmap.find(p) == vmap.end()) {
                vmap[p] = np++;
                new_pids.push_back(p);
            }
            ne++;
            w = w->next;
        }
    }

    fprintf(fh, "POINTS\n");
    fprintf(fh, SF_D_ "\n", np);
    for (const auto &v : new_pids) {
        fprintf(fh, "%f %f %f\n", v->x, v->y, v->z);
    }
    
    fprintf(fh, "INDPG\n");
    fprintf(fh, SF_D_ "\n", ne+1);
    E_Int sizeNGon = 0;
    fprintf(fh, SF_D_ " ", sizeNGon);
    for (E_Int i = 0; i < ne; i++) {
        sizeNGon += 2;
        fprintf(fh, SF_D_ " ", sizeNGon);
    }
    assert(sizeNGon == 2*ne);
    fprintf(fh, "\n");

    fprintf(fh, "NGON\n");
    fprintf(fh, SF_D_ "\n", sizeNGon);
    for (Face *f : faces) {
        Hedge *h = f->rep;
        Vertex *p = h->orig;
        Vertex *q = h->twin->orig;
        fprintf(fh, SF_D_ " "  SF_D_ " ", vmap[p], vmap[q]);
        Hedge *w = h->next;
        while (w != h) {
            p = w->orig;
            q = w->twin->orig;
            fprintf(fh, SF_D_ " " SF_D_ " ", vmap[p], vmap[q]);
            w = w->next;
        }
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, SF_D_ "\n", nf+1);
    E_Int sizeNFace = 0;
    fprintf(fh, SF_D_ " ", sizeNFace);
    for (Face *f : faces) {
        Hedge *h = f->rep;
        sizeNFace += 1;
        Hedge *w = h->next;
        while (w != h) {
            assert(w->left == f);
            sizeNFace += 1;
            w = w->next;
        }
        fprintf(fh, SF_D_ " ", sizeNFace);
    }
    fprintf(fh, "\n");

    fprintf(fh, "NFACE\n");
    fprintf(fh, SF_D_ "\n", sizeNFace);
    for (E_Int i = 0; i < sizeNFace; i++)
        fprintf(fh, SF_D_ " ", i);

    fclose(fh);
}

std::vector<Face *> Dcel::extract_faces_of_indices(
    const std::vector<E_Int> &indices)
{
    std::vector<Face *> ret;
    ret.reserve(indices.size());

    for (E_Int index : indices) ret.push_back(F[index]);

    return ret;
}

std::vector<E_Int> Dcel::extract_indices_of_type(E_Int type)
{
    std::vector<E_Int> ret;

    for (size_t i = 0; i < C.size(); i++) {
        if (C[i]->inout == type)
            ret.push_back(i);
    }

    return ret;
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
            Hedge *REP = (R != NULL) ? R : B;
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

Dcel::Dcel(Smesh &M0, Smesh &M1)
{
    init_vertices(M0, M1);
    Q.inorder(V);
    for (size_t i = 0; i < V.size(); i++) {
        V[i]->id = i;
    }

    init_hedges_and_faces(M0, RED);
    
    init_hedges_and_faces(M1, BLACK);

    assert(check_hedges(H));

    assert(check_faces(H, F));
}

void mat3_mult(E_Float A[3][3], E_Float B[3][3], E_Float C[3][3])
{
    for (E_Int i = 0; i < 3; i++) {
        for (E_Int j = 0; j < 3; j++) {
            C[i][j] = 0;

            for (E_Int k = 0; k < 3; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void mat3_vec(E_Float A[3][3], E_Float x[3], E_Float b[3])
{
    for (E_Int i = 0; i < 3; i++) {
        b[i] = 0;
        for (E_Int j = 0; j < 3; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
}

void Dcel::init_hedges_and_faces(Smesh &M, E_Int color)
{
    printf("Doing color %d\n", color);
    size_t nh = H.size();
    size_t nhh = nh + 2 * M.E.size();

    H.reserve(nhh);

    std::vector<std::vector<Hedge *>> list(M.np);

    for (E_Int i = 0; i < M.ne; i++) {
        const auto &e = M.E[i];

        E_Int p = e.p;
        E_Int q = e.q;

        Event *xit = Q.lookup(M.X[p], M.Y[p], M.Z[p]);
        assert(xit);

        Vertex *P = xit->key;

        Hedge *h = new Hedge(P);
        h->eid = i;

        list[p].push_back(h);

        xit = Q.lookup(M.X[q], M.Y[q], M.Z[q]);
        assert(xit);        

        Vertex *V = xit->key;
        
        Hedge *t = new Hedge(V);
        t->eid = i;

        list[q].push_back(t);

        h->twin = t;
        t->twin = h;

        h->color = color;
        t->color = color;

        H.push_back(h);
        H.push_back(t);
    }
    
    // Pair-up hedges

    const auto &pnormals = M.pnormals;

    for (E_Int pid = 0; pid < M.np; pid++) {
        auto &hedges = list[pid];

        assert(!hedges.empty());

        const E_Float *N = &pnormals[3*pid];
        assert(Sign(K_MATH::norm(N, 3)-1) == 0);

        sort_leaving_hedges(hedges, N, M);

        for (size_t i = 0; i < hedges.size(); i++) {
            Hedge *h = hedges[i];
            Hedge *w = hedges[(i+1)%hedges.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }

        Event *xit = Q.lookup(hedges[0]->orig);

        xit->key->rep = hedges[0];
    }



    for (E_Int i = 0; i < M.nf; i++) {
        const auto &edges = M.F2E[i];
        E_Int first_edge = edges[0];
        E_Int where = nh + 2 * first_edge;
        Hedge *h = H[where];
        Hedge *t = H[where + 1];
        assert(h->twin == t);
        assert(t->twin == h);

        Face *f = new Face;
        f->oid[color] = i;

        assert(M.E2F[first_edge][0] == (E_Int)i || M.E2F[first_edge][1] == E_Int(i));
        Hedge *REP = (M.E2F[first_edge][0] == (E_Int)i) ? h : t;

        assert(REP->left == NULL);

        f->rep = REP;
        REP->left = f;
        Hedge *w = REP->next;
        while (w != REP) { w->left = f; w = w->next; }
        
        F.push_back(f);
    }

    for (Face *f : F) {
        Hedge *h = f->rep;
        assert(h->left == f);
        Hedge *w = h->next;
        while (w != h) {
            assert(w->left == f);
            w = w->next;
        }
    }
    
    // Create the unbounded faces
    f_unbounded[color] = new Face;
    f_unbounded[color]->oid[color] = -1;

    // Set it as the left face for hedges without a left face
    for (size_t i = nh; i < nhh; i++) {
        if (H[i]->left == NULL)
            H[i]->left = f_unbounded[color];
    }
}

E_Int Dcel::check_hedges(const std::vector<Hedge *> &H)
{
    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];
        if (h->prev->next != h) { assert(0); return 0; }
        if (h->next->prev != h) { assert(0); return 0; }
        if (h->twin->twin != h) { assert(0); return 0; }
        if (h->twin->next->orig != h->orig) { assert(0); return 0; }
        if (h->prev->twin->orig != h->orig) { assert(0); return 0; }
    }

    puts("CHECK: EDGES OK.");

    return 1;
}


E_Int Dcel::check_faces(const std::vector<Hedge *> &H,
    const std::vector<Face *> &F)
{
    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];
        if (h->prev->left != h->left) { assert(0); return 0; }
        if (h->next->left != h->left) { assert(0); return 0; }
    }

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];
        if (f->rep->left != f) { assert(0); return 0; }
    }

    puts("CHECK: FACES OK.");

    return 1;
}

Dcel::~Dcel()
{
    delete f_unbounded[0];
    delete f_unbounded[1];

    for (size_t i = 0; i < V.size(); i++) delete V[i];
    for (size_t i = 0; i < H.size(); i++) delete H[i];
    for (size_t i = 0; i < F.size(); i++) delete F[i];
    for (size_t i = 0; i < C.size(); i++) delete C[i];

    Q.drop();
}

void Dcel::set_cycles_inout(const Smesh &M, const Smesh &S)
{
    E_Int inner = 0;
    E_Int outer = 0;
    E_Int degen = 0;

    for (Cycle *c : C) {
        // Get the leftmost vertex in the cycle
        Hedge *h = c->rep;
        Vertex *v = h->orig;

        Hedge *e2 = h; // Half-edge starting at v
        Hedge *e1 = h->prev; // Half-edge ending at v

        Hedge *w = h->next;
        while (w != h) {
            Vertex *p = w->orig;
            E_Int cmp = cmp_vtx(p, v);
            if (cmp < 0) {
                v = p;
                e2 = w;
                e1 = w->prev;
            }

            w = w->next;
        }

        assert(e2->orig == v);
        assert(e1->twin->orig == v);

        c->left = v;

        Vertex *a = e1->orig;
        Vertex *b = e2->twin->orig;

        // If the angle from e1 to e2 is less than 180°, c is an outer cycle.
        // Else, c is an inner cycle.
        
        E_Float px = v->x - a->x;
        E_Float py = v->y - a->y;
        E_Float pz = v->z - a->z;
        E_Float nx = b->x - v->x;
        E_Float ny = b->y - v->y;
        E_Float nz = b->z - v->z;

        E_Float cp[3] = {py*nz - pz*ny, pz*nx - px*nz, px*ny - py*nx};


        E_Float N[3]= { };

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

            Hedge *h = v->xhedge;
            assert(h);

            Face *f1 = h->left;
            Face *f2 = h->twin->left;

            E_Int mf1 = f1->oid[0];
            E_Int mf2 = f2->oid[0];
    
            const E_Float *fN1 = &M.fnormals[3*mf1];
            const E_Float *fN2 = &M.fnormals[3*mf2];

            for (E_Int i = 0; i < 3; i++) {
                N[i] += fN1[i];
                N[i] += fN2[i];
            }

            E_Float NORM = K_MATH::norm(N, 3);
            for (E_Int i = 0; i < 3; i++) N[i] /= NORM;
        }

        E_Float NORM = K_MATH::norm(N, 3);
        assert(Sign(NORM -1) == 0);

        E_Float cmp = K_MATH::dot(N, cp, 3);

        if (cmp < 0) {
            c->inout = Cycle::INNER;
            inner++;
        } else if (cmp == 0) {
            c->inout = Cycle::DEGEN;
            degen++;
        } else {
            c->inout = Cycle::OUTER;
            outer++;
        }
    }

    printf("Inner cycles: " SF_D_ "\n", inner);
    printf("Outer cycles: " SF_D_ "\n", outer);
    printf("Degen cycles: " SF_D_ "\n", degen);
    
    printf("Total faces: " SF_D_ "\n", outer);
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

void Dcel::locate_spoints(const Smesh &M, const Smesh &S)
{
    for (E_Int sp = 0; sp < S.np; sp++) {

        Event *xit = Q.lookup(S.X[sp], S.Y[sp], S.Z[sp]);
        assert(xit);

        Vertex *V = xit->key;
        auto &ploc = V->loc;

        E_Int found = 0;

        E_Int voxel_x = floor((S.X[sp] - M.xmin) / M.HX);
        E_Int voxel_y = floor((S.Y[sp] - M.ymin) / M.HY);
        E_Int voxel_z = floor((S.Z[sp] - M.zmin) / M.HZ);
        E_Int sp_bin = voxel_x + M.NX * voxel_y + M.NXY * voxel_z;

        auto it = M.fmap.find(sp_bin);

        assert(it != M.fmap.end());

        const auto &pf = it->second;

        for (size_t mf = 0; mf < pf.size() && !found; mf++) {

            const auto &pn = M.F[pf[mf]];

            E_Float o[3] = {0, 0, 0};

            for (E_Int p : pn) {
                o[0] += M.X[p];
                o[1] += M.Y[p];
                o[2] += M.Z[p];
            }
            for (E_Int i = 0; i < 3; i++) o[i] /= pn.size(); 
    
            for (size_t i = 0; i < pn.size(); i++) {
                E_Int p = pn[i];
                E_Int q = pn[(i+1)%pn.size()];

                E_Float u, v, w;

                if (Triangle::is_point_inside(
                    S.X[sp], S.Y[sp], S.Z[sp],
                    M.X[p], M.Y[p], M.Z[p],
                    M.X[q], M.Y[q], M.Z[q],
                    o[0], o[1], o[2],
                    u, v, w)) {

                    found = 1;

                    ploc.fid = pf[mf];

                    if (Sign(v) == 0) ploc.e_idx = i;
                    else if (Sign(1-u) == 0) ploc.v_idx = (i+1)%pn.size();
                    else if (Sign(1-w) == 0) ploc.v_idx = i;

                    break;
                }
            }
        }

        assert(found);
    }
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

void Dcel::handle_intersecting_endpoint(Vertex *v, const Smesh &M)
{
    if (!Cp[v].empty()) return;

    const auto &vloc = v->loc;

    if (vloc.e_idx == -1) return;

    E_Int fid = vloc.fid;

    const auto &pe = M.F2E[fid]; 

    E_Int me = pe[vloc.e_idx];

    Hedge *start = H[2*me];

    Face *face = F[fid];

    if (start->left != face) start = start->twin;

    Hedge *h = start;

    E_Int done = 0;

    while (1) {

        if (hedge_contains_vertex(h, v)) {
            done = 1;
            v->xhedge = h;
            cut_hedge_at_vertex(h, v);
            break;
        }

        h = h->next;
        if (h == start) break;
    }

    assert(done == 1);
}

E_Int Dcel::get_next_face(const Smesh &M, E_Float px, E_Float py, E_Float pz,
    const std::vector<E_Int> &pf, E_Float dir[3])
{
    E_Int next_face = -1;
    E_Float t_min = EFLOATMAX;

    for (size_t i = 0; i < pf.size(); i++) {

        E_Int fid = pf[i];

        Face *face = F[fid];

        const E_Float *fN = &M.fnormals[3*fid];

        E_Float proj[3] = { };
        E_Float dp = K_MATH::dot(fN, dir, 3);
        for (E_Int j = 0; j < 3; j++) proj[j] = dir[j] - dp * fN[j];

        Hedge *h = face->rep;

        E_Int hit = 0;

        while (1) {

            Vertex *a = h->orig;
            Vertex *b = h->twin->orig;

            E_Float dx = px + 10000 * proj[0];
            E_Float dy = py + 10000 * proj[1];
            E_Float dz = pz + 10000 * proj[2];

            E_Float t;

            hit = EdgeEdgeIntersect(
                px, py, pz,
                dx, dy, dz,
                a->x, a->y, a->z,
                b->x, b->y, b->z,
                t);

            if (hit) {
                if (t < t_min) {
                    next_face = fid;
                    t_min = t;
                }

                break;
            }
    
            h = h->next;
            if (h == face->rep) break;
        }
    }

    /*
    if (next_face == -1) {

        point_write("test_point", px, py, pz);

        for (E_Int fid : pf) {
            char fname[128] = {};
            sprintf(fname, "test_face_%d", fid);
            face_write(fname, F[fid]);
        }
    }
    */

    return next_face;
}

void Dcel::trace_hedge(Hedge *sh, const Smesh &M, const Smesh &S, E_Int hid)
{
    Vertex *p = sh->orig;
    Vertex *q = sh->twin->orig;

    E_Float dir[3] = {q->x-p->x, q->y-p->y, q->z-p->z};

    const auto &ploc = p->loc; 

    std::vector<E_Int> test_faces;

    E_Int mfid = ploc.fid;

    // Get the potential starting faces

    if (ploc.e_idx != -1) {
        assert(ploc.v_idx == -1);
        const auto &pe = M.F2E[mfid];
        E_Int eid = pe[ploc.e_idx];
        const auto &pf = M.E2F[eid];
        for (E_Int fid : pf) test_faces.push_back(fid);
    } else if (ploc.v_idx != -1) {
        assert(ploc.e_idx == -1);
        const auto &pn = M.F[mfid];
        E_Int pid = pn[ploc.v_idx];
        const auto &pf = M.P2F[pid];
        for (E_Int fid : pf) test_faces.push_back(fid);
    } else {
        test_faces.push_back(mfid);
    }

    // Handle potential intersection of starting point

    handle_intersecting_endpoint(p, M);
    handle_intersecting_endpoint(q, M);

    // Determine the starting face
    E_Int start_face = get_next_face(M, p->x, p->y, p->z, test_faces, dir);

    assert(start_face != -1);

    // Trace
    
    E_Int found = 0;
    E_Int walk = 0;
    E_Int max_walk = 10;

    Face *current_face = F[start_face];

    E_Float px = p->x, py = p->y, pz = p->z;

    Hedge *start_hedge = current_face->rep;

    Hedge *current_hedge = sh;

    // Pinpoint the endpoint
    std::vector<E_Int> end_faces;
    const auto &qloc = q->loc;
    E_Int qfid = qloc.fid;

    if (qloc.e_idx != -1) {
        assert(qloc.v_idx == -1);
        const auto &pe = M.F2E[qfid];
        E_Int eid = pe[qloc.e_idx];
        const auto &pf = M.E2F[eid];
        for (E_Int fid : pf) end_faces.push_back(fid);
    } else if (qloc.v_idx != -1) {
        assert(qloc.e_idx == -1);
        const auto &pn = M.F[qfid];
        E_Int pid = pn[qloc.v_idx];
        const auto &pf = M.P2F[pid];
        for (E_Int fid : pf) end_faces.push_back(fid);
    } else {
        end_faces.push_back(qfid);
    }


    while (!found && walk < max_walk) {

        // Check if we reached q

        for (E_Int fid : end_faces) {
            if (F[fid] == current_face) {
                found = 1;
                break;
            }
        }

        if (found) break;

        E_Int current_fid = current_face->oid[0];

        const E_Float *fN = &M.fnormals[3*current_fid];

        E_Float proj[3] = { };
        E_Float dp = K_MATH::dot(fN, dir, 3);
        for (E_Int i = 0; i < 3; i++) proj[i] = dir[i] - dp * fN[i];

        E_Float dx = px + 2*proj[0];
        E_Float dy = py + 2*proj[1];
        E_Float dz = pz + 2*proj[2];

        Hedge *h = current_face->rep;
        E_Int reached = 0;
        E_Int hit = 0;

        E_Float ix, iy, iz;
        ix = iy = iz = -10000;
    
        while (!reached && !found) {
    
            Vertex *a = h->orig;
            Vertex *b = h->twin->orig;

            hit = EdgeEdgeIntersect(
                px, py, pz,
                dx, dy, dz,
                a->x, a->y, a->z,
                b->x, b->y, b->z,
                ix, iy, iz);

            if (hit) {

                Vertex *x = NULL;

                E_Int hit_a = cmp_points(ix, iy, iz, a->x, a->y, a->z) == 0;
                E_Int hit_b = cmp_points(ix, iy, iz, b->x, b->y, b->z) == 0;

                // Hit a vertex: original m vertex, or intersection

                if (hit_a) x = a;
                else if (hit_b) x = b;

                if (x != NULL) {

                    // Stop if reached destination
                    if (x->oid[1] != -1) {
                        assert(x == q);
                        found = 1;
                    }

                    // M point, get the next face
                    else if (x->oid[0] != -1) {

                        E_Int mpid = x->oid[0];
                        const auto &pf = M.P2F[mpid];
                        E_Int next_fid = get_next_face(M, x->x, x->y, x->z, pf, dir);
                        assert(next_fid != -1);
                        assert(next_fid != current_fid);
                        current_face = F[next_fid];

                    } else {
                        
                        // E_Intersection, move

                        current_face = h->twin->left;
    
                    }
                } else {

                    // Hit the inside of an edge

                    // Must be a new intersection
                    Event *xit = Q.lookup(ix, iy, iz);

                    assert(xit == NULL);

                    x = new Vertex(ix, iy, iz);
                    x->id = V.size();
                    V.push_back(x);
                    x->xhedge = h;

                    cut_hedge_at_vertex(h, x);

                    current_face = h->twin->left;

                }

                if (found) break;
                
                assert(x);

                cut_hedge_at_vertex(current_hedge, x);
                current_hedge = current_hedge->next;

                px = ix;
                py = iy;
                pz = iz;

                break;
            }

            h = h->next;
            if (h == start_hedge) {
                reached = 1;
            }
        }

        assert(reached == 0);
        walk++;
    }

    assert(walk < max_walk);
}

void Dcel::find_intersections_3D(const Smesh &M, const Smesh &S)
{
    puts("Isolating s_hedges...");

    std::vector<Hedge *> s_hedges;

    for (E_Int i = 2*M.ne; i < 2*(M.ne + S.ne); i += 2) {
        Hedge *h = H[i];
        assert(h->twin == H[i+1]);
        assert(h->color == Dcel::BLACK);
        Hedge *t = h->twin;
        Vertex *p = h->orig;
        Vertex *q = t->orig;
        if (cmp_vtx(p, q) <= 0) {
            s_hedges.push_back(h);
        } else {
            s_hedges.push_back(t);
        }
    }

    puts("Sorting s_hedges...");

    std::sort(s_hedges.begin(), s_hedges.end(), [&] (Hedge *h, Hedge *w)
    {
        return cmp_vtx(h->orig, w->orig) <= 0;
    });

    puts("Tracing edges...");

    for (size_t hid = 0; hid < s_hedges.size(); hid++) {
        Hedge *sh = s_hedges[hid];

        //printf("Tracing hedge %d / %zu\n", hid+1, s_hedges.size());

        trace_hedge(sh, M, S, hid);
    }
}

void Dcel::sort_leaving_hedges(std::vector<Hedge *> &leaving,
    const E_Float N[3],
    const Smesh &M) const
{
    // Choose a vector that is not parallel to N

    E_Float ref_vec[3] = {0, N[2], -N[1]};
    
    if (Sign(K_MATH::norm(ref_vec, 3)) == 0) {
        ref_vec[0] = -N[2];
        ref_vec[1] = 0;
        ref_vec[2] =  N[0];
        assert(Sign(K_MATH::norm(ref_vec, 3)) != 0);
    }

    E_Float dp = K_MATH::dot(ref_vec, N, 3);

    for (E_Int i = 0; i < 3; i++) ref_vec[i] = ref_vec[i] - dp * N[i];

    E_Float NORM = K_MATH::norm(ref_vec, 3);

    for (E_Int i = 0; i < 3; i++) ref_vec[i] /= NORM;

    assert(Sign(K_MATH::norm(ref_vec, 3) - 1) == 0);

    std::vector<E_Float> angles;

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

        for (E_Int j = 0; j < 3; j++) {
            PQ_proj[j] = PQ[j] - dp * N[j];
        }

        E_Float costheta = K_MATH::dot(ref_vec, PQ_proj, 3) / K_MATH::norm(PQ_proj, 3);
        
        costheta = std::min(costheta, 1.0);
        
        costheta = std::max(costheta, -1.0);
        
        assert(costheta >= -1 && costheta <= 1);
        
        E_Float angle = acos(costheta);
        
        // Determine the direction of the angle
        E_Float C[3] = {};
        
        K_MATH::cross(ref_vec, PQ_proj, C);
        
        if (K_MATH::dot(N, C, 3) > 0)
            angle = 2*K_MATH::PI - angle;
        
        //angle = angle * 180 / K_MATH::PI;
        
        angles.push_back(angle);
    }

    std::vector<E_Int> indices(leaving.size());
    for (size_t i = 0; i < leaving.size(); i++)
        indices[i] = i;

    std::sort(indices.begin(), indices.end(), [&](E_Int i, E_Int j)
    {
        if (angles[i] < angles[j]) return true;
        
        else if (angles[i] > angles[j]) return false;
        
        else {
            Hedge *h = leaving[i];
            Hedge *w = leaving[j];
            
            assert(h->color != w->color);

            Vertex *P = h->orig;
            Vertex *Q = h->twin->orig;
            if (cmp_vtx(P, Q) < 0) return true;
            
            return false;
        }
    });
    
    std::vector<Hedge *> tmp(leaving);
    
    for (size_t i = 0; i < leaving.size(); i++)
        leaving[i] = tmp[indices[i]];
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

        E_Int do_sort = 0;

        for (size_t i = 1; i < leaving.size(); i++) {
            if (leaving[i]->color != leaving[0]->color) {
                do_sort = 1;
                break;
            }
        }

        if (!do_sort) continue;
        

        E_Float N[3]= { };

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

            Hedge *h = v->xhedge;
            assert(h);

            Face *f1 = h->left;
            Face *f2 = h->twin->left;

            E_Int mf1 = f1->oid[0];
            E_Int mf2 = f2->oid[0];
    
            const E_Float *fN1 = &M.fnormals[3*mf1];
            const E_Float *fN2 = &M.fnormals[3*mf2];

            for (E_Int i = 0; i < 3; i++) {
                N[i] += fN1[i];
                N[i] += fN2[i];
            }

            E_Float NORM = K_MATH::norm(N, 3);
            for (E_Int i = 0; i < 3; i++) N[i] /= NORM;
        }

        E_Float NORM = K_MATH::norm(N, 3);
        assert(Sign(NORM -1) == 0);

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

    set_cycles_inout(M, S);

    auto new_F = make_cycle_faces(C);
    
    set_face_labels(new_F);

    update_hedge_faces(new_F);

    for (Face *f : F) delete f;
    F = new_F;

    check_faces(H, F);

    //write_degen_faces("degen");
    //write_inner_faces("inner");
}

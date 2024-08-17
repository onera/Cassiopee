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
#include "sweep.h"
#include "io.h"
#include "hedge.h"
#include "smesh.h"
#include "event.h"
#include "face.h"
#include "cycle.h"
#include "triangle.h"

Int Dcel::RED = 0;
Int Dcel::BLACK = 1;
Int Dcel::NO_IDEA = 2;

void Dcel::resolve(Vertex *p, const std::vector<Segment *> &L,
    const std::vector<Segment *> &C, const std::vector<Segment *> &U,
    std::vector<Hedge *> &H)
{
    // The half-edges starting at p
    std::vector<Hedge *> leaving;

    for (Segment *s : U) {
        Hedge *h = s->rep;
        assert(h->orig == p);
        leaving.push_back(h);
    }

    for (Segment *s : L) {
        Hedge *h = s->rep;
        Hedge *t = h->twin;
        assert(t->orig == p);
        leaving.push_back(t);
    }

    for (Segment *s : C) {
        // Create two new half-edges records with p as their origin
        Hedge *e1 = new Hedge(p);
        Hedge *e2 = new Hedge(p);

        e1->color = s->color;
        e2->color = s->color;

        // Half-edge connected to segment s
        Hedge *e = s->rep;
        Hedge *t = e->twin;

        // Copy the face data
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

        leaving.push_back(e1);
        leaving.push_back(e2);

        H.push_back(e1);
        H.push_back(e2);

        // Segment s is represented by the new half-edge extending to the right
        // of the sweep line
        s->rep = e1;
    }

    // Correct the situation around p
    // Order leaving half-edges in clockwise order around p

    Hedge::sort_cwise(leaving, 0, leaving.size() - 1);

    // Link-up consecutive half-edges from different sets
    for (size_t i = 0; i < leaving.size(); i++) {
        Hedge *h = leaving[i];
        Hedge *w = leaving[(i+1)%leaving.size()];

        h->twin->next = w;
        w->prev = h->twin;
    }
}

void Dcel::find_intersections()
{
    Status T;

    Float BIG = 1;

    for (Vertex *p : V) {
        while (fabs(p->x) >= BIG || fabs(p->y) >= BIG) {
            BIG *= 2;
        }
    }

    T.rx = T.ry = -BIG;
    
    std::vector<Segment *> S;

    for (size_t i = 0; i < H.size(); i += 2) {
        Hedge *h = H[i];
        assert(h->twin = H[i+1]);

        Segment *s = new Segment(h, i >> 1);
        
        S.push_back(s);
    }

    Segment::sort(S, 0, S.size()-1, cmp_mySeg);
    //std::reverse(S.begin(), S.end());

    Vertex *lowerLeft = new Vertex(-BIG, -BIG, 0);
    Vertex *lowerRight = new Vertex(BIG, -BIG, 0);
    Vertex *upperLeft = new Vertex(-BIG, BIG, 0);
    Vertex *upperRight = new Vertex(BIG, BIG, 0);

    Segment *lowerSentinel = new Segment(lowerLeft, lowerRight, S.size());
    Segment *upperSentinel = new Segment(upperLeft, upperRight, S.size()+1);

    T.insert(lowerSentinel);
    T.insert(upperSentinel);

    sweep(Q, T, S, V, H);

    check_hedges(H);

    make_cycles();

    set_cycles_inout();

    auto new_F = make_cycle_faces(C);
    
    set_face_labels(new_F);

    update_hedge_faces(new_F);

    for (Face *f : F) delete f;
    F = new_F;

    check_faces(H, F);

    // Clean up
    T.erase(lowerSentinel);
    T.erase(upperSentinel);

    delete lowerLeft;
    delete lowerRight;
    delete upperLeft;
    delete upperRight;
    delete upperSentinel;
    delete lowerSentinel;
    for (Segment *s : S) delete s;
}

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

    Int np = 0;
    Int ne = 0;
    Int nf = (Int)faces.size();

    std::map<Vertex *, Int> vmap;
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
    Int sizeNGon = 0;
    fprintf(fh, SF_D_ " ", sizeNGon);
    for (Int i = 0; i < ne; i++) {
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
    Int sizeNFace = 0;
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
    for (Int i = 0; i < sizeNFace; i++)
        fprintf(fh, SF_D_ " ", i);

    fclose(fh);
}

std::vector<Face *> Dcel::extract_faces_of_indices(
    const std::vector<Int> &indices)
{
    std::vector<Face *> ret;
    ret.reserve(indices.size());

    for (Int index : indices) ret.push_back(F[index]);

    return ret;
}

std::vector<Int> Dcel::extract_indices_of_type(Int type)
{
    std::vector<Int> ret;

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
        Int RB = 0;

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

Hedge *Dcel::get_hedge_of_color(Face *f, Int color)
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

    for (Int i = 0; i < M0.np; i++) {
        //Q.insert(M0.X[i], M0.Y[i], M0.Z[i], M0.l2gp.at(i), Dcel::RED);
        Q.insert(M0.X[i], M0.Y[i], M0.Z[i], i, Dcel::RED);
    }

    for (Int i = 0; i < M1.np; i++) {
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

void mat3_mult(Float A[3][3], Float B[3][3], Float C[3][3])
{
    for (Int i = 0; i < 3; i++) {
        for (Int j = 0; j < 3; j++) {
            C[i][j] = 0;

            for (Int k = 0; k < 3; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void mat3_vec(Float A[3][3], Float x[3], Float b[3])
{
    for (Int i = 0; i < 3; i++) {
        b[i] = 0;
        for (Int j = 0; j < 3; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
}

void Dcel::init_hedges_and_faces(Smesh &M, Int color)
{
    //printf("Doing color %d\n", color);
    size_t nh = H.size();
    size_t nhh = nh + 2 * M.E.size();

    H.reserve(nhh);

    std::vector<std::vector<Hedge *>> list(M.np);

    for (Int i = 0; i < M.ne; i++) {
        const auto &e = M.E[i];

        Int p = e.p;
        Int q = e.q;

        Event *xit = Q.lookup(M.X[p], M.Y[p], M.Z[p]);
        assert(xit);

        Vertex *P = xit->key;

        Hedge *h = new Hedge(P);

        list[p].push_back(h);

        xit = Q.lookup(M.X[q], M.Y[q], M.Z[q]);
        assert(xit);        

        Vertex *V = xit->key;
        
        Hedge *t = new Hedge(V);

        list[q].push_back(t);

        h->twin = t;
        t->twin = h;

        h->color = color;
        t->color = color;

        H.push_back(h);
        H.push_back(t);
    }
    
    // Sort hedges

    const auto &pnormals = M.pnormals;

    Float ez[3] = {0, 0, 1};
    Float I[3][3] = {}; I[0][0] = I[1][1] = I[2][2] = 1;

    for (Int pid = 0; pid < M.np; pid++) {
        auto &hedges = list[pid];
        const Float *N = &pnormals[3*pid];

        // Project the hedges tails onto the plane (i, N)

        for (size_t i = 0; i < hedges.size(); i++) {
            Hedge *h = hedges[i];
            Hedge *t = h->twin;

            Vertex *tail = t->orig;

            Float dp = tail->x*N[0] + tail->y*N[1] + tail->z*N[2];

            h->proj_tx = tail->x - dp * N[0]; 
            h->proj_ty = tail->y - dp * N[1]; 
            h->proj_tz = tail->z - dp * N[2]; 
        }

        // Rotation axis
        Float r[3];
        K_MATH::cross(N, ez, r);
        Float normr = K_MATH::norm(r, 3);
        if (Sign(normr == 0)) normr = 1;
        for (int i = 0; i < 3; i++) r[i] /= normr;

        // Rotation angle
        Float costheta = K_MATH::dot(N, ez, 3) / (K_MATH::norm(N, 3) * K_MATH::norm(ez, 3));
        assert(costheta > -1-TOL && costheta < 1+TOL);

        Float sintheta = sqrt(1 - costheta*costheta); 

        // Create the rotation matrix
        Float K[3][3];
        K[0][0] = 0;     K[0][1] = -r[2]; K[0][2] = r[1];
        K[1][0] = r[2];  K[1][1] = 0;     K[1][2] = -r[0];
        K[2][0] = -r[1]; K[2][1] = r[0];  K[2][2] = 0;

        Float K2[3][3];
        mat3_mult(K, K, K2);

        // Find the rotation matrix
        Float R[3][3] = {};
        
        for (Int i = 0; i < 3; i++) {
            for (Int j = 0; j < 3; j++) {
                R[i][j] = I[i][j] + sintheta * K[i][j] + (1-costheta)*K2[i][j];
                //printf("%.3f ", R[i][j]);
            }
        }
    
        // Project the origin
        Float o[3] = {M.X[pid], M.Y[pid], M.Z[pid]};
        Float proj_o[3] = {};
        mat3_vec(R, o, proj_o);

        for (size_t i = 0; i < hedges.size(); i++) {
            Hedge *h = hedges[i];

            h->proj_ox = proj_o[0];
            h->proj_oy = proj_o[1];
            h->proj_oz = proj_o[2];

            Float t[3] = {h->proj_tx, h->proj_ty, h->proj_tz};
            Float proj_t[3];
            mat3_vec(R, t, proj_t);

            h->proj_tx = proj_t[0];
            h->proj_ty = proj_t[1];
            h->proj_tz = proj_t[2];
        }
            
        Hedge::sort_cwise(hedges, 0, hedges.size()-1);

        for (size_t i = 0; i < hedges.size(); i++) {
            Hedge *h = hedges[i];
            Hedge *w = hedges[(i+1)%hedges.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }

        assert(!hedges.empty());

        Event *xit = Q.lookup(hedges[0]->orig);

        xit->key->rep = hedges[0];
    }

    for (Int i = 0; i < M.nf; i++) {
        const auto &edges = M.F2E[i];
        Int first_edge = edges[0];
        Int where = nh + 2 * first_edge;
        Hedge *h = H[where];
        Hedge *t = H[where + 1];
        assert(h->twin == t);
        assert(t->twin == h);

        Face *f = new Face;
        f->oid[color] = i;
        //f->oid[color] = M.l2gf.at(i);


        assert(M.E2F[first_edge][0] == (Int)i || M.E2F[first_edge][1] == Int(i));
        Hedge *REP = (M.E2F[first_edge][0] == (Int)i) ? h : t;

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

Int Dcel::check_hedges(const std::vector<Hedge *> &H)
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


Int Dcel::check_faces(const std::vector<Hedge *> &H,
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
}

void Dcel::set_cycles_inout()
{
    Int inner = 0;
    Int outer = 0;
    Int degen = 0;

    for (Cycle *c : C) {
        // Get the leftmost vertex in the cycle
        Hedge *h = c->rep;
        Vertex *v = h->orig;

        Hedge *e2 = h; // Half-edge starting at v
        Hedge *e1 = h->prev; // Half-edge ending at v

        Hedge *w = h->next;
        while (w != h) {
            Vertex *p = w->orig;
            Int cmp = cmp_vtx(p, v);
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
        
        Float px = v->x - a->x;
        Float py = v->y - a->y;
        Float pz = v->z - a->z;
        Float nx = b->x - v->x;
        Float ny = b->y - v->y;
        Float nz = b->z - v->z;

        Float cp[3] = {py*nz - pz*ny, pz*nx - px*nz, px*ny - py*nx};

        // TODO(Imad): compute r more rigorously...
        Float r[3] = {0, 0, 1};
        Float cmp = Sign(K_MATH::dot(r, cp, 3)); 

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

    //printf("Inner cycles: " SF_D_ "\n", inner);
    //printf("Outer cycles: " SF_D_ "\n", outer);
    //printf("Degen cycles: " SF_D_ "\n", degen);
    
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
    for (Int sp = 0; sp < S.np; sp++) {

        Event *xit = Q.lookup(S.X[sp], S.Y[sp], S.Z[sp]);
        assert(xit);

        Vertex *V = xit->key;
        auto &ploc = V->loc;

        Int found = 0;

        for (Int mf = 0; mf < M.nf && !found; mf++) {

            const auto &pn = M.F[mf];

            Float o[3] = {0, 0, 0};

            for (Int p : pn) {
                o[0] += M.X[p];
                o[1] += M.Y[p];
                o[2] += M.Z[p];
            }
            for (Int i = 0; i < 3; i++) o[i] /= pn.size(); 
    
            for (size_t i = 0; i < pn.size(); i++) {
                Int p = pn[i];
                Int q = pn[(i+1)%pn.size()];

                Float u, v, w;

                if (Triangle::is_point_inside(
                    S.X[sp], S.Y[sp], S.Z[sp],
                    M.X[p], M.Y[p], M.Z[p],
                    M.X[q], M.Y[q], M.Z[q],
                    o[0], o[1], o[2],
                    u, v, w)) {

                    found = 1;

                    ploc.fid = mf;

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

    Int fid = vloc.fid;

    const auto &pe = M.F2E[fid]; 

    Int me = pe[vloc.e_idx];

    Hedge *start = H[2*me];

    Face *face = F[fid];

    if (start->left != face) start = start->twin;

    Hedge *h = start;

    Int done = 0;

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

Int Dcel::get_next_face(const Smesh &M, Float px, Float py, Float pz,
    const std::vector<Int> &pf, Float dir[3])
{
    Int next_face = -1;

    for (size_t i = 0; i < pf.size(); i++) {

        Int fid = pf[i];

        Face *face = F[fid];

        const Float *fN = &M.fnormals[3*fid];

        Float proj[3] = { };
        Float dp = K_MATH::dot(fN, dir, 3);
        for (Int j = 0; j < 3; j++) proj[j] = dir[j] - dp * fN[j];

        Hedge *h = face->rep;

        Int hit = 0;

        while (1) {

            Vertex *a = h->orig;
            Vertex *b = h->twin->orig;

            Float dx = px + proj[0];
            Float dy = py + proj[1];
            Float dz = pz + proj[2];

            hit = EdgeEdgeIntersect(
                px, py, pz,
                dx, dy, dz,
                a->x, a->y, a->z,
                b->x, b->y, b->z);

            if (hit) break;
    
            h = h->next;
            if (h == face->rep) break;
        }

        if (hit) {
            next_face = fid;
            break;
        }
    }

    assert(next_face != -1);

    return next_face;
}

void Dcel::trace_hedge(Hedge *sh, const Smesh &M, const Smesh &S)
{
    Vertex *p = sh->orig;
    Vertex *q = sh->twin->orig;

    Float dir[3] = {q->x-p->x, q->y-p->y, q->z-p->z};

    const auto &ploc = p->loc; 

    std::vector<Int> test_faces;

    Int mfid = ploc.fid;

    // Get the potential starting faces

    if (ploc.e_idx != -1) {
        assert(ploc.v_idx == -1);
        const auto &pe = M.F2E[mfid];
        Int eid = pe[ploc.e_idx];
        const auto &pf = M.E2F[eid];
        for (Int fid : pf) test_faces.push_back(fid);
    } else if (ploc.v_idx != -1) {
        assert(ploc.e_idx == -1);
        const auto &pn = M.F[mfid];
        Int pid = pn[ploc.v_idx];
        const auto &pf = M.P2F[pid];
        for (Int fid : pf) test_faces.push_back(fid);
    } else {
        test_faces.push_back(mfid);
    }

    // Handle potential intersection of starting point

    handle_intersecting_endpoint(p, M);
    handle_intersecting_endpoint(q, M);

    // Determine the starting face
    Int start_face = get_next_face(M, p->x, p->y, p->z, test_faces, dir);

    // Trace
    
    Int found = 0;
    Int walk = 0;
    Int max_walk = 10;

    Face *current_face = F[start_face];

    Float px = p->x, py = p->y, pz = p->z;

    Hedge *start_hedge = current_face->rep;

    Hedge *current_hedge = sh;

    while (!found && walk < max_walk) {

        Int qfid = q->loc.fid;
        
        if (F[qfid] == current_face) {
            found = 1;
            break;
        }

        Int current_fid = current_face->oid[0];

        const Float *fN = &M.fnormals[3*current_fid];

        Float proj[3] = { };
        Float dp = K_MATH::dot(fN, dir, 3);
        for (Int i = 0; i < 3; i++) proj[i] = dir[i] - dp * fN[i];

        Float dx = px + proj[0];
        Float dy = py + proj[1];
        Float dz = pz + proj[2];

        Hedge *h = current_face->rep;
        Int reached = 0;
        Int hit = 0;

        Float ix, iy, iz;
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

                Int hit_a = cmp_points(ix, iy, iz, a->x, a->y, a->z) == 0;
                Int hit_b = cmp_points(ix, iy, iz, b->x, b->y, b->z) == 0;

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

                        Int mpid = x->oid[0];
                        const auto &pf = M.P2F[mpid];
                        Int next_fid = get_next_face(M, x->x, x->y, x->z, pf, dir);
                        assert(next_fid != current_fid);
                        current_face = F[next_fid];

                    } else {
                        
                        // Intersection, move

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
}

void Dcel::find_intersections_3D(const Smesh &M, const Smesh &S)
{
    std::vector<Hedge *> s_hedges;

    for (Int i = 2*M.ne; i < 2*(M.ne + S.ne); i += 2) {
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

    std::sort(s_hedges.begin(), s_hedges.end(), [&] (Hedge *h, Hedge *w)
    {
        return cmp_vtx(h->orig, w->orig) <= 0;
    });

    for (size_t hid = 0; hid < s_hedges.size(); hid++) {
        Hedge *sh = s_hedges[hid];

        trace_hedge(sh, M, S);
    }
}

void Dcel::resolve_hedges(const Smesh &M, const Smesh &S)
{
    E_Float ez[3] = {0, 0, 1};
    E_Float I[3][3] = {};
    I[0][0] = I[1][1] = I[2][2] = 1;

    assert(Up.empty());
    assert(Lp.empty());

    Up.clear();
    Lp.clear();
    //Cp.clear();

    for (Hedge *h : H) {
        
        Hedge *t = h->twin;

        Vertex *p = h->orig;
        Vertex *q = t->orig;

        if (cmp_vtx(p, q) <= 0) {

            Up[p].push_back(h);
            Up[q].push_back(t);

        }
    }

    for (size_t i = 0; i < V.size(); i++) {

        Vertex *v = V[i];
        
        std::vector<Hedge *> leaving;

        for (Hedge *h : Up[v]) {
            assert(h->orig == v);
            leaving.push_back(h);
        }

        Int do_sort = 0;

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

            Int mpid = v->oid[0];

            const Float *pN = &M.pnormals[3*mpid];
            for (Int i = 0; i < 3; i++) N[i] = pN[i];

        }
        
        // S point
        else if (v->oid[1] != -1) {

            const auto &loc = v->loc;
            
            Int mfid = loc.fid;

            if (loc.e_idx != -1) {

                const auto &pe = M.F2E[mfid];
                Int eid = pe[loc.e_idx];
                const auto &pf = M.E2F[eid];
                assert(mfid == pf[0] || mfid == pf[1]);

                Int mf1 = pf[0];
                Int mf2 = pf[1];

                const Float *fN1 = &M.fnormals[3*mf1];
                const Float *fN2 = &M.fnormals[3*mf2];

                for (Int i = 0; i < 3; i++) {
                    N[i] += fN1[i];
                    N[i] += fN2[i];
                }

                Float NORM = K_MATH::norm(N, 3);
                for (Int i = 0; i < 3; i++) N[i] /= NORM;

            } else if (loc.v_idx != -1) {
                
                const auto &pn = M.F[mfid];
                Int mpid = pn[loc.v_idx];
                const Float *pN = &M.pnormals[3*mpid];
                for (Int i = 0; i < 3; i++) N[i] = pN[i];

            } else {

                const Float *fN = &M.fnormals[3*mfid];

                for (Int i = 0; i < 3; i++) N[i] = fN[i];

            }

        }

        // Intersection
        else {

            Hedge *h = v->xhedge;
            assert(h);

            Face *f1 = h->left;
            Face *f2 = h->twin->left;

            Int mf1 = f1->oid[0];
            Int mf2 = f2->oid[0];
    
            const Float *fN1 = &M.fnormals[3*mf1];
            const Float *fN2 = &M.fnormals[3*mf2];

            for (Int i = 0; i < 3; i++) {
                N[i] += fN1[i];
                N[i] += fN2[i];
            }

            Float NORM = K_MATH::norm(N, 3);
            for (Int i = 0; i < 3; i++) N[i] /= NORM;
        }

        Float NORM = K_MATH::norm(N, 3);
        assert(Sign(NORM -1) == 0);
        
        for (Hedge *h : leaving) {

            // Project h on the plane with normal pnormal
            
            Hedge *w = h->twin;
            assert(w->twin == h);

            Vertex *tail = w->orig;

            Float dp = tail->x*N[0] + tail->y*N[1] + tail->z*N[2];

            h->proj_tx = tail->x - dp * N[0]; 
            h->proj_ty = tail->y - dp * N[1]; 
            h->proj_tz = tail->z - dp * N[2]; 

            // Rotation axis
            Float r[3];
            K_MATH::cross(N, ez, r);
            Float normr = K_MATH::norm(r, 3);
            if (Sign(normr == 0)) normr = 1;
            for (int i = 0; i < 3; i++) r[i] /= normr;

            // Rotation angle
            Float costheta = K_MATH::dot(N, ez, 3) / (K_MATH::norm(N, 3) * K_MATH::norm(ez, 3));
            assert(costheta > -1-TOL && costheta < 1+TOL);

            Float sintheta = sqrt(1 - costheta*costheta); 

            // Create the rotation matrix
            Float K[3][3];
            K[0][0] = 0;     K[0][1] = -r[2]; K[0][2] = r[1];
            K[1][0] = r[2];  K[1][1] = 0;     K[1][2] = -r[0];
            K[2][0] = -r[1]; K[2][1] = r[0];  K[2][2] = 0;

            Float K2[3][3];
            mat3_mult(K, K, K2);

            // Find the rotation matrix
            Float R[3][3] = {};
        
            for (Int i = 0; i < 3; i++) {
                for (Int j = 0; j < 3; j++) {
                    R[i][j] = I[i][j] + sintheta * K[i][j] + (1-costheta)*K2[i][j];
                }
            }
     
            // Project the origin
            Vertex *orig = h->orig;
            assert(orig == v);
            Float o[3] = {orig->x, orig->y, orig->z};
            Float proj_o[3] = {};
            mat3_vec(R, o, proj_o);

            h->proj_ox = proj_o[0];
            h->proj_oy = proj_o[1];
            h->proj_oz = proj_o[2];

            Float t[3] = {h->proj_tx, h->proj_ty, h->proj_tz};
            Float proj_t[3];
            mat3_vec(R, t, proj_t);

            h->proj_tx = proj_t[0];
            h->proj_ty = proj_t[1];
            h->proj_tz = proj_t[2];
        }

        Hedge::sort_cwise(leaving, 0, leaving.size()-1);

        for (size_t i = 0; i < leaving.size(); i++) {
            Hedge *h = leaving[i];
            Hedge *w = leaving[(i+1)%leaving.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }

    }
}

void Dcel::reconstruct()
{
    check_hedges(H);

    make_cycles();

    set_cycles_inout();

    auto new_F = make_cycle_faces(C);
    
    set_face_labels(new_F);

    update_hedge_faces(new_F);

    for (Face *f : F) delete f;
    F = new_F;

    check_faces(H, F);
}

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

E_Int Dcel::RED = 0;
E_Int Dcel::BLACK = 1;
E_Int Dcel::NO_IDEA = 2;

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

    E_Float BIG = 1;

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

    Vertex *lowerLeft = new Vertex(-BIG, -BIG);
    Vertex *lowerRight = new Vertex(BIG, -BIG);
    Vertex *upperLeft = new Vertex(-BIG, BIG);
    Vertex *upperRight = new Vertex(BIG, BIG);

    Segment *lowerSentinel = new Segment(lowerLeft, lowerRight, S.size());
    Segment *upperSentinel = new Segment(upperLeft, upperRight, S.size()+1);

    T.insert(lowerSentinel);
    T.insert(upperSentinel);

    sweep(Q, T, S, V, H);

    printf("Intersections: %lu\n", V.size());

    point_write("ipoints", V);

    check_hedges(H);

    make_cycles();

    set_cycles_inout();

    auto new_F = make_cycle_faces(C);
    
    set_face_labels(new_F);

    update_hedge_faces(new_F);

    for (Face *f : F) delete f;
    F = new_F;

    check_faces(H, F);

    /*
    write_degen_faces("degen");
    write_outer_faces("outer");
    write_inner_faces("inner");
    */


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
        fprintf(fh, "%f %f 0\n", v->x, v->y);
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
            f->oid[Dcel::RED] = R->left->oid[Dcel::RED];
            f->oid[Dcel::BLACK] = B->left->oid[Dcel::BLACK];
        } else {
            // Second case: the face is single color
            Hedge *REP = (R != NULL) ? R : B;
            assert(REP->left);
            //assert(REP->left->oid[REP->color] != -1);
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
        Q.insert(M0.X[i], M0.Y[i], M0.l2gp.at(i), Dcel::RED);
    }

    for (E_Int i = 0; i < M1.np; i++) {
        Q.insert(M1.X[i], M1.Y[i], M1.l2gp.at(i), Dcel::BLACK);
    }
}

Dcel::Dcel(const Smesh &M0, const Smesh &M1)
{
    init_vertices(M0, M1);
    Q.inorder(V);
    for (size_t i = 0; i < V.size(); i++) {
        Vertex *v = V[i];
        if (v->oid[0] != -1 && v->oid[1] != -1) {
            printf("Vertex %lu is dup (" SF_D_ " " SF_D_ ")\n", i,
                v->oid[0], v->oid[1]);
        }
        V[i]->id = i;
    }

    init_hedges_and_faces(M0, RED);

    init_hedges_and_faces(M1, BLACK);

    assert(check_hedges(H));

    assert(check_faces(H, F));
}

void Dcel::init_hedges_and_faces(const Smesh &M, E_Int color)
{
    size_t nh = H.size();
    size_t nhh = nh + 2 * M.E.size();

    H.reserve(nhh);

    std::vector<std::vector<Hedge *>> list(M.np);

    for (E_Int i = 0; i < M.ne; i++) {
        const auto &e = M.E[i];

        E_Int p = e.p;
        E_Int q = e.q;

        Event *xit = Q.lookup(M.X[p], M.Y[p]);
        assert(xit);

        Hedge *h = new Hedge(xit->key);

        list[p].push_back(h);

        xit = Q.lookup(M.X[q], M.Y[q]);
        assert(xit);        
        Hedge *t = new Hedge(xit->key);

        list[q].push_back(t);

        h->twin = t;
        t->twin = h;

        h->color = color;
        t->color = color;

        H.push_back(h);
        H.push_back(t);
    }

    for (size_t i = 0; i < list.size(); i++) {
        auto &hedges = list[i];

        Hedge::sort_cwise(hedges, 0, hedges.size()-1);

        for (size_t j = 0; j < hedges.size(); j++) {
            Hedge *h = hedges[j];
            Hedge *w = hedges[(j + 1) % hedges.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }

        assert(!hedges.empty());

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
        f->oid[color] = M.l2gf.at(i);

        assert(M.E2F[first_edge][0] == (E_Int)i || M.E2F[first_edge][1] == E_Int(i));

        Hedge *REP = (M.E2F[first_edge][0] == (E_Int)i) ? h : t;

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

    puts("DCEL HEDGES OK.");

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

    puts("DCEL FACES OK.");

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
            E_Int cmp = compare(*p, *v);
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
        E_Float nx = b->x - v->x;
        E_Float ny = b->y - v->y;
        E_Int cmp = Sign(px * ny - py * nx);

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

#include "proto.h"
#include <cstdio>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <unordered_map>

void dcel_set_face_labels(std::vector<face *> &F)
{
    // Label each face with the ids of the original faces containing it.

    for (size_t i = 0; i < F.size(); i++) {
        face *f = F[i];

        // Get the first RED and BLACK half-edges in the face cycle.
        hedge *h = f->rep;
        assert(h);

        hedge *R = NULL;
        hedge *B = NULL;
        int RB = 0;

        if (h->color == RED) {
            R = h;
            B = dcel_get_hedge_of_color_from_face(f, BLACK);
            if (B) RB = 1;
        } else if (h->color == BLACK) {
            B = h;
            R = dcel_get_hedge_of_color_from_face(f, RED);
            if (R) RB = 1;
        } else {
            assert(0);
        }

        // First case: R and B both exist
        if (RB) {
            if (R->left) f->ofp[RED] = R->left->id;
            if (B->left) f->ofp[BLACK] = B->left->id;
        } else {
            
            // The face cycle is single color.
            hedge *rep = (R != NULL) ? R : B;
            assert(rep);
            int COLOR = rep->color;

            if (rep->left) f->ofp[COLOR] = rep->left->id;

            vertex *v = dcel_get_leftmost_vertex_of_outer_cycle(f);
            hedge *h = v->left;

            if (h == NULL)
                continue;

            cycle *c = h->cycl;
            //cycle *cc = c->prev;
            //cycle *cc = dcel_get_outer_cycle_from_cycle(c);
            face *fp = dcel_get_face_from_cycle(c);

            if (fp == NULL) {
                // fp is unlimited, nothing to do.
                continue;
            }

            int OTHER = abs(COLOR-1);

            hedge *hp = dcel_get_hedge_of_color_from_face(fp, OTHER);

            if (hp) {
                if (hp->left)
                    f->ofp[OTHER] = hp->left->id;
            }
        }
    }
}

void face::print_parents()
{
    printf("<");
    printf("%d,", ofp[RED]);
    printf("%d", ofp[BLACK]);
    puts(">");
}

void dcel_make_faces_from_connected_components(const std::vector<cycle *> &comps,
    std::vector<face *> &new_faces)
{ 
    assert(new_faces.empty()); 

    for (cycle *comp : comps) {
        // Create a face record
        face *f = new face;

        // Set its rep hedge to some edge of comp
        hedge *h = comp->rep;
        f->rep = h;

        // Set the left face record of all hedges to f
        //h->left = f;
        hedge *w = h->next;
        while (w != h) {
            //w->left = f;
            w = w->next;
        }

        // Construct the list of holes that exist within the face
        cycle *hole = comp->next;
        while (hole != NULL) {
            //assert("UNIMPLEMENTED holes in faces" && 0);
            hedge *h = hole->rep;
            f->inner.push_back(h);

            //h->left = f;
            hedge *w = h->next;
            while (w != h) {
                //w->left = f;
                w = w->next;
            }

            hole = hole->next;
        }

        new_faces.push_back(f);
    }
}

void dcel_make_connected_components(const std::vector<cycle *> &C,
    std::vector<cycle *> &connected_components)
{
    for (cycle *c : C) {
        if (c->inout == DEGEN)
            continue;
        
        if (c->inout == OUTER) {
            connected_components.push_back(c);
            continue;
        }

        // Leftmost vertex in the cycle
        vertex *p = c->left;

        // Half-edge immediately to the left of p
        hedge *h = p->left;

        // Cycle corresponding to h
        if (h != NULL) {    
            cycle *cc = h->cycl;
            assert(cc);

            // Link cycles c and cc
            c->prev = cc;
            cc->next = c;
        }
    }
}

void dcel::_init_from_smesh(const smesh &M, int color)
{
    const auto& X = M.X;
    const auto& Y = M.Y;
    const auto& E = M.E;

    for (size_t i = 0; i < X.size(); i++) {
        Q.insert(X[i], Y[i], i);
    }

    std::unordered_map<vertex *, std::vector<hedge *>> v2h;

    size_t nh = H.size();

    for (size_t i = 0; i < E.size(); i++) {
        vertex *p = Q.locate_v(X[E[i].p], Y[E[i].p]);
        vertex *q = Q.locate_v(X[E[i].q], Y[E[i].q]);
        assert(p && q);

        hedge *h = new hedge(p);
        hedge *t = new hedge(q);
        h->twin = t;
        t->twin = h;
        h->color = t->color = color;

        H.push_back(h);
        H.push_back(t);

        v2h[p].push_back(h);
        v2h[q].push_back(t);
    }

    for (auto &it : v2h) {
        auto &v = it.first;
        auto &hedges = it.second;

        hedge_sort_cwise(hedges, 0, hedges.size()-1);

        for (size_t j = 0; j < hedges.size(); j++) {
            hedge *e = hedges[j];
            hedge *w = hedges[(j+1)%hedges.size()]; 
            e->twin->next = w;
            w->prev = e->twin;
        }

        if (hedges.size())
            v->inc = hedges[0];
    }

    const auto &E2F = M.E2F;
    const auto &F2E = M.F2E;

    for (size_t i = 0; i < M.F.size(); i++) {
        const auto &edgs = F2E[i];
        int first_edge = edgs[0];
        hedge *h = H[nh+2*first_edge];
        hedge *t = H[nh+2*first_edge+1];
        assert(h->twin == t);
        assert(t->twin == h);

        face *f = new face;
        f->id = i;
        
        assert(E2F[first_edge][0] == (int)i || E2F[first_edge][1] == (int)i);

        hedge *rep = E2F[first_edge][0] == (int)i ? h : t;

        f->rep = rep;
        rep->left = f;
        hedge *w = rep->next;
        while (w != rep) {
            w->left = f;
            w = w->next;
        }

        F.push_back(f);
    }
}

dcel::dcel(const smesh &M0, const smesh &M1)
{
    _init_from_smesh(M0, RED);
    _init_from_smesh(M1, BLACK);

    assert(dcel_check_hedges_without_faces(H));

    assert(dcel_check_faces(H, F));

    // Init vertex ids
    Q.inorder(V);

    for (size_t i = 0; i < V.size(); i++)
        V[i]->nid = i;
}

int dcel_check_faces(const std::vector<hedge *> &H,
    const std::vector<face *> &F)
{
    printf("\nDCEL faces: %zu\n", F.size());
    
    for (hedge *h : H) {
        if (h->prev->left != h->left) { assert(0); return 0; };
        if (h->next->left != h->left) { assert(0); return 0; };
    }

    for (face *f : F) {
        if (f->rep->left != f) { assert(0); return 0; };
    }

    puts("DCEL FACES OK.\n");

    return 1;
}

const char *color_to_str(int color)
{
    if (color == RED) return "R";
    if (color == BLACK) return "B";
    return "NO_IDEA";
}

std::vector<hedge *> dcel_get_incident_edges(vertex *v)
{
    std::vector<hedge *> list;
    
    hedge *h = v->inc;
    if (h == NULL)
        return list;


    list.push_back(h);
    hedge *w = h->twin->next;
    while (w != h) {
        list.push_back(w);
        w = w->twin->next;
    }
    return list;
}

void dcel_resolve(vertex *v, std::vector<segment *> &L,
    std::vector<segment *> &C, std::vector<segment *> &U,
    std::vector<hedge *> &H)
{
    // The half-edges starting at v
    std::vector<hedge *> leaving;
    
    for (segment *s : U) {
        hedge *h = s->rep;
        assert(h->orig == v);
        leaving.push_back(h);
    }

    for (segment *s : L) {
        hedge *h = s->rep;
        hedge *t = h->twin;
        assert(t->orig == v);
        leaving.push_back(t);
    } 

    for (segment *s : C) {
        // Create two new half-edge records with v as their origin
        hedge *e1 = new hedge(v);
        hedge *e2 = new hedge(v);

        e1->color = s->color();
        e2->color = s->color();

        // Half-edge connected to segment s
        hedge *e = s->rep;
        hedge *t = e->twin;

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
        e1->s = e->s;
        e2->s = t->s;
    }
    
    // Correct the situation around v
    // Order leaving half-edges in clockwise order around v

    hedge_sort_cwise(leaving, 0, leaving.size()-1); 

    // Link-up consecutive half-edges from different sets
    for (size_t i = 0; i < leaving.size(); i++) {
        hedge *e = leaving[i];
        hedge *w = leaving[(i+1)%leaving.size()];

        e->twin->next = w;
        w->prev = e->twin;
    }
}

int dcel_check_hedges_without_faces(const std::vector<hedge *> &H)
{
    printf("\nDCEL half-edges: %zu\n", H.size());

    for (size_t i = 0; i < H.size(); i++) {
        hedge *h = H[i];
        if (h->prev->next != h) { assert(0); return 0; };
        if (h->next->prev != h) { assert(0); return 0; };
        if (h->twin->twin != h) { assert(0); return 0; };
        if (h->twin->next->orig != h->orig) { assert(0); return 0; };
        if (h->prev->twin->orig != h->orig) { assert(0); return 0; };
    }

    puts("DCEL HALF-EDGES OK.\n");
    return 1;
}

cycle::cycle(hedge *h)
: rep(h), inout(NO_IDEA), left(NULL), prev(NULL), next(NULL)
{}

void cycle::print()
{
    if (inout == INNER) printf("INNER: ");
    else if (inout == OUTER) printf("OUTER: ");
    else if (inout == DEGEN) printf("DEGEN: ");
    else printf("NO IDEA: ");
    hedge *h = rep;
    vertex *v = h->orig;
    v->print();
    hedge *w = h->next;
    while(w != h) {
        v = w->orig;
        v->print();
        w = w->next;
    }
    puts("\n");
}

std::vector<cycle *> dcel_make_cycles(const std::vector<hedge *> &H)
{
    std::vector<cycle *> C;

    // Overwrite existing cycles!
    for (hedge *h : H)
        h->cycl = NULL;

    for (size_t i = 0; i < H.size(); i++) {
        hedge *h = H[i];

        if (h->cycl)
            continue;

        cycle *c = new cycle(h);
        C.push_back(c);

        h->cycl = c;

        hedge *w = h->next;
        while (w != h) {
            w->cycl = c;
            w = w->next;
        }
    }

    return C;
}

void dcel_set_cycles_inout(std::vector<cycle *> &C)
{
    // Distinguish outer and inner cycles
    for (cycle *c : C) {
        // Get the leftmost vertex in the cycle
        hedge *h = c->rep;
        vertex *v = h->orig;

        hedge *e2 = h;       // half-edge starting at v
        hedge *e1 = h->prev; // half-edge ending at v

        hedge *w = h->next;
        while (w != h) {
            vertex *p = w->orig;
            int cmp = vertex_cmp(p, v);
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

        vertex *a = e1->orig;
        vertex *b = e2->twin->orig; 

        // If the angle from e1 to e2 is less that 180°, c is an outer cycle.
        // Else, c is an inner cycle.
        double px = v->x - a->x;
        double py = v->y - a->y;
        double nx = b->x - v->x;
        double ny = b->y - v->y;
        int cmp = sign(px*ny - py*nx);

        // Note(Imad): the angle can be zero if the cycle is degenerate.
        // In this case we consider the cycle to be an outer cycle.
        if (cmp < 0) c->inout = INNER;
        else if (cmp == 0) c->inout = DEGEN;
        else c->inout = OUTER;
    }
}

dcel::dcel()
: V(), H(), F(), C(), Q()
{}

void dcel::find_intersections()
{
    double BIG = 1;
    for (hedge *h : H) {
        vertex *p = h->orig;
        vertex *q = h->twin->orig;

        while (fabs(p->x) >= BIG || fabs(p->y) >= BIG || fabs(q->x) >= BIG ||
            fabs(q->y) >= BIG)
            BIG *= 2;
    }

    printf("-BIG: %f\n", -BIG);

    vertex *vinf0 = new vertex(-BIG, -BIG, -1);
    vertex *vinf1 = new vertex( BIG, -BIG, -2);
    vertex *vinf2 = new vertex(-BIG,  BIG, -3);
    vertex *vinf3 = new vertex( BIG,  BIG, -4);
    segment *sinf0 = new segment(vinf0, vinf1, -1);
    segment *sinf1 = new segment(vinf2, vinf3, -2);

    status T;
    T.xs = T.ys = -BIG;
    T.insert(sinf0);
    T.insert(sinf1);

    std::vector<segment *> S;
    queue Q;

    for (size_t i = 0; i < H.size(); i += 2) {
        hedge *h = H[i];
        hedge *t = H[i+1];
        assert(h->twin = t);
        assert(t->twin = h);
        
        // Does the swap if necessary
        segment *s = new segment(h, i>>1);

        event *it = Q.insert(s->p);
        s->p = it->key;

        S.push_back(s);
    }

    segment_sort(S, segment_cmp_lexico, 0, S.size()-1);
    
    for (size_t i = 0; i < S.size()-1; i++)
        assert(segment_cmp_lexico(S[i], S[i+1]) < 0);

    std::reverse(S.begin(), S.end());

    sweep(S, V, H, Q, T);

    assert(dcel_check_hedges_without_faces(H));
    
    C = dcel_make_cycles(H);
    dcel_set_cycles_inout(C);

    //for (cycle *c : C) c->print();

    std::vector<cycle *> connected_components;
    dcel_make_connected_components(C, connected_components);
    printf("Connected components: %zu\n", connected_components.size());

    // Make the new faces
    std::vector<face *> new_faces;
    dcel_make_faces_from_connected_components(connected_components, new_faces);
    printf("New faces: %zu\n", new_faces.size());

    dcel_set_face_labels(new_faces);

    F = new_faces;
}

void dcel::set_color(int color)
{
    assert(color == RED || color == BLACK);
    for (hedge *h : H)
        h->color = color;
}

void dcel::write_edges(const char *fname)
{
    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "Points\n");
    fprintf(fh, "%zu\n", V.size());
    for (size_t i = 0; i < V.size(); i++)
        fprintf(fh, "%f %f 0\n", V[i]->x, V[i]->y);
    fprintf(fh, "Edges\n");
    fprintf(fh, "%zu\n", H.size());
    for (size_t i = 0; i < H.size(); i++) {
        hedge *h = H[i];
        //if (H[i]->color == DEGEN)
        //    continue;
        hedge *t = h->twin;
        h->orig->fprint(fh);
        t->orig->fprint(fh);
        //fprintf(fh, "%d %d ", h->orig->nid, t->orig->nid);
    }
    fprintf(fh, "\n");

    fclose(fh);
}

void dcel::write_ngon(const char *fname)
{
    assert(0 && "UNIMPLEMENTED.");
    FILE *fh = fopen(fname, "w");
    assert(fh);

    fprintf(fh, "Points\n");
    fprintf(fh, "%zu\n", V.size());
    for (size_t i = 0; i < V.size(); i++)
        fprintf(fh, "%f %f 0\n", V[i]->x, V[i]->y);
   
    fprintf(fh, "INDPG\n");
    fprintf(fh, "%zu\n", H.size()+1);
    fprintf(fh, "0 ");
    int indpg = 0;
    for (size_t i = 0; i < H.size(); i++) {
        indpg += 2;
        fprintf(fh, "%d ", indpg);
    }
    fprintf(fh, "\n");
 
    fprintf(fh, "NGON\n");
    fprintf(fh, "%zu\n", H.size());
    for (size_t i = 0; i < H.size(); i++) {
        hedge *h = H[i];
        hedge *t = h->twin;
        fprintf(fh, "%d %d ", h->orig->nid, t->orig->nid);
    }
    fprintf(fh, "\n");

    fprintf(fh, "INDPH\n");
    fprintf(fh, "%zu\n", C.size()+1);
    fprintf(fh, "0 ");
    int indph = 0;
    for (size_t i = 0; i < C.size(); i++) {
        indph += C[i]->get_size(); 
        fprintf(fh, "%d ", indph);
    }
    fprintf(fh, "\n"); 
    
    fprintf(fh, "NFACE\n");
    fprintf(fh, "%d\n", indph);
    for (size_t i = 0; i < C.size(); i++) {
        hedge *h = C[i]->rep;
        vertex *v = h->orig;
        fprintf(fh, "%d ", v->nid);
        hedge *w = h->next;
        while (w != h) {
           v = w->orig;
           fprintf(fh, "%d ", v->nid);
           w = w->next;
        }
    }
    fprintf(fh, "\n"); 

    fclose(fh);
}

int cycle::get_size()
{
    hedge *h = rep;
    assert(h);
    int ret = 1;
    hedge *w = h->next;
    while (w != h) {
        ret++;
        w = w->next;
    }
    return ret;
}

face::face()
{
    rep = NULL;
    inner = std::vector<hedge *>();
    id = -1;
    ofp[0] = -1;
    ofp[1] = -1;
}

vertex *dcel_get_leftmost_vertex_of_outer_cycle(face *f)
{
    hedge *h = f->rep;
    vertex *v = h->orig;
    hedge *w = h->next;
    while (w != h) {
        vertex *p = w->orig;
        int cmp = vertex_cmp(p, v);
        if (cmp < 0)
            v = p;
        w = w->next;
    }
    return v;
}

face *dcel_get_face_from_cycle(cycle *c)
{
    if (!c) return NULL;
    return c->rep->left;
}

cycle *dcel_get_outer_cycle_from_cycle(cycle *c)
{
    if (!c) return NULL;
    cycle *outer = c;
    while (outer->prev) {
        outer = outer->prev;
    }
    assert(outer);
    return outer;
}

int dcel_faces_are_adjacent(face *f0, face *f1)
{
    assert(f0);
    assert(f1);

    auto H0 = dcel_get_face_outer_loop(f0);
    auto H1 = dcel_get_face_outer_loop(f1);

    for (size_t i = 0; i < H0.size(); i++) {
        hedge *h0 = H0[i];
        hedge *t0 = h0->twin;

        for (size_t j = 0; j < H1.size(); j++) {
            hedge *h1 = H1[j];
            if (t0 == h1)
                return 1;
        }
    }

    return 0;
}

std::vector<hedge *> dcel_get_face_outer_loop(face *f)
{
    std::vector<hedge *> ret;
    hedge *h = f->rep;
    ret.push_back(h);
    hedge *w = h->next;
    while (w != h) {
        ret.push_back(w);
        w = w->next;
    }
    return ret;
}

hedge *dcel_get_hedge_of_color_from_loop(hedge *h, int color)
{
    if (h->color == color)
        return h;

    hedge *w = h->next;
    while (w != h) {
        if (w->color == color)
            return w;
        w = w->next;
    }

    return NULL;
}

hedge *dcel_get_hedge_of_color_from_face(face *f, int color)
{
    hedge *ret = NULL;

    // Outer loop
    hedge *h = f->rep;
    ret = dcel_get_hedge_of_color_from_loop(h, color);
    if (ret)
        return ret;

    // Inner loop
    for (hedge *h : f->inner) {
        ret = dcel_get_hedge_of_color_from_loop(h, color);
        if (ret)
            return ret;
    }

    return ret;
}

void face::print_vertices()
{
    hedge *h = rep;
    printf("%d ", h->orig->nid);
    hedge *w = h->next;
    while (w != h) {
        printf("%d ", w->orig->nid);
        w = w->next;
    }
    puts("");
}

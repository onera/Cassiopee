#include "proto.h"
#include <cstddef>
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cmath>

static
void _find_new_event(queue &Q, status &T, std::vector<vertex *> &V,
    snode *sit0, snode *sit1)
{
    segment *s0 = sit0->s;
    segment *s1 = sit1->s;

    vertex a = *(s0->p);
    vertex b = *(s0->q);
    vertex c = *(s1->p);
    vertex d = *(s1->q);
    
    double denom = a.x * (d.y - c.y) +
                   b.x * (c.y - d.y) + 
                   d.x * (b.y - a.y) +
                   c.x * (a.y - b.y);

    if (sign(denom) == 0) return;

    double num = a.x * (d.y - c.y) +
                 c.x * (a.y - d.y) +
                 d.x * (c.y - a.y);

    double s = num / denom;

    if (s < -TOL || s > 1.0+TOL) return;

    num = - ( a.x * (c.y - b.y) +
              b.x * (a.y - c.y) +
              c.x * (b.y - a.y) );

    double t = num / denom;

    if (t < -TOL || t > 1.0+TOL) return;

    double x = a.x + s * (b.x - a.x);
    double y = a.y + s * (b.y - a.y);

    // Where does the intersection lie wrt to sweep vertex
    int cmp = xy_cmp(T.xs, T.ys, x, y);

    if (cmp == 1) return;

    // Intersection is to the right of sweep vertex

    // Does a vertex with the same coordinate already exist in Q?
    event *I = Q.locate(x, y);

    if (I == NULL) {
        // No: insert new event and set its sit to S1
        vertex *v = new vertex(x, y, V.size(), true);
        I = Q.insert(v);
        V.push_back(v);
    }

    // Set intersection info of S1 to I
    I->inf = s0;
    sit0->inf = I->key;
}


static
int _report_intersection(const std::vector<segment *> &L,
    const std::vector<segment *> &C, const std::vector<segment *> &U)
{
    int ncolors[2] = {0, 0};

    for (segment *s : L) ncolors[s->color()]++;
    for (segment *s : U) ncolors[s->color()]++;
    for (segment *s : C) ncolors[s->color()]++;

    // We have an intersection if both colors exist
    if (ncolors[RED] == 0 || ncolors[BLACK] == 0)
        return 0;

    /*
    printf("Intersection at P%d (%f %f)\n", v->id, v->x, v->y);
    printf("Ending: ");
    for (segment *s : L) printf("S%d ", s->id);
    puts("");
    printf("Passing: ");
    for (segment *s : C) printf("S%d ", s->id);
    puts("");
    printf("Starting: ");
    for (segment *s : U) printf("S%d ", s->id);
    puts("\n");
    */
    return 1;
}

static
void _handle_event(event *E, queue &Q, status &T, std::vector<vertex *> &V,
    std::vector<segment *> &S, std::vector<hedge *> &H)
{
    // Get starting segments
    segment *s = E->inf;
    vertex *v = E->key;
    std::vector<segment *> U;

    while (S.size() && S.back()->p == v) {
        U.push_back(S.back());
        S.pop_back();
    }

    snode *sit = NULL;
    
    if (s == NULL) {
        segment dummy(v, v, -1);
        sit = T.locate(&dummy);
        if (sit)
            s = sit->s;
    }

    // Handle passing and ending segments
    std::vector<segment *> L, C;
    std::vector<void *> C_info;
    
    segment *s_succ = NULL;
    segment *s_pred = NULL;

    if (s != NULL) {
        sit = T.lookup(s);
        assert(sit);

        while (sit->inf == v || sit->inf == T.succ(sit)->s)
            sit = T.succ(sit);

        snode *sit_succ = T.succ(sit);    
        s_succ = sit_succ->s;

        do {
            s = sit->s;

            if (s->q == v) {
                L.push_back(s);

                snode *sit_pred = T.pred(sit);

                if (sit_pred->inf == s) {
                    sit_pred->inf = s->q;
                    //sit_pred->inf = sit->inf;
                }

                sit = sit_pred;
            } else {
                C.push_back(s);
                
                sit = T.pred(sit);
            }
        } while (sit->inf == v || sit->inf == T.succ(sit)->s);

        s_pred = sit->s;
    }

    if (L.size() + C.size() + U.size() > 1) {
        int intersect = _report_intersection(L, C, U);
        if (intersect)
            dcel_resolve(v, L, C, U, H);
    }

    for (segment *s : C) {
        sit = T.lookup(s);
        if (sit->inf == T.succ(sit)->s)
            C_info.push_back(sit->inf);
        else
            C_info.push_back(NULL);
        T.erase(s);
    }

    assert(C_info.size() == C.size());
    
    for (segment *s : L)
        T.erase(s);
    
    T.update_sweep_position(v->x, v->y);

    for (size_t i = 0; i < C.size(); i++) {
        segment *s = C[i];
        sit = T.insert(s);
        sit->inf = C_info[i];
    }

    for (size_t i = 0; i < U.size(); i++) {
        s = U[i];

        event *xit = Q.insert(s->q);
        xit->inf = s;
        s->q = xit->key;
        
        sit = T.insert(s);
        
        snode *sit_succ = T.succ(s);
        snode *sit_pred = T.pred(s);

        if (s_succ == NULL) {
            s_succ = sit_succ->s;
            s_pred = sit_pred->s;
        }

        assert(!segments_are_colli(sit_succ->s, s));

        if (segments_are_colli(sit_pred->s, s))
            sit_pred->inf = s;
    }


    if (s_succ) {
        snode *sit_succ = T.lookup(s_succ);
        snode *sit_pred = T.lookup(s_pred);
        
        snode *sit_first = T.succ(sit_pred);
        _find_new_event(Q, T, V, sit_pred, sit_first);

        snode *sit_last = T.pred(sit_succ);
        if (sit_last != sit_pred)
            _find_new_event(Q, T, V, sit_last, sit_succ);


        // Store the half-edge immediately to the left of v on the sweep line
        v->left = s_pred->rep;
    }
}

void sweep(std::vector<segment *> &S, std::vector<vertex *> &V,
    std::vector<hedge *> &H, queue &Q, status &T)
{
    puts("SWEEP START...");

    while (!Q.empty()) {
        event *E = Q.min();
        _handle_event(E, Q, T, V, S, H);
        Q.erase(E);
    }

    puts("SWEEP DONE.");
}

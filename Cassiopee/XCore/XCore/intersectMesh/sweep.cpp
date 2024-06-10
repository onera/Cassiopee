#include <vector>
#include <cassert>

#include "queue.h"
#include "status.h"
#include "primitives.h"
#include "dcel.h"
#include "event.h"

static
E_Int report_E_Intersection_(const std::vector<Segment *> &L,
    const std::vector<Segment *> &C, const std::vector<Segment *> &U)
{
    E_Int ncolors[2] = {0, 0};

    for (Segment *s : L) ncolors[s->color]++;
    for (Segment *s : C) ncolors[s->color]++;
    for (Segment *s : U) ncolors[s->color]++;

    // We have an E_Intersection if both colors exist
    if (ncolors[Dcel::RED] == 0 || ncolors[Dcel::BLACK] == 0) return 0;

    return 1;
}

void sweep(Queue &Q, Status &T, std::vector<Segment *> &S,
    std::vector<Vertex *> &I, std::vector<Hedge *> &H)
{
    size_t s_pos = 0;

    while (!Q.empty()) {
        Event *event = Q.min();
        Segment *seg = event->inf;
        Vertex *p = event->key;

        // Get starting segments
        std::vector<Segment *> U;

        while (s_pos < S.size() && S[s_pos]->p == p) {
            U.push_back(S[s_pos]);
            s_pos++;
            //S.pop_back();
        }

        Snode *sit = NULL;

        if (seg == NULL) {
            Segment s(p);
            sit = T.locate(&s);
            if (sit) {
                seg = sit->key;
            }
        }

        // Handle passing and ending segments
        std::vector<Segment *> L, C;
        std::vector<void *> C_info;

        Segment *seg_succ = NULL;
        Segment *seg_pred = NULL;

        if (seg != NULL) {
            sit = T.lookup(seg);
            assert(sit);

            while (sit->inf == p || sit->inf == T.succ(sit)->key) {
                sit = T.succ(sit);
            }

            Snode *sit_succ = T.succ(sit);
            seg_succ = sit_succ->key;

            do {
                seg = sit->key;

                if (seg->q == p) {
                    L.push_back(seg);

                    Snode *sit_pred = T.pred(sit);

                    if (sit_pred->inf == seg) {
                        //assert(0);
                        sit_pred->inf = seg->q;
                    }

                    sit = sit_pred;
                } else {
                    C.push_back(seg);

                    sit = T.pred(sit);
                }
            } while (sit->inf == p || sit->inf == T.succ(sit)->key);

            seg_pred = sit->key;
            assert(T.lookup(seg_pred) == sit);
        }

        // Resolve dcel E_Intersections
        if (L.size() + C.size() + U.size() > 1) {
            E_Int E_Intersect = report_E_Intersection_(L, C, U);
            if (E_Intersect) Dcel::resolve(p, L, C, U, H);
        }

        // Cache collinear segments info and delete passing segments

        for (Segment *seg : C) {
            sit = T.lookup(seg);
            if (sit->inf == T.succ(sit)->key) {
                C_info.push_back(sit->inf);
            } else {
                C_info.push_back(NULL);
            }
            T.erase(seg);
        }

        assert(C_info.size() == C.size());

        // Delete ending segments
        for (Segment *seg : L) T.erase(seg);
    
        // Advance sweep line
        T.rx = p->x, T.ry = p->y;

        // Insert passing segments
        for (size_t i = 0; i < C.size(); i++) {
            Segment *s = C[i];
            T.insert(s, C_info[i]);
        }

        // Insert starting segments and compute new E_Intersections
        for (size_t i = 0; i < U.size(); i++) {
            seg = U[i];

            // TODO: optimize
            Event *xit = Q.lookup(seg->q->x, seg->q->y);
            assert(xit);
            xit->inf = seg;

            //Event *xit = Q.insert(seg->q->x, seg->q->y, seg);
            //seg->q = xit->key;
            
            sit = T.insert(seg);

            Snode *sit_succ = T.succ(seg);
            Snode *sit_pred = T.pred(seg);

            if (seg_succ == NULL) {
                seg_succ = sit_succ->key;
                seg_pred = sit_pred->key;
            }

            assert(!seg->overlaps(*sit_succ->key));

            if (seg->overlaps(*sit_pred->key)) {
                sit_pred->inf = seg;
            }
        }

        if (seg_succ) {
            assert(seg_pred);
            Snode *sit_succ = T.lookup(seg_succ);
            Snode *sit_pred = T.lookup(seg_pred);

            assert(sit_succ);
            assert(sit_pred);
            
            Snode *sit_first = T.succ(sit_pred);
            compute_E_Intersection(Q, sit_pred, sit_first, I); 

            Snode *sit_last = T.pred(sit_succ);
            if (sit_last != sit_pred) {
                compute_E_Intersection(Q, sit_last, sit_succ, I);
            }

            // Store the half-edge immediately to the left of v on the sweep line
            p->left = seg_pred->rep;
        }

        // Done
        Q.erase(event);
    }
}

#include "dcel.h"
#include "hedge.h"
#include "vertex.h"
#include "primitives.h"
#include "smesh.h"
#include "io.h"

void Dcel::init_mh_sv_intersections(const Smesh &M)
{
    for (Vertex *v : V) {

        // Strictly S points

        if (v->oid[0] == -1 && v->oid[1] != -1) {
            const auto &loc = v->loc;
            if (loc.e_idx == -1) continue;

            // Set the reference edge for normal computation
            E_Int fid = loc.fid;
            const auto &pe = M.F2E[fid];
            E_Int e = pe[loc.e_idx];
            v->meid = e;

            // Add the vertex to hedge intersections
            // TODO(Imad): do we need to take the correctly-oriented hedge?
            Hedge *h = H[2*e];
            if (cmp_vtx(h->orig, h->twin->orig) > 0)
                h = h->twin;
            hedge_intersections[h].push_back(v);
        }
    }
}

/*
void Dcel::cut_hedges_intersecting_with_endpoint(Vertex *v, const Smesh &M)
{
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
*/

void Dcel::cut_hedges()
{
    size_t hid = 0;

    for (auto &hdata : hedge_intersections) {
        Hedge *h = hdata.first;
        auto &xs = hdata.second;

        // Distances (squared) of intersections to sh origin

        Vertex *o = h->orig;

        for (size_t i = 0; i < xs.size()-1; i++) {
            assert(cmp_vtx(xs[i], xs[i+1]) != 0);
        }

        for (Vertex *x : xs) {
            E_Float D[3] = {x->x-o->x, x->y-o->y, x->z-o->z};
            x->d = K_MATH::dot(D, D, 3);
        }

        std::sort(xs.begin(), xs.end(), [&] (const Vertex *a, const Vertex *b)
        {
            if (Sign(a->d - b->d) == 0) {
                hedge_write("hedge", h);
                point_write("a", a->x, a->y, a->z);
                point_write("b", b->x, b->y, b->z);
                assert(Sign(a->d - b->d) != 0);
            }
            return a->d < b->d;
        });

        Hedge *current_h = h;
        Hedge *current_t = h->twin;

        // Cache original connectivity
        Vertex *tail = current_t->orig;
        //Hedge *h_next = h->next;
        //Hedge *h_prev = h->prev;

        assert(current_h->color == current_t->color);

        for (Vertex *x : xs) {
            // Create two new hedges with x as their origin
            Hedge *e1 = new Hedge(x);
            Hedge *e2 = new Hedge(x);

            // Add to hedges collection
            H.push_back(e1);
            H.push_back(e2);
            
            // Copy the color
            e1->color = current_h->color;
            e2->color = current_t->color;

            // Copy the face record
            e1->left = current_h->left;
            e2->left = current_t->left;

            // Pair-up the new half-edges
            current_h->twin = e2;
            e2->twin = current_h;
            current_t->twin = e1;
            e1->twin = current_t;

            // Set prev and next pointers at the endpoints
            //e1->next = current_h->next;
            //e2->next = current_t->next;
            //current_h->next->prev = e1;
            //current_t->next->prev = e2;

            // Set prev and next pointers at x
            current_h->next = e1;
            e1->prev = current_h;
            current_t->next = e2;
            e2->prev = current_t;

            // Keep going
            current_h = e1;
        }

        assert(current_h->twin->orig == tail);

        hid++;
    }
}


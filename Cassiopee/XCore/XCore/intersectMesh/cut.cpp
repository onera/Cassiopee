#include "dcel.h"
#include "hedge.h"
#include "vertex.h"
#include "primitives.h"

void Dcel::cut_hedges(std::map<Hedge *, std::vector<Vertex *>> &hmap)
{
    size_t hid = 0;

    for (auto &hdata : hmap) {
        Hedge *h = hdata.first;
        auto &xs = hdata.second;

        // Distances (squared) of intersections to sh origin

        Vertex *o = h->orig;

        E_Float D[3];

        for (Vertex *x : xs) {
            E_Float D[3] = {x->x-o->x, x->y-o->y, x->z-o->z};
            x->d = K_MATH::dot(D, D, 3);
        }

        std::sort(xs.begin(), xs.end(), [&] (const Vertex *a, const Vertex *b)
        {
            assert(Sign(a->d - b->d) != 0);
            return a->d < b->d;
        });

        Hedge *current_h = h;
        Hedge *current_t = h->twin;

        // Cache the tail
        Vertex *tail = current_t->orig;

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

            // Set prev and next pointers at x
            current_h->next = e1;
            e1->prev = current_h;
            current_t->next = e2;
            e2->prev = current_t;

            // NOTE(Imad): is this part necessary?
            // Set prev and next pointers at the endpoints
            e1->next = current_h->next;
            e2->next = current_t->next;
            current_h->next->prev = e1;
            current_t->next->prev = e2;

            // Keep going
            current_h = e1;
        }

        assert(current_h->twin->orig == tail);

        hid++;
    }
}


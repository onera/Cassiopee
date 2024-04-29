#include "proto.h"
#include <algorithm>
#include <cassert>

hedge::hedge(vertex *v)
: orig(v), twin(NULL), next(NULL), prev(NULL), left(NULL), color(RED),
  cycl(NULL), s(0)
{}

static
int _partition(std::vector<hedge *> &H, int low, int high)
{
    hedge *pivot = H[high];
    int i = low-1;

    for (int j = low; j < high; j++) {
        if (hedge_cmp_cwise(H[j], pivot) <= 0) {
            i++;
            std::swap(H[i], H[j]);
        }
    }

    i++;
    std::swap(H[i], H[high]);
    return i;
}

void hedge_sort_cwise(std::vector<hedge *> &H, int low, int high)
{
    if (low >= high)
        return;

    int p = _partition(H, low, high);

    hedge_sort_cwise(H, low, p-1);
    hedge_sort_cwise(H, p+1, high);
}

void hedge_sort_ccwise(std::vector<hedge *> &H, int low, int high)
{
    hedge_sort_cwise(H, low, high);
    std::reverse(H.begin(), H.end());
}

int hedge_cmp_cwise(hedge *h, hedge *w)
{
    assert(h->orig == w->orig);
    vertex *c = h->orig;
    vertex *a = h->twin->orig;
    vertex *b = w->twin->orig;

    double ax = a->x;
    double ay = a->y;
    double bx = b->x;
    double by = b->y;
    double cx = c->x;
    double cy = c->y;

    double acx = ax - cx;
    double acy = ay - cy;
    double bcx = bx - cx;
    double bcy = by - cy;

    int sign_acx = sign(acx);
    int sign_acy = sign(acy);
    int sign_bcx = sign(bcx);
    int sign_bcy = sign(bcy);

    if (sign_acx >= 0 && sign_bcx < 0)
        return -1;
    if (sign_acx < 0 && sign_bcx >= 0)
        return 1;
    if (sign_acx == 0 && sign_bcx == 0) {
        if (sign_acy >= 0 || sign_bcy >= 0)
            return sign(by - ay);
        return sign(ay - by);
    }
    
    double det = acx*bcy - bcx*acy;
    int cmp = sign(det);

    if (cmp < 0)
        return -1;
    else if (cmp > 0)
        return 1;

    // Overlapping segments
    
    assert(h->color != w->color);

    // If right half, red before black
    // Otherwise, black before red

    cmp = sign(h->color - w->color);

    if (sign(acx) >= 0) {
        assert(sign(bcx) >= 0);
        return cmp;
    } else {
        return -cmp;
    }
}

int hedges_are_collinear(hedge *h, hedge *w)
{
    vertex *hp = h->orig;
    vertex *hq = h->twin->orig;
    vertex *wp = w->orig;
    vertex *wq = w->twin->orig;

    return vertex_orient(hp, hq, wp) == 0 &&
           vertex_orient(hp, hq, wq) == 0;
}

int hedges_incident_overlap(hedge *h, hedge *w)
{
    assert(h->orig == w->orig);
    vertex *a = h->orig;
    vertex *b = h->twin->orig;
    vertex *c = w->twin->orig;

    if (vertex_orient(a, b, c) != 0)
        return 0;

    // Dot product should be positive.
    double dp = (b->x-a->x)*(c->x-a->x) + (b->y-a->y)*(c->y-a->y);
    
    return sign(dp) >= 0;
}

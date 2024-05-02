#include "proto.h"
#include <cassert>
#include <cstddef>
#include <cstdio>

vertex::vertex(double X, double Y, int Oid, bool is_new)
{
    x = X;
    y = Y;
    if (!is_new) {
        oid[0] = Oid;
        nid = -1;
    } else {
        oid[0] = -1;
        nid = Oid;
    }
    oid[1] = -1;
    inc = NULL;
    left = NULL;
}

void vertex::print()
{
    printf(SF_D_ " ", nid);
    if (nid != -1) printf(SF_D_ " ", nid);
    else if (oid[0] != -1 && oid[1] != -1) printf(SF_D_ "/" SF_D_ " ", oid[0], oid[1]);
    else if (oid[0] != -1) printf(SF_D_ " ", oid[0]);
    else printf(SF_D_ " ", oid[1]);
}

void vertex::fprint(FILE *fh)
{
    if (nid != -1) fprintf(fh, SF_D_ " ", nid);
    else fprintf(fh, SF_D_ " ", oid[0]);
}

edge::edge(vertex *P, vertex *Q)
: p(P), q(Q), color(RED)
{}

int vertex_cmp_xy(vertex *p, vertex *q)
{
    return xy_cmp(p->x, p->y, q->x, q->y);
}

int vertex_cmp(vertex *p, vertex *q)
{
    if (p == q)
        return 0;

    int cmp = xy_cmp(p->x, p->y, q->x, q->y);

    assert(cmp);

    return cmp;
}

// Position of c wrt directed segment ab
int vertex_orient(vertex *a, vertex *b, vertex *c)
{
    double det = (b->x - a->x)*(c->y - a->y) - (b->y - a->y)*(c->x - a->x);
    return sign(det);
}

int vertex_cmp_cwise(vertex *a, vertex *b, vertex *c)
{
    if (a == b)
        return 0;

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

    double da = acx*acx + acy*acy;
    double db = bcx*bcx + bcy*bcy;
    cmp = sign(da - db);
    assert(cmp);
    return cmp;
}

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
#include <algorithm>
#include <cassert>

#include "hedge.h"
#include "vertex.h"
#include "primitives.h"
#include "dcel.h"

Hedge::Hedge(Vertex *v)
: orig(v), twin(NULL), prev(NULL), next(NULL), left(NULL), color(Dcel::NO_IDEA),
  cycle(NULL)
{}

static
Int _partition(std::vector<Hedge *> &H, Int low, Int high)
{
    Hedge *pivot = H[high];
    Int i = low-1;

    for (Int j = low; j < high; j++) {
        if (Hedge::cmp_cwise(H[j], pivot) <= 0) {
            i++;
            std::swap(H[i], H[j]);
        }
    }

    i++;
    std::swap(H[i], H[high]);
    return i;
}

void Hedge::sort_cwise(std::vector<Hedge *> &H, Int low, Int high)
{
    if (low >= high)
        return;

    Int p = _partition(H, low, high);

    sort_cwise(H, low, p - 1);
    sort_cwise(H, p + 1, high);
}

void Hedge::sort_ccwise(std::vector<Hedge *> &H, Int low, Int high)
{
    sort_cwise(H, low, high);
    std::reverse(H.begin(), H.end());
}

Int Hedge::cmp_cwise(const Hedge *h, const Hedge *w)
{
    assert(h != w);
    assert(h->orig == w->orig);
    Vertex *c = h->orig;
    Vertex *a = h->twin->orig;
    Vertex *b = w->twin->orig;

    // a is tail of h
    // b is tail of w
    // c is orig

    Float ax = h->proj_tx;
    Float ay = h->proj_ty;
    Float bx = w->proj_tx;
    Float by = w->proj_ty;
    Float cx = h->proj_ox;
    Float cy = h->proj_oy;

    assert(Sign(cx-w->proj_ox) == 0);
    assert(Sign(cy-w->proj_oy) == 0);

    long double acx = (long double)ax - (long double)cx;
    long double acy = (long double)ay - (long double)cy;
    long double bcx = (long double)bx - (long double)cx;
    long double bcy = (long double)by - (long double)cy;

    Int sign_acx = Sign(acx);
    Int sign_acy = Sign(acy);
    Int sign_bcx = Sign(bcx);
    Int sign_bcy = Sign(bcy);

    if (sign_acx >= 0 && sign_bcx < 0)
        return -1;
    if (sign_acx < 0 && sign_bcx >= 0)
        return 1;
    if (sign_acx == 0 && sign_bcx == 0) {
        if (sign_acy >= 0 || sign_bcy >= 0) {
            long double diff = (long double)ay - (long double)by;
            if (Sign(diff) > 0) return -1;
            else return 1;
        }

        long double diff = (long double)by - (long double)ay;
        if (Sign(diff) > 0) return -1;
        else return 1;
    }
    
    Float det = acx * bcy - bcx * acy;
    //Float det = DifferenceOfProducts(acx, bcy, bcx, acy);
    Int cmp = Sign(det);

    if (cmp < 0)
        return -1;
    else if (cmp > 0)
        return 1;

    // Overlapping segments
    /*
    if (h->color == w->color) {
        assert(h->orig == w->orig);
        FILE *fh = fopen("bad_h", "w");
        fprintf(fh, "POINTS\n");
        fprintf(fh, "2\n");
        fprintf(fh, "%f %f %f\n", h->orig->x, h->orig->y, h->orig->z);
        Vertex *tail = h->twin->orig;
        fprintf(fh, "%f %f %f\n", tail->x, tail->y, tail->z);
        fprintf(fh, "EDGES\n");
        fprintf(fh, "1\n");
        fprintf(fh, "0 1\n");
        fclose(fh);


        fh = fopen("bad_w", "w");
        fprintf(fh, "POINTS\n");
        fprintf(fh, "2\n");
        fprintf(fh, "%f %f %f\n", w->orig->x, w->orig->y, w->orig->z);
        tail = w->twin->orig;
        fprintf(fh, "%f %f %f\n", tail->x, tail->y, tail->z);
        fprintf(fh, "EDGES\n");
        fprintf(fh, "1\n");
        fprintf(fh, "0 1\n");
        fclose(fh);
    }
    */
    assert(h->color != w->color);

    // If right half, red before black
    // Otherwise, black before red

    cmp = Sign(h->color - w->color);

    if (sign_acx >= 0) {
        assert(sign_bcx >= 0);
        return cmp;
    } else {
        return -cmp;
    }
}

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

#include "primitives.h"
#include "event.h"

Float TOL = 1e-8;

Int Sign(Float x)
{
    if (x > TOL) return 1;
    if (x < -TOL) return -1;
    //assert(x == 0.0);
    return 0;
}

Int cmp_points(Float x1, Float y1, Float z1, Float x2, Float y2, Float z2)
{
    Float t = x1 - x2;
    Int s = Sign(t);
    if (s) return s;
    
    t = y1 - y2;
    s = Sign(t);
    if (s) return s;

    t = z1 - z2;
    return Sign(t);
}

/*
Int cmp_points(Float x1, Float y1, Float x2, Float y2)
{
    Float t = x1 - x2;
    Int s = Sign(t);
    if (s) return s;
    
    t = y1 - y2;
    return Sign(t);
}
*/

Int cmp_segments
(
    Float px0, Float py0,
    Float px1, Float py1,
    Float qx1, Float qy1,
    Float rx, Float ry,
    Float dx0, Float dy0,
    Float dx1, Float dy1
)
{
    Float T1 = dy0 * dx1 - dy1 * dx0;

    Int sign1 = Sign(T1);

    if (sign1 == 0) {

        Float mdx = qx1 - px0;
        Float mdy = qy1 - py0;

        Int sign2 = Sign(dy0 * mdx - mdy * dx0);
        
        if (sign2 == 0) {
            Int sign3 = Sign(dy1 * mdx - mdy * dx1);

            assert(sign3 == 0);

            if (sign3 == 0) {
                return 0;
            }
        }
    }

    if (Sign(dx0) == 0) {
        Float T2 = (py1 * dx1 - px1 * dy1) + (dy1 * rx - ry * dx1);
        Int sign2 = Sign(T2);
        return (sign2 <= 0) ? 1 : -1;
    }

    if (Sign(dx1) == 0) {
        Float T2 = (py0 * dx0 - px0 * dy0) + (dy0 * rx - ry * dx0);
        Int sign2 = Sign(T2);
        return (sign2 <= 0) ? -1 : 1;
    }

    Float T2 = dx1 * (py0 * dx0 + dy0 * (rx - px0)) - dx0
        * (py1 * dx1 + dy1 * (rx - px1));

    Int sign2 = Sign(T2);
    if (sign2 != 0) return sign2;

    Float T3 = (py0 * dx0 - px0 * dy0) + (dy0 * rx - ry * dx0);
    Int sign3 = Sign(T3);
    return (sign3 <= 0) ? sign1 : -sign1;
}

// We define the compare functions for Vertexs and Segments by first calling
// cmp_points and cmp_segments on the floating point filter coordinates of the
// corresponding points and segments. In the case that these calls do not
// return a reliable result (i.e. return NO_IDEA) we call them again with the
// exact routines.
Int compare(const Vertex &a, const Vertex &b)
{
    return cmp_points(a.x, a.y, a.z, b.x, b.y, b.z);
}

Int compare(const Segment &s1, const Segment &s2, Float rx, Float ry)
{
    return cmp_segments(s1.p->x, s1.p->y,
                        s2.p->x, s2.p->y,
                        s2.q->x, s2.q->y,
                        rx, ry,
                        s1.dx, s1.dy, s2.dx, s2.dy);
}

Int cmp_mySeg(const Segment &s1, const Segment &s2)
{
    Int cmp = cmp_points(s1.p->x, s1.p->y, s1.p->z, s2.p->x, s2.p->y,
        s2.p->z);
    if (cmp) return cmp;

    cmp = Sign(s1.color - s2.color);
    if (cmp) return cmp;

    cmp = Sign(s1.id - s2.id);

    assert(cmp);

    return cmp;
}

void compute_intersection(Queue &Q, Snode *sit0, Snode *sit1,
    std::vector<Vertex *> &I)
{
    Segment s0 = *sit0->key;
    Segment s1 = *sit1->key;

    Float w = s0.dy * s1.dx - s1.dy * s0.dx;
    Int i = Sign(w);
    if (i == -1 || i == 0) return;

    Float c1 = s0.X2() * s0.Y1() - s0.X1() * s0.Y2();
    Float c2 = s1.X2() * s1.Y1() - s1.X1() * s1.Y2();

    Float x = c2 * s0.dx - c1 * s1.dx;
    Float d0 = x - s0.X2() * w;
    if (Sign(d0) > 0) return;
    if (Sign(x - s1.X2() * w) > 0) return;

    Float y = c2 * s0.dy - c1 * s1.dy;
    if (Sign(d0) == 0 &&
        Sign(y - s0.Y2() * w) > 0) return;


    x /= w;
    y /= w;

    Event *xit = Q.lookup(x, y, 0);

    if (xit == NULL) {
        xit = Q.insert(x, y, 0);
        xit->key->id = I.size();
        I.push_back(xit->key);
    }
    
    xit->inf = sit0->key;
    sit0->inf = xit->key;
}

Float DifferenceOfProducts(Float a, Float b, Float c, Float d)
{
    Float cd = c * d;
    Float differenceOfProducts = std::fma(a, b, -cd);
    Float err = std::fma(-c, d, cd);
    return differenceOfProducts + err;
}

Float dRand(Float dMin, Float dMax)
{
    Float d = (Float) rand() / RAND_MAX;
    return dMin + d * (dMax - dMin);
}

Int is_point_on_segment(Float px, Float py, Float pz, Float ax, Float ay,
    Float az, Float bx, Float by, Float bz)
{
    Float Vab[3] = {bx-ax, by-ay, bz-az};
    Float Vap[3] = {px-ax, py-ay, pz-az};
    Float N[3];
    K_MATH::cross(Vab, Vap, N);
    if (Sign(K_MATH::norm(N, 3)) != 0) return 0;

    Float Vbp[3] = {px-bx, py-by, pz-bz};

    Float dp = K_MATH::dot(Vap, Vab, 3);
    if (dp < -TOL) return 0;

    dp = K_MATH::dot(Vbp, Vab, 3);
    if (dp > TOL) return 0;

    return 1;
}

Int EdgeEdgeIntersect(Float ax, Float ay, Float az, Float bx, Float by,
    Float bz, Float px, Float py, Float pz, Float qx, Float qy, Float qz,
    Float &ix, Float &iy, Float &iz)
{
    Float d1[3] = {bx-ax, by-ay, bz-az};
    Float d2[3] = {qx-px, qy-py, qz-pz};
    Float r[3] = {px-ax, py-ay, pz-az};

    Float d1d2[3];
    K_MATH::cross(d1, d2, d1d2);
    Float denom = K_MATH::dot(d1d2, d1d2, 3);

    if (Sign(denom) == 0) {

        Float colli[3];
        K_MATH::cross(d1, r, colli);
        Float NORM = K_MATH::norm(colli, 3);
        if (Sign(NORM) == 0) {
            assert("collinear!" && 0);
        } else {
            return 0;
        }
    }

    Float tmp[3];
    K_MATH::cross(r, d2, tmp);
    Float t = K_MATH::dot(tmp, d1d2, 3);
    t /= denom;
    if (t < -TOL) return 0;

    K_MATH::cross(r, d1, tmp);
    Float u = K_MATH::dot(tmp, d1d2, 3);
    u /= denom;
    if (u < -TOL || u > 1 + TOL) return 0;

    ix = px + u*(qx - px);
    iy = py + u*(qy - py);
    iz = pz + u*(qz - pz);

    return 1;
}



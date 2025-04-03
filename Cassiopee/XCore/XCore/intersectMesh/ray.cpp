/*    
    Copyright 2013-2025 Onera.

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
#include "ray.h"
#include "primitives.h"
#include "AABB.h"
#include "point.h"

Ray::Ray(E_Float px, E_Float py, E_Float pz, E_Float dx, E_Float dy, E_Float dz,
    Policy p)
{
    o[0] = px, o[1] = py, o[2] = pz;
    d[0] = dx, d[1] = dy, d[2] = dz;
    inv_dx = (d[0] != 0.0) ? 1.0/d[0] : EFLOATMAX;
    inv_dy = (d[1] != 0.0) ? 1.0/d[1] : EFLOATMAX;
    inv_dz = (d[2] != 0.0) ? 1.0/d[2] : EFLOATMAX;

    E_Float adx = fabs(d[0]);
    E_Float ady = fabs(d[1]);
    E_Float adz = fabs(d[2]);

    if (adx > ady && adx > adz) kz = 0;
    else if (ady > adz) kz = 1;
    else kz = 2;

    kx = kz+1; if (kx == 3) kx = 0;
    ky = kx+1; if (ky == 3) ky = 0;

    if (d[kz] < 0.0) std::swap(kx, ky);

    Sx = d[kx] / d[kz];
    Sy = d[ky] / d[kz];
    Sz = 1.0 / d[kz];

    policy = p;
}

E_Float diff_of_prod(E_Float a, E_Float b, E_Float c, E_Float d)
{
    E_Float cd = c*d;
    E_Float dp = std::fma(a, b, -cd);
    E_Float err = std::fma(-c, d, cd);
    return dp + err;
}

bool Ray::intersect_triangle(const Point &tA, const Point &tB,
    const Point &tC, E_Float &u, E_Float &v, E_Float &w, E_Float &t,
    E_Float &x, E_Float &y, E_Float &z) const
{
    const Point O = (Point){o[0], o[1], o[2]};
    const Point A = tA - O;
    const Point B = tB - O;
    const Point C = tC - O;

    const E_Float Ax = A[kx] - Sx*A[kz];
    const E_Float Ay = A[ky] - Sy*A[kz];
    const E_Float Bx = B[kx] - Sx*B[kz];
    const E_Float By = B[ky] - Sy*B[kz];
    const E_Float Cx = C[kx] - Sx*C[kz];
    const E_Float Cy = C[ky] - Sy*C[kz];

    u = diff_of_prod(Cx, By, Cy, Bx);
    v = diff_of_prod(Ax, Cy, Ay, Cx);
    w = diff_of_prod(Bx, Ay, By, Ax);

    //u = Cx*By - Cy*Bx;
    //v = Ax*Cy - Ay*Cx;
    //w = Bx*Ay - By*Ax;

    // WARNING: THIS IS NOT ROBUST!
    // Forward-error analysis + arbitray precision needed.
    if ((u < -TOL || v < -TOL || w < -TOL) && (u > TOL || v > TOL || w > TOL))
        return false;
    
    E_Float det = u + v + w;

    if (det == 0.0) return false;

    const E_Float Az = Sz*A[kz];
    const E_Float Bz = Sz*B[kz];
    const E_Float Cz = Sz*C[kz];
    t = u*Az + v*Bz + w*Cz;

    if (policy == FORWARD) {
        if (det < 0.0 && t >= 0.0) return false;
        else if (det > 0.0 && t <= 0.0) return false;
    }

    const E_Float inv_det = 1.0 / det;
    u = u * inv_det;
    v = v * inv_det;
    w = w * inv_det;
    t = t * inv_det;

    x = u * tA.x + v * tB.x + w * tC.x;
    y = u * tA.y + v * tB.y + w * tC.y;
    z = u * tA.z + v * tB.z + w * tC.z;

    return true;
}

bool Ray::intersect_AABB(const AABB &box) const
{
    E_Float tmin = 0;
    E_Float tmax = EFLOATMAX;

    E_Float boxMin[3] = { box.xmin, box.ymin, box.zmin };
    E_Float boxMax[3] = { box.xmax, box.ymax, box.zmax };

    for (int i = 0; i < 3; i++) {
        E_Float O = o[i];
        E_Float D = d[i];
        E_Float bmin = boxMin[i];
        E_Float bmax = boxMax[i];

        if (d != 0) {
            E_Float t1 = (bmin - O) / D;
            E_Float t2 = (bmax - O) / D;

            if (t1 > t2) { std::swap(t1, t2); }

            tmin = (t1 > tmin) ? t1 : tmin;
            tmax = (t2 < tmax) ? t2 : tmax;

            if (tmin > tmax) return false;  // No intersection
        } else {
            if (O < bmin || O > bmax) return false;  // Parallel and outside slab
        }
    }

    return true;
}

E_Int ray_point_orient(const E_Float o[3], const E_Float d[3],
    const E_Float fN[3], E_Float px, E_Float py, E_Float pz)
{
    E_Float w[3] = {px-o[0], py-o[1], pz-o[2]};
    E_Float c[3];
    K_MATH::cross(d, w, c);
    E_Float dp = K_MATH::dot(c, fN, 3);
    E_Int cmp = Sign(dp);
    if (cmp > 0) return 1;
    if (cmp < 0) return -1;
    return 0;
}

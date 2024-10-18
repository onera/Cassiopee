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
#include "ray.h"
#include "primitives.h"
#include "AABB.h"

E_Int ray_point_orient(const E_Float o[3], const E_Float d[3],
    const E_Float fN[3], E_Float px, E_Float py, E_Float pz)
{
    E_Float w[3] = {px-o[0], py-o[1], pz-o[2]};
    E_Float c[3];
    K_MATH::cross(d, w, c);
    E_Float dp = K_MATH::dot(c, fN, 3);
    // TODO(Imad): needs FEA + float128
    E_Int cmp = Sign(dp);
    if (cmp > 0) return 1;
    if (cmp < 0) return -1;
    return 0;
}

bool ray_AABB_intersect(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz,
    const AABB &box)
{
    E_Float inv_dx = (dx != 0.0) ? 1.0/dx : EFLOATMAX;
    E_Float inv_dy = (dy != 0.0) ? 1.0/dy : EFLOATMAX;
    E_Float inv_dz = (dz != 0.0) ? 1.0/dz : EFLOATMAX;

    E_Float t_min = EFLOATMIN;
    E_Float t_max = EFLOATMAX;

    // X
    if (dx == 0.0) {
        if (ox < box.xmin || ox > box.xmax) {
            return false;
        }
    } else {
        E_Float t1 = (box.xmin - ox) * inv_dx;
        E_Float t2 = (box.xmax - ox) * inv_dx;

        if (t1 > t2) std::swap(t1, t2);

        t_min = std::max(t_min, t1);
        t_max = std::max(t_max, t2);

        if (t_min > t_max)
            return false;
    }

    // Y
    if (dy == 0.0) {
        if (oy < box.ymin || oy > box.ymax) {
            return false;
        }
    } else {
        E_Float t1 = (box.ymin - oy) * inv_dy;
        E_Float t2 = (box.ymax - oy) * inv_dy;

        if (t1 > t2) std::swap(t1, t2);

        t_min = std::max(t_min, t1);
        t_max = std::max(t_max, t2);

        if (t_min > t_max)
            return false;
    }

    // Z
    if (dz == 0.0) {
        if (oz < box.zmin || oz > box.zmax) {
            return false;
        }
    } else {
        E_Float t1 = (box.zmin - oz) * inv_dz;
        E_Float t2 = (box.zmax - oz) * inv_dz;

        if (t1 > t2) std::swap(t1, t2);

        t_min = std::max(t_min, t1);
        t_max = std::max(t_max, t2);

        if (t_min > t_max)
            return false;
    }

    return true;
}

bool MollerTrumboreAnyDir(
    E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz,
    E_Float ax, E_Float ay, E_Float az,
    E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz,
    E_Float &u, E_Float &v, E_Float &w, E_Float &t,
    E_Float &x, E_Float &y, E_Float &z)
{
    E_Float v1[3] = {bx-ax, by-ay, bz-az};
    E_Float v2[3] = {cx-ax, cy-ay, cz-az};

    E_Float d[3] = {dx, dy, dz};
    E_Float h[3];
    K_MATH::cross(d, v2, h);
    E_Float det = K_MATH::dot(v1, h, 3);
    if (Sign(det) == 0) return false;

    E_Float inv_det = 1.0 / det;

    E_Float s[3] = {ox-ax, oy-ay, oz-az};
    v = K_MATH::dot(s, h, 3) * inv_det;
    if (v < -TOL || v > 1+TOL) return false;

    E_Float q[3];
    K_MATH::cross(s, v1, q);
    w = K_MATH::dot(d, q, 3) * inv_det;
    if (w < -TOL || w > 1+TOL) return false;

    u = 1 - v - w;
    if (u < -TOL || u > 1+TOL) return false;

    t = K_MATH::dot(v2, q, 3) * inv_det;
    x = ox + t * dx;
    y = oy + t * dy;
    z = oz + t * dz;

    return true;
}

E_Int MollerTrumbore(E_Float px, E_Float py, E_Float pz, E_Float dx, E_Float dy,
    E_Float dz, E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz, TriangleIntersection &TI)
{
    E_Float e1x = bx - ax;
    E_Float e1y = by - ay;
    E_Float e1z = bz - az;

    E_Float e2x = cx - ax;
    E_Float e2y = cy - ay;
    E_Float e2z = cz - az;

    E_Float pvecx = dy * e2z - dz * e2y;
    E_Float pvecy = dz * e2x - dx * e2z;
    E_Float pvecz = dx * e2y - dy * e2x;

    E_Float det = e1x * pvecx + e1y * pvecy + e1z * pvecz;

    if (det > -TOL && det < TOL) return 0;

    E_Float inv_det = 1.0 / det;

    E_Float tvecx = px - ax;
    E_Float tvecy = py - ay;
    E_Float tvecz = pz - az;

    TI.u = inv_det * (tvecx * pvecx + tvecy * pvecy + tvecz * pvecz);

    if (TI.u < -TOL || TI.u > 1+TOL) return 0;

    E_Float qvecx = tvecy * e1z - tvecz * e1y;
    E_Float qvecy = tvecz * e1x - tvecx * e1z;
    E_Float qvecz = tvecx * e1y - tvecy * e1x;

    TI.v = inv_det * (dx * qvecx + dy * qvecy + dz * qvecz);

    if (TI.v < -TOL || TI.v > 1+TOL) return 0;

    TI.w = 1 - TI.u - TI.v;

    if (TI.w < -TOL || TI.w  > 1+TOL) return 0;

    TI.t = inv_det * (e2x * qvecx + e2y * qvecy + e2z * qvecz);

    if (TI.t >= 0) {
        TI.x = px + TI.t * dx;
        TI.y = py + TI.t * dy;
        TI.z = pz + TI.t * dz;

        return 1;
    }
    
    return 0;
}

/*
Ray::Ray(Point O, E_Float D[3])
: org(O), dir(D)
{
    org.x = O.x;
    org.y = O.y;
    org.z = O.z;


    E_Float adx = fabs(dir.x);
    E_Float ady = fabs(dir.y);
    E_Float adz = fabs(dir.z);

    if (adx > ady && adx > adz) kz = 0;
    else if (ady > adz) kz = 1;
    else kz = 2;

    kx = kz+1; if (kx == 3) kx = 0;
    ky = kx+1; if (ky == 3) ky = 0;

    if (dir[kz] < 0.0) std::swap(kx, ky);

    Sx = dir[kx] / dir[kz];
    Sy = dir[ky] / dir[kz];
    Sz = 1.0 / dir[kz];
}

E_Int Ray::intersect_triangle(const Point &a, const Point &b, const Point &c,
    TriangleIntersection &TI)
{
    // Calculate vertices relative to ray origin
    const Point A = a - org;
    const Point B = b - org;
    const Point C = c - org;

    // Perform shear and scale of vertices
    const E_Float Ax = A[kx] - Sx*A[kz];
    const E_Float Ay = A[ky] - Sy*A[kz];
    const E_Float Bx = B[kx] - Sx*B[kz];
    const E_Float By = B[ky] - Sy*B[kz];
    const E_Float Cx = C[kx] - Sx*C[kz];
    const E_Float Cy = C[ky] - Sy*C[kz];

    // Calculate scale barycentric coordinates
    E_Float U = Cx*By - Cy*Bx;
    E_Float V = Ax*Cy - Ay*Cx;
    E_Float W = Bx*Ay - By*Ax;

    // Fall back to test against edges using quad precision
    if (U == 0.0 || V == 0.0 || W == 0.0) {
        __float128 CxBy = (__float128)Cx * (__float128)By;
        __float128 CyBx = (__float128)Cy * (__float128)Bx;
        U = (E_Float)(CxBy - CyBx);
        
        __float128 AxCy = (__float128)Ax * (__float128)Cy;
        __float128 AyCx = (__float128)Ay * (__float128)Cx;
        V = (E_Float)(AxCy - AyCx);
        
        __float128 BxAy = (__float128)Bx * (__float128)Ay;
        __float128 ByAx = (__float128)By * (__float128)Ax;
        W = (E_Float)(BxAy - ByAx);
    }

    // No backface culling
    if ((U < 0.0 || V < 0.0 || W < 0.0) &&
        (U > 0.0 || V > 0.0 || W > 0.0))
        return 0;
    
    // Calculate determinant
    E_Float det = U + V + W;

    if (det == 0.0) return 0;

    // Calculate scaled z-coordinates of vertices
    // and use them to calculate the hit distance.

    const E_Float Az = Sz * A[kz];
    const E_Float Bz = Sz * B[kz];
    const E_Float Cz = Sz * C[kz];
    const E_Float T = U*Az + V*Bz + W*Cz;

    if (det < 0.0 && T >= 0) return 0;
    else if (det > 0.0 && T <= 0) return 0;

    // Normalize U, V, W and T
    const E_Float inv_det = 1.0 / det;

    TI.u = U * inv_det;
    TI.v = V * inv_det;
    TI.w = W * inv_det;
    TI.t = T * inv_det;

    return 1;
}
*/

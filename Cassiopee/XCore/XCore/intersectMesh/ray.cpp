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

Ray::Ray(Point O, Vec3 D)
: org(O), dir(D)
{
    E_Float adx = fabs(dir[0]);
    E_Float ady = fabs(dir[1]);
    E_Float adz = fabs(dir[2]);

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

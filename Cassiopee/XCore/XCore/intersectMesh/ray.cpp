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

Int MollerTrumbore(Float px, Float py, Float pz, Float dx, Float dy,
    Float dz, Float ax, Float ay, Float az, Float bx, Float by, Float bz,
    Float cx, Float cy, Float cz, TriangleIntersection &TI)
{
    Float e1x = bx - ax;
    Float e1y = by - ay;
    Float e1z = bz - az;

    Float e2x = cx - ax;
    Float e2y = cy - ay;
    Float e2z = cz - az;

    Float pvecx = dy * e2z - dz * e2y;
    Float pvecy = dz * e2x - dx * e2z;
    Float pvecz = dx * e2y - dy * e2x;

    Float det = e1x * pvecx + e1y * pvecy + e1z * pvecz;

    if (det > -TOL && det < TOL)
        return 0;

    Float inv_det = 1.0 / det;

    Float tvecx = px - ax;
    Float tvecy = py - ay;
    Float tvecz = pz - az;

    TI.u = inv_det * (tvecx * pvecx + tvecy * pvecy + tvecz * pvecz);

    if (TI.u < 0 || TI.u > 1)
        return 0;

    Float qvecx = tvecy * e1z - tvecz * e1y;
    Float qvecy = tvecz * e1x - tvecx * e1z;
    Float qvecz = tvecx * e1y - tvecy * e1x;

    TI.v = inv_det * (dx * qvecx + dy * qvecy + dz * qvecz);

    if (TI.v < 0 || TI.u + TI.v > 1)
        return 0;
    
    TI.t = inv_det * (e2x * qvecx + e2y * qvecy + e2z * qvecz);

    if (TI.t > TOL) {
        TI.x = px + TI.t * dx;
        TI.y = py + TI.t * dy;
        TI.z = pz + TI.t * dz;
        return 1;
    }
    
    return 0;
}